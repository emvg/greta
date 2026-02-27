import argparse
import mudata as mu
import pandas as pd
import time
import scipy.sparse as sp
import os

from LingerGRN.pseudo_bulk import *
from LingerGRN.preprocess import *  
import LingerGRN.LINGER_tr as LINGER_tr
import LingerGRN.LL_net as LL_net

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--linger_GRN', required=True)
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-v','--version', required=True)
args = vars(parser.parse_args())

GRNdir = args['linger_GRN'] + "/"
out_dir = args['out_dir'] + "/"
path_mdata = args['path_mdata']
version = args['version']

# Read mdata
mdata = mu.read(path_mdata)
adata_RNA = mdata["rna"].copy()
adata_ATAC = mdata["atac"].copy()
log(f"RNA adata loaded:\n{adata_RNA}\n")
log(f"ATAC adata loaded:\n{adata_ATAC}\n")

# Add barcode col
adata_RNA.obs['barcode'] = adata_RNA.obs_names
adata_ATAC.obs['barcode'] = adata_ATAC.obs_names

# Add gene_ids col
adata_RNA.var['gene_ids'] = adata_RNA.var_names
adata_ATAC.var['gene_ids'] = adata_ATAC.var_names.str.replace('-', ':', n=1)

# Extract cell type annotation from mdata
label = pd.DataFrame({'barcode_use': mdata.obs['celltype'].index, 'label': mdata.obs['celltype'].values})

# Get the input data for LINGER
matrix = sp.vstack([adata_RNA.X.T, adata_ATAC.X.T])
features = pd.DataFrame(adata_RNA.var['gene_ids'].values.tolist() + adata_ATAC.var['gene_ids'].values.tolist(), columns=[1])
K = adata_RNA.shape[1]
N = K + adata_ATAC.shape[1]
types = ['Gene Expression' if i <= K-1 else 'Peaks' for i in range(0, N)]
features[2] = types
barcodes = pd.DataFrame(adata_RNA.obs['barcode'].values,columns=[0])

adata_RNA, adata_ATAC = get_adata(matrix, features, barcodes, label) 
log(f"adata_RNA for Linger:\n{adata_RNA}\n")
log(f"adata_ATAC for Linger:\n{adata_ATAC}\n")


# 3) Generate pseudo-bulk
log("Generating pseudo-bulk / metacells...")
samplelist = list(set(adata_ATAC.obs['sample'].values))
TG_pseudobulk = pd.DataFrame([])        # TG x metacells
RE_pseudobulk = pd.DataFrame([])        # RE x metacells
singlepseudobulk = adata_RNA.obs['sample'].nunique() > 10

for tempsample in samplelist:
    adata_RNAtemp = adata_RNA[adata_RNA.obs['sample'] == tempsample]
    adata_ATACtemp = adata_ATAC[adata_ATAC.obs['sample'] == tempsample]
    TG_temp, RE_temp = pseudo_bulk(adata_RNAtemp, adata_ATACtemp, singlepseudobulk)
    TG_pseudobulk = pd.concat([TG_pseudobulk, TG_temp], axis=1)
    RE_pseudobulk = pd.concat([RE_pseudobulk, RE_temp], axis=1)
    RE_pseudobulk[RE_pseudobulk > 100] = 100


# Save pseudobulk and adata matrices
log("Saving pseudobulk data...")
os.makedirs(out_dir + 'data/', exist_ok=True)
adata_ATAC.write(out_dir + 'data/adata_ATAC.h5ad')
adata_RNA.write(out_dir + 'data/adata_RNA.h5ad')
TG_pseudobulk.fillna(0).to_csv(out_dir + 'data/TG_pseudobulk.tsv')
RE_pseudobulk.fillna(0).to_csv(out_dir + 'data/RE_pseudobulk.tsv')
pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv(out_dir + 'data/Peaks.txt', header=None, index=None)


# Train the model
log("Preprocessing and training LINGER model...")
genome = 'hg38'

# in Linger's source code (preprocess.py) l_144 `extract_overlap_regions` forgot to prefix peaks path with outdir
# so we have to chdir to outdir
GRNdir = os.path.abspath(GRNdir) + "/"
out_dir = os.path.abspath(out_dir) + "/"
os.chdir(out_dir)

print(GRNdir)
print(out_dir)

preprocess(TG_pseudobulk, RE_pseudobulk, GRNdir, genome, version, out_dir)

activef='ReLU'
LINGER_tr.training(GRNdir, version, out_dir, activef, species='Human')

# Generate regulatory networks
log("Generating cell population GRNs...")
LL_net.TF_RE_binding(GRNdir, adata_RNA, adata_ATAC, genome, version, out_dir)
LL_net.cis_reg(GRNdir, adata_RNA, adata_ATAC, genome, version, out_dir)
LL_net.trans_reg(GRNdir, version, out_dir, genome)
log("GRNs generation done")
