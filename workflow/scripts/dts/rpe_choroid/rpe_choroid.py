import gzip
import scipy.io
import scipy.sparse
import pandas as pd
import anndata as ad
import mudata as md
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--path_barcodes', required=True)
parser.add_argument('-b', '--path_features', required=True)
parser.add_argument('-c', '--path_matrix', required=True)
parser.add_argument('-d', '--path_peaks', required=True)
parser.add_argument('-e', '--path_annot', required=True)
parser.add_argument('-f', '--path_gid', required=True)
parser.add_argument('-o', '--path_output', required=True)
args = vars(parser.parse_args())

path_barcodes = args['path_barcodes']
path_features = args['path_features']
path_matrix = args['path_matrix']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_geneids = args['path_gid']
path_output = args['path_output']

# barcode :  {celltype}_{plain_bc}
# batch:     'donor2', 'donor11', ...
# celltype:  'Endothelial', ...
obs = pd.read_csv(path_annot, index_col=0)
obs['plain_bc'] = obs.index.str.split('_').str[1]

# Open rna matrix
with gzip.open(path_barcodes, 'rt') as f:
    rna_barcodes = [l.strip() for l in f if l.strip()]
with gzip.open(path_features, 'rt') as f:
    feature_names = [l.strip().split('\t')[1] for l in f if l.strip()]

matrix = scipy.io.mmread(path_matrix).T     
matrix = scipy.sparse.csr_matrix(matrix)            
rna = ad.AnnData(
    X=matrix,
    obs=pd.DataFrame(index=rna_barcodes),
    var=pd.DataFrame(index=feature_names),
)

# Filter faulty gene symbols
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var['gene_ids'] = ensmbls[msk]

# Align rna to annot
rna_bc_set = set(rna.obs_names)
donor_id = obs['batch'].str.replace('donor', '', regex=False)
obs['rna_bc'] = donor_id + '_' + obs['plain_bc'] + '-1'
obs = obs[obs['rna_bc'].isin(rna_bc_set)].copy()

rna = rna[obs['rna_bc']].copy()
rna.obs_names = obs.index        # reindex to '{celltype}_{plain_bc}'
obs = obs.drop(columns=['plain_bc', 'rna_bc'])

atac = ad.read_h5ad(path_peaks)
atac = atac[obs.index].copy()

mdata = md.MuData({'rna': rna, 'atac': atac}, obs=obs)
mdata.write(path_output)