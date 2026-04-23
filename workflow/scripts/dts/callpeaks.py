import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import snapatac2 as snap
import anndata as ad
from snapatac2.datasets import _datasets, datasets
from pathlib import Path
import pandas as pd
import numpy as np
import scipy
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f','--path_frags', nargs='+', required=True)
parser.add_argument('-a','--path_annot', required=True)
parser.add_argument('-t','--path_tmp', required=True)
parser.add_argument('-n','--n_jobs', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())
path_frags = args['path_frags']
path_annot = args['path_annot']
path_tmp = args['path_tmp']
n_jobs = int(args['n_jobs'])
path_output = args['path_output']

if __name__ == '__main__':
    print(f'[1/8] Starting callpeaks.py')
    print(f'      frags={path_frags}')
    print(f'      annot={path_annot}')
    print(f'      tmp={path_tmp}')
    print(f'      n_jobs={n_jobs}')
    print(f'      output={path_output}', flush=True)

    """
    #path_tmp="/workdir/vangysel"
    # Change default cache dir
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    """

    os.environ["TMPDIR"] = path_tmp
    _datasets = datasets()
    _datasets.path = Path(path_tmp)
    print(f'[2/8] Temp dir ready: {path_tmp}', flush=True)

    # Find sample_ids
    sample_ids = [os.path.basename(p).split('.')[0].replace('_atac_fragments', '') for p in path_frags]
    tmp_files = [os.path.join(path_tmp, p + '.frags.h5ad') for p in sample_ids]
    print(f'[3/8] Sample IDs: {sample_ids}', flush=True)

    # Read and create h5ad fragment files
    print(f'[4/8] Importing fragments (this may take a while)...', flush=True)
    _ = snap.pp.import_data(
        path_frags,
        chrom_sizes=snap.genome.hg38,
        file=tmp_files,
        tempdir=path_tmp,
        sorted_by_barcode=False,
        n_jobs=n_jobs
    )
    del _
    print(f'[4/8] Done importing fragments', flush=True)

    # Filter by annotation
    print(f'[5/8] Filtering by annotation...', flush=True)
    annot = pd.read_csv(path_annot, index_col=0)
    uns = None
    type_frags = None
    lst_obs = []
    lst_frags = []
    for sample_id, tmp_file in zip(sample_ids, tmp_files):
        print(f'      Processing sample: {sample_id}', flush=True)
        tmp = ad.read_h5ad(tmp_file, backed='r')
        tmp.obs.index = [barcode.split('-1')[0].replace('_atac_fragments', '') for barcode in tmp.obs.index]
        obs = pd.merge(tmp.obs, annot, left_index=True, right_index=True)
        print(f'      {len(obs)} barcodes matched annotation', flush=True)
        uns = tmp.uns['reference_sequences']
        type_frags = list(tmp.obsm.keys())[0]
        lst_frags.append(tmp[obs.index, :].obsm[type_frags])
        lst_obs.append(obs)

    atac = ad.AnnData(obs=pd.concat(lst_obs))
    atac.uns = {'reference_sequences': uns}
    atac.obsm[type_frags] = scipy.sparse.vstack(lst_frags)
    print(f'[5/8] AnnData assembled: {atac.shape}', flush=True)

    # Call and merge peaks
    print(f'[6/8] Calling peaks with macs3 (n_jobs={n_jobs})...', flush=True)

    #atac.obs['batch'] = atac.obs['batch'].astype(str)
    snap.tl.macs3(atac, groupby='celltype', replicate='batch', n_jobs=n_jobs, tempdir=path_tmp)
    peaks = snap.tl.merge_peaks(atac.uns['macs3'], snap.genome.hg38)
    print(f'[6/8] Peaks called and merged: {len(peaks)} peaks', flush=True)

    print(f'[7/8] Building peak matrix...', flush=True)
    atac = snap.pp.make_peak_matrix(atac, use_rep=peaks['Peaks'])
    print(f'[7/8] Peak matrix shape: {atac.shape}', flush=True)

    # Clean
    del atac.obs
    del atac.var

    # Update format peaks
    new_var_names = []
    for p in atac.var_names:
        p = p.replace(':', '-')
        seq, start, end = p.split('-')
        end = int(end) - 1
        p = '{0}-{1}-{2}'.format(seq, start, end)
        new_var_names.append(p)
    atac.var_names = new_var_names

    # Write
    print(f'[8/8] Writing output to {path_output}...', flush=True)
    atac.write(path_output)
    print(f'[8/8] Done!', flush=True)

    os._exit(0)  # Add this else it gets stuck