import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

out_dir = args['out_dir']
path_out = args['path_out']

df = pd.read_csv(f'{out_dir}/cell_population_trans_regulatory.txt', sep='\t', index_col=0)

grn = df.stack().reset_index()
grn.columns = ['target', 'source', 'score']
grn['pval'] = 0.01

# Reorder columns
grn = grn[['source', 'target', 'score', 'pval']]
grn.to_csv(path_out, index=False)