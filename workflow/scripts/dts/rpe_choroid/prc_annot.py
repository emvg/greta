import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--path_frags', nargs='+', required=True)
parser.add_argument('-b', '--path_rna_barcodes', required=True)
parser.add_argument('-o', '--path_output', required=True)
args = vars(parser.parse_args())

path_frags = args['path_frags']
path_rna_barcodes = args['path_rna_barcodes']
path_output = args['path_output']

raw_rna_barcodes = pd.read_csv(path_rna_barcodes, header=None)[0]
rna_barcode = {bc.split('_')[1].split('-')[0]: bc.split('_')[0] for bc in raw_rna_barcodes}

records = []
for fp in path_frags:
    barcodes = pd.read_csv(fp, sep='\t', header=None, comment='#')[3].unique()
    for bc in barcodes:
        celltype, plain_bc = bc.split('_', 1)
        records.append({'celltype': celltype, 'plain_bc': plain_bc})

annot = pd.DataFrame(records)
annot = annot[annot['plain_bc'].isin(rna_barcode)]
annot['batch'] = 'donor' + annot['plain_bc'].map(rna_barcode).astype(str)

annot['barcode'] = annot['celltype'] + '_' + annot['plain_bc']
annot = annot.drop(columns='plain_bc').set_index('barcode')
annot.to_csv(path_output, header=True)