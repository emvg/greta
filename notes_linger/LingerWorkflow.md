# Linger Workflow

## 1. File structure

```
.
├── LINGER_data/
│   └── data_bulk/
│
├── LINGER_output/
│   ├── cell_population_TF_RE_binding.txt
│   ├── cell_population_cis_regulatory.txt
│   ├── cell_population_trans_regulatory.txt
│   └── ...
│
└── data/
    ├── Peaks.txt
    ├── TG_pseudobulk.tsv
    ├── RE_pseudobulk.tsv
    ├── adata_RNA.h5ad
    └── adata_ATAC.h5ad
```

### 2 LINGER_data
Contains the prior knowledge, the bulk GRN trained on bulk data. 

| Files | Used by |
|---|---|
| `all_models_{chr}.pt` | `LINGER_tr.training` — bulk pre-trained model weights loaded as starting point |
| `fisher_{chr}.pt` | `LINGER_tr.training` — Fisher information for EWC regularization |
| `Primary_TF_RE_{chr}.txt` | `LL_net.TF_RE_binding` (baseline mode) |
| `Primary_RE_TG_{chr}.txt` | `LL_net.cis_reg` (baseline mode) |
| `Primary_TF_TG_{chr}.txt` | `LL_net.trans_reg` (baseline mode) |
| `hg19_Peaks_{chr}.bed`, `hg38_Peaks_{chr}.bed` | `preprocess.extract_overlap_regions` — maps user peaks to reference |
| `MotifTarget_Matrix_{chr}.txt`, `MotifTarget_matrix_{chr}.bed` | `preprocess.load_TFbinding` |
| `TF_binding_{chr}.txt` | `preprocess.load_TFbinding` |
| `RE_TG_distance_{chr}.txt` | `LL_net.cis_reg` (baseline mode) |
| `TFName.txt`, `Match2.txt`, `motifWeight.txt` | `preprocess` — TF motif matching |
| `bulk_gene_all.txt`, `all_hg19.txt` | `preprocess.gene_expression` |
| `{chr}_index.txt`, `{chr}_index_all.txt`, `{chr}_gene.txt` | `LINGER_tr.training` — per-chromosome gene/TF/RE indices into bulk model |


## 3 LINGER_output

### 3.1 Preprocessing

| File | Created by | Content |
|---|---|---|
| `Exp.txt` | `preprocess.gene_expression` | Log2 normalized expression for genes overlapping bulk |
| `Symbol.txt` | `preprocess.gene_expression` | Gene symbols corresponding to `Exp.txt` rows |
| `Col.txt` | `preprocess.gene_expression` | Metacell column names |
| `TFexp.txt` | `preprocess.TF_expression` | Expression matrix for TFs only |
| `TFName.txt` | `preprocess.TF_expression` | TF names corresponding to `TFexp.txt` rows |
| `Openness.txt` | `preprocess` | Raw RE pseudobulk chromatin accessibility |
| `TF_binding.txt` | `preprocess.load_TFbinding` | TF-to-RE motif binding scores |
| `index.txt` | `preprocess` — `index_generate` | Per-gene index mapping to RE and TF indices |
| `Region.bed` | `preprocess.extract_overlap_regions` | User peaks reformatted as BED |
| `Region_overlap_{chr}.bed` | `preprocess.extract_overlap_regions` | Intersection of user peaks with reference peaks per chromosome |
| `hg19_Peak_hg19_gene_u.txt` | `preprocess.extract_overlap_regions` | hg19 peak-to-gene mapping after liftover |
| `MotifTarget_hg19_hg38_{chr}.txt` | `preprocess.extract_overlap_regions` | Motif target mapping between assemblies |

### 3.2 Training outputs (`LINGER_tr.training`)

| File | Content |
|---|---|
| `data_merge.txt` | Merged gene table mapping sc genes to bulk gene indices and chromosomes |
| `net_{chr}.pt` | Trained per-gene neural networks for this chromosome |
| `shap_{chr}.pt` | SHAP values for each gene's network — used by `cis_reg` and `trans_reg` |
| `result_{chr}.txt` | Per-gene training result scores |
| `Loss_{chr}.txt` | Training loss curves per gene |


### 3.3 Training outputs (`LL_net`)

| File | Created by | Content |
|---|---|---|
| `chr{N}_cell_population_TF_RE_binding.txt` | `LL_net.TF_RE_binding` | Per-chromosome TF→RE binding matrix (REs × TFs) |
| `cell_population_TF_RE_binding.txt` | `LL_net.TF_RE_binding` | Concatenated TF→RE binding across all chromosomes |
| `cell_population_cis_regulatory.txt` | `LL_net.cis_reg` | RE→gene links (long format: RE, gene, score) |
| `cell_population_trans_regulatory.txt` | `LL_net.trans_reg` | TF→gene GRN matrix (genes × TFs) |



## 1.3 data
| File | Created by | Content |
|---|---|---|
| `data/Peaks.txt` | pseudobulk step | List of ATAC peak IDs in `chr:start-end` format |
| `data/TG_pseudobulk.tsv` | pseudobulk step | Gene expression metacells (genes × metacells) |
| `data/RE_pseudobulk.tsv` | pseudobulk step | Chromatin accessibility metacells (peaks × metacells) |
| `data/adata_RNA.h5ad` | pseudobulk step | RNA AnnData with `barcode`, `label`, `sample`, `gene_ids` fields |
| `data/adata_ATAC.h5ad` | pseudobulk step | ATAC AnnData with same fields |