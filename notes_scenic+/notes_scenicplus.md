
## Greta complete Scenic+ workflow

Going through rule *mdl_o_scenicplus* steps

### 1. Inputs

(all files are located in `rules` folder)

- From `dts/general.smk`
    - `mdata` : mdata.h5mu (multiomics data) 
- From `dbs/cre.smk`:
    - `blist` : hg38_cre_blacklist.bed
- From `dbs/gen.smk`:
    - `rnk`: regions_vs_motifs.rankings
    - `man` : motif annotations
    - `scr`: regions_vs_motifs.scores
    - `ann`: annotation.tsv
    - `csz`: chromsizes.tsv

### 2. Outputs 

- dir `dts/{org}/{dat}/cases/{case}/runs/scenicplus/`
- output file `dts/.../runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv`

### 3. Script execution

From `workflow/scripts/mth/scenicplus/o_mdl.sh`, divided in 4 parts.

#### 3.1 *pre* (preprocessing)

From `workflow/scripts/mth/scenicplus/topics.py`, it will process the multiome data `madata`

```bash
mdata.mod['rna']   → RNA matrix
mdata.mod['atac']  → ATAC matrix
mdata.obs          → cells
```

Both modalities get split into 
- `rna.h5ad` : single file for RNA  [cells x gex]
- `cistopic_obj` : made of ATAC data only   [cells x regions]

Topic = Group of co-accessible peaks, we can then find which peaks belong to which topics and which topics are active in which cell.

**Flow** : Raw multiome data → builds ATAC topics → merges with RNA → outputs a SCENIC+ ready `mdata.h5mu` file (with peaks)



#### 3.2 `p2g` - CRE (peak) to gene

`scenicplus prepare_data search_space ... `

Links peaks to genes, creates a search space file `space.tsv`. <br>
**Flow** : Takes peak positions + gene locations → build a list of which peaks are near which genes

| Peak             | Gene | Distance |
|------------------|------|----------|
| chr1:1000-1200   | TP53 | -500     |
| chr1:2000-2200   | TP53 | +1500    |




#### 3.3 `tfb` - TF to CRE binding prediction

##### 3.3.1 `pycistopic topic_modeling mallet binarize`
This step turns fuzzy ATAC topics into concrete peak sets using automatic thresholding.

Gives `TopicX_regions.txt`, where each file = list of peaks active in topic X

##### 3.3.2 `scenicplus grn_inference motif_enrichment_cistarget`

We then connect TFs to topics : which TF motifs (TFs bind to specific DNA sequences called motifs) are enriched in each topic's regions ? 
This is done by taking each topic's region and scanning them for TF motifs (using Pycistarget : discovers potential TFBSs in candidates enhancers)

Gives `cistarget.hdf5` : Topic × TF enrichment scores

**Flow** : This step scans the topic-specific peaks for TF motifs using a reference motif database and outputs which TFs are enriched in each topic.


##### 3.3.3 `scenicplus grn_inference motif_enrichment_dem`

DEM = Differential Enrichment Motifs. 
- Pycistarget finds known top motifs (using precomputed db `rnk`)
- DEM finds additional motifs enriched (using our dataset)


##### 3.3.4 `scenicplus prepare_data prepare_menr`

Creates the “map of TFs to their candidate regulatory regions”, which is the input for TF → gene inference.

Gives :

- `tfs.txt` : List of TFs that passed motif enrichment thresholds
- `direct.h5ad` :	Annotated dataset with direct TF → peak associations
- `extended.h5ad`	: Annotated dataset with extended TF → peak associations (orthology or additional motifs)



#### 3.4 `mdl` - Modeling (final TF-Gene links)

From `scenicplus grn_inference TF_to_gene`

So far we have : 
 
- Topics (groups of co-accessible peaks) - 3.1
- Peak &rarr; TF motif annotations (`direct.h5ad` + `extended.h5ad`) - 3.3.4
- Peak &rarr; gene candidate space (`space.tsv`) - 3.2

Now, we want to link TFs to TGs :

- For each TF in `tfs.txt` we can get its assigned regions from `direct.h5ad`. We can then find genes potentially regulated by the found regions with `p2g` in 3.2

Gives : `tg_adg.tsv`, TF-TG-Score matrix

| TF    | TG   | score |
|-------|------|-------|
| GATA1 | HBB  | 0.82  |
| CEBPA | CSF1 | 0.75  |




#### 3.5 Fine tuned `mdl`

From `grn_inference eGRN` 

eGRN takes all previous steps — topics, motif enrichment, peak → gene mappings, TF → gene regression — and produces a high-confidence, final TF → gene regulatory network

Gives `egrn.tsv`, the final network :

| TF    | TG         | Score | Source_region  |
|-------|------------|-------|----------------|
| GATA1 | HBB        | 0.82  | chr11:…        |
| CEBPA | CSF1       | 0.75  | chr20:…        |


Source_region = peak through which TF regulates the gene


Then `python workflow/scripts/mth/scenicplus/egrn.py $new_dir/egrn.tsv $path_out` to output the final flat csv : 

| source | cre      | target | score | pval |
|--------|----------|--------|-------|------|
| TF1    | chr1-23 | GeneA  | 0.82  | 0.01 |
| TF2    | chr2-45 | GeneB  | 0.75  | 0.01 |

Gives : `dts/hg38/pbmc10k/cases/all/runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv`



### 4. Flow diagram 


Raw multiome data (RNA + ATAC, .h5mu)  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 1: Topics (topics.py) │  
└──────────────────────────┘  
        │  
        ▼  
$new_dir/mdata.h5mu   ← processed multiome  
$new_dir/rna.h5ad     ← RNA data  
$new_dir/bool_mat.mtx ← binary peak matrix  
$new_dir/cres.txt     ← peak names  
$new_dir/brcs.txt     ← cell names  
$new_dir/cistopic_obj.pkl ← cisTopic object  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 2: P2G (search_space) │  
└──────────────────────────┘  
        │  
        ▼  
$new_dir/space.tsv  ← peak → gene candidate regions  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 3: Motif Enrichment │  
└──────────────────────────┘  
        │  
        ▼  
cistarget.hdf5  ← database-based motif enrichment  
dem.hdf5        ← differential enrichment motif  
$new_dir/topics/otsu/...  ← topic-specific peak binarization  
        │  
        ▼  
┌───────────────────────────────┐  
│ Step 4: prepare_menr           │  
└───────────────────────────────┘  
        │  
        ▼  
direct.h5ad   ← high-confidence TF → region links  
extended.h5ad ← extended TF → region links  
tfs.txt       ← list of TFs to consider  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 5: TF → Gene         │  
└──────────────────────────┘  
        │  
        ▼  
tg_adj.tsv  ← TF → gene adjacency table (scores)  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 6: eGRN              │  
└──────────────────────────┘  
        │  
        ▼  
egrn.tsv   ← final TF → gene network (raw)  
        │  
        ▼  
┌──────────────────────────┐  
│ Step 7: eGRN formatting  │  
│ (egrn.py script)         │  
└──────────────────────────┘  
        │  
        ▼  
o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv  
    ← final, ranked, scored CSV for downstream analysis

