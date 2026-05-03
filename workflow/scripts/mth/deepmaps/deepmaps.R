#!/usr/bin/env Rscript

library(hdf5r)
library(Matrix)
library(spatstat.core)
library(Signac)
library(Seurat)

# Skip conda environment (required inside SIF)
Sys.unsetenv("CONDA_PREFIX")
Sys.unsetenv("CONDA_DEFAULT_ENV")
Sys.setenv(CONDA_EXE = "")
Sys.setenv(RETICULATE_PYTHON = "/home/user/miniconda/bin/python")

library(reticulate)
assignInNamespace(
  "python_munge_path",
  function(python) Sys.getenv("PATH"),
  ns = "reticulate"
)
use_python("/home/user/miniconda/bin/python", required = TRUE)
py_config()

# load deepMaps source code
source("workflow/scripts/mth/deepmaps/scRNA_scATAC1.R")

# ---- Argument parsing -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, required = TRUE) {
  pos <- which(args == flag)
  if (length(pos) == 0 || pos == length(args)) {
    if (required) stop(paste("Missing required argument:", flag))
    return(NULL)
  }
  args[pos + 1]
}

path_mdata <- get_arg("-m")
out_dir    <- get_arg("-d")
path_out   <- get_arg("-o")

# ---- Directory setup --------------------------------------------------------
obj_dir     <- file.path(out_dir, "objects")
jaspar_path <- "dbs/jaspar_data/"
lisa_path   <- paste0(file.path(out_dir, "lisa_data"), "/")

for (d in c(obj_dir, jaspar_path, lisa_path)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Config -----------------------------------------------------------------
species    <- "hg38"
res_clust  <- 1.0
hgt_lr     <- 0.1
hgt_epoch  <- 100
hgt_n_hid  <- 128
hgt_n_heads <- 16
hgt_cuda   <- 0

# =============================================================================
# 1. Read h5mu → Seurat object with RNA + ATAC assays
# =============================================================================
message("[1/9] Reading h5mu and building Seurat object")
f <- H5File$new(path_mdata, mode = "r")

read_sparse <- function(group_path) {
  data    <- f[[paste0(group_path, "/data")]][]
  indices <- f[[paste0(group_path, "/indices")]][]
  indptr  <- f[[paste0(group_path, "/indptr")]][]
  list(data = data, indices = indices, indptr = indptr)
}

# RNA
rna_sp    <- read_sparse("mod/rna/layers/counts")
rna_genes <- f[["mod/rna/var/_index"]][]
rna_cells <- f[["mod/rna/obs/_index"]][]
rna_counts <- sparseMatrix(
  i        = rna_sp$indices + 1,
  p        = rna_sp$indptr,
  x        = rna_sp$data,
  dims     = c(length(rna_genes), length(rna_cells)),
  dimnames = list(rna_genes, rna_cells)
)

# ATAC
atac_sp    <- read_sparse("mod/atac/layers/counts")
atac_peaks <- f[["mod/atac/var/_index"]][]
atac_cells <- f[["mod/atac/obs/_index"]][]
atac_counts <- sparseMatrix(
  i        = atac_sp$indices + 1,
  p        = atac_sp$indptr,
  x        = atac_sp$data,
  dims     = c(length(atac_peaks), length(atac_cells)),
  dimnames = list(atac_peaks, atac_cells)
)

f$close_all()

# Align cells: keep the intersection of both modalities
common_cells <- intersect(rna_cells, atac_cells)
rna_counts   <- rna_counts[,  common_cells]
atac_counts  <- atac_counts[, common_cells]

# Build Seurat object
chrom_assay <- CreateChromatinAssay(
  counts       = atac_counts,
  sep          = c(":", "-"),   # peaks as "chrN:start-end"
  min.cells    = 0,
  min.features = 0
)
obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

exp_assay <- CreateAssayObject(counts = rna_counts)
obj[["RNA"]] <- exp_assay
DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Cell QC filter (matches filterCell() logic from scRNA_scATAC1.r)
obj <- filterCell(obj, data_type = "scRNA_scATAC")

message("Cells after QC: ", ncol(obj))

# =============================================================================
# 2. Peak-gene regulatory potential (MAESTRO)
# =============================================================================
message("[2/9] Computing peak-to-gene regulatory potential")
ATAC_gene_peak <- CalGenePeakScore(
  peak_count_matrix = obj@assays$ATAC@counts,
  organism          = "GRCh38"
)

# =============================================================================
# 3. Gene Activity Score (GAS)
# =============================================================================
message("[3/9] Computing Gene Activity Score (GAS)")
gas_res   <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak,
                               obj            = obj,
                               method         = "wnn")
GAS <- gas_res[[1]]
obj <- gas_res[[2]]

# =============================================================================
# 4. HGT model (Python via reticulate)
# =============================================================================
hgt_cache <- file.path(obj_dir, "HGT_result.rds")
if (file.exists(hgt_cache)) {
  message("[4/9] Loading cached HGT result")
  HGT_result <- readRDS(hgt_cache)
} else {
  message("[4/9] Running HGT model")
  HGT_result <- run_HGT(
    GAS        = as.matrix(GAS),
    result_dir = obj_dir,
    data_type  = "scRNA_scATAC",
    lr         = hgt_lr,
    epoch      = hgt_epoch,
    n_hid      = hgt_n_hid,
    n_heads    = hgt_n_heads,
    cuda       = hgt_cuda
  )
  saveRDS(HGT_result, hgt_cache)
}

# =============================================================================
# 5. Cell clustering on HGT embedding
# =============================================================================
message("[5/9] Clustering cells")
cell_hgt_matrix <- HGT_result[["cell_hgt_matrix"]]
rownames(cell_hgt_matrix) <- colnames(GAS)
obj             <- obj[, colnames(GAS)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS), ]

obj@reductions[["HGT"]] <- CreateDimReducObject(
  embeddings = cell_hgt_matrix, key = "HGT_", assay = "RNA"
)

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunUMAP(obj,
               reduction = "HGT", dims = 1:ncol(cell_hgt_matrix),
               reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
obj <- FindNeighbors(obj,
                     reduction = "HGT", graph.name = "HGT_snn",
                     dims = 1:ncol(cell_hgt_matrix))
obj <- FindClusters(obj, graph.name = "HGT_snn", resolution = res_clust)

graph.out <- as.factor(obj$seurat_clusters)
names(graph.out) <- colnames(obj)

saveRDS(obj, file.path(obj_dir, "seurat_obj.rds"))

# =============================================================================
# 6. Cell-type active gene modules
# =============================================================================
message("[6/9] Building gene modules")
co_cache <- file.path(obj_dir, "co.rds")
if (file.exists(co_cache)) {
  co <- readRDS(co_cache)
} else {
  co <- get_gene_module(obj = obj, GAS = GAS,
                        att = HGT_result[["attention"]], method = "SFP") #SFP
  saveRDS(co, co_cache)
}

write_GM(co, lisa_path)
message("    Gene modules written → running LISA")

system(
  paste0(
    "/home/user/miniconda/bin/python workflow/scripts/mth/deepmaps/run_lisa.py",
    " --path ", lisa_path,
    " --lisa_data dbs/lisa_data/",
    " --species hg38"
  )
)

# =============================================================================
# 7. Promoter filter + JASPAR/LISA regulon inference
# =============================================================================
message("[7/9] Filtering promoter accessibility and inferring regulons")
gene_peak_pro <- AccPromoter(obj       = obj,
                             gene_peak = ATAC_gene_peak,
                             GAS       = GAS,
                             species   = species)

pre_reg    <- Calregulon(GAS           = GAS,
                         co            = co,
                         gene_peak_pro = gene_peak_pro,
                         species       = species,
                         jaspar_path   = jaspar_path,
                         lisa_path     = lisa_path)
BA_score      <- pre_reg[[1]]
ct_regulon_v1 <- pre_reg[[2]]
TFinGAS       <- pre_reg[[3]]

peak_TF <- uni(gene_peak_pro = gene_peak_pro, BA_score = BA_score)

# =============================================================================
# 8. RI scores + RAS + master TFs
# =============================================================================
message("[8/9] Computing RI, RAS, and master TFs")
RI_C <- RI_cell(obj           = obj,
                ct_regulon    = ct_regulon_v1,
                GAS           = GAS,
                gene_peak_pro = gene_peak_pro,
                peak_TF       = peak_TF,
                graph.out     = graph.out)

reg_res    <- calRAS(RI_C       = RI_C,
                     ct_regulon = ct_regulon_v1,
                     graph.out  = graph.out,
                     TFinGAS    = TFinGAS)
RAS_CT     <- reg_res[[1]]
RI_CT      <- reg_res[[2]]
ct_regulon <- reg_res[[3]]
RAS_C1     <- reg_res[[4]]

masterTF <- masterFac(ct_regulon = ct_regulon, RI_CT = RI_CT)
TF_cen   <- masterTF[[1]]
gene_cen <- masterTF[[2]]
network  <- masterTF[[3]]

# Persist all objects so grn.r can run independently if needed
saveRDS(ct_regulon,    file.path(obj_dir, "ct_regulon.rds"))
saveRDS(RI_CT,         file.path(obj_dir, "RI_CT.rds"))
saveRDS(TF_cen,        file.path(obj_dir, "TF_cen.rds"))
saveRDS(gene_cen,      file.path(obj_dir, "gene_cen.rds"))
saveRDS(network,       file.path(obj_dir, "network.rds"))
saveRDS(peak_TF,       file.path(obj_dir, "peak_TF.rds"))
saveRDS(gene_peak_pro, file.path(obj_dir, "gene_peak_pro.rds"))
