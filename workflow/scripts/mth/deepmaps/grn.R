library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, required = TRUE) {
  pos <- which(args == flag)
  if (length(pos) == 0 || pos == length(args)) {
    if (required) stop(paste("Missing required argument:", flag))
    return(NULL)
  }
  args[pos + 1]
}

out_dir  <- get_arg("-d")
path_out <- get_arg("-o")
obj_dir  <- file.path(out_dir, "objects")

# Load required objects from DeepMaps run
message("[grn.R] Loading objects from ", obj_dir)
ct_regulon    <- readRDS(file.path(obj_dir, "ct_regulon.rds"))
RI_CT         <- readRDS(file.path(obj_dir, "RI_CT.rds"))
gene_peak_pro <- readRDS(file.path(obj_dir, "gene_peak_pro.rds"))
peak_TF       <- readRDS(file.path(obj_dir, "peak_TF.rds"))

message("  ", length(ct_regulon), " regulons, ",
        nrow(RI_CT), " (TF,TG) pairs x ", ncol(RI_CT), " clusters")

# Build TF–RE–TG triplets 
message("[grn.R] Building TF-RE-TG triplets")
triplets <- list()
idx <- 0L

for (reg_name in names(ct_regulon)) {
  parts <- strsplit(reg_name, "_")[[1]]
  tf <- parts[1]; ct <- parts[2]

  for (target in ct_regulon[[reg_name]]) {
    pair_key <- paste0(tf, "_", target)

    # cluster-level RI score
    ri <- if (pair_key %in% rownames(RI_CT) && ct %in% colnames(RI_CT)) {
      RI_CT[pair_key, ct]
    } else { NA_real_ }
    if (is.na(ri) || ri <= 0) next

    # Supporting CREs: peaks in C_ti = {k : w_ik > 0 and b_kt > 0}.
    # Guaranteed non-empty since geneTF_it > 0 is required for ct_regulon membership.
    gene_peaks <- names(which(gene_peak_pro[target, ] > 0))
    tf_peaks   <- names(which(peak_TF[, tf] > 0))
    cres       <- intersect(gene_peaks, tf_peaks)

    for (cre in cres) {
      idx <- idx + 1L
      triplets[[idx]] <- data.frame(
        source                 = tf,
        cre                    = cre,
        target                 = target,
        cluster                = ct,
        score                  = ri,
        peak_to_gene_potential = gene_peak_pro[target, cre],
        peak_TF_affinity       = peak_TF[cre, tf],
        stringsAsFactors       = FALSE
      )
    }
  }
}

grn_celltype <- bind_rows(triplets)
message("  ", nrow(grn_celltype), " cluster-resolved triplets")

# Cell-population GRN: max score across clusters per (TF, RE, TG) 
message("[grn.R] Pooling to cell-population level")
grn <- grn_celltype %>%
  group_by(source, cre, target) %>%
  summarise(
    score                  = max(score),
    peak_to_gene_potential = first(peak_to_gene_potential),
    peak_TF_affinity       = first(peak_TF_affinity),
    .groups                = "drop"
  )

write.csv(grn, path_out, row.names = FALSE)
message("  Population GRN -> ", path_out, " (", nrow(grn), " triplets)")

ct_out <- sub("\\.csv$", "_celltype.csv", path_out)
write.csv(grn_celltype, ct_out, row.names = FALSE)
message("  Cell-type GRN  -> ", ct_out,   " (", nrow(grn_celltype), " triplets)")

message("[grn.R] Summary")
message("  TFs:     ", length(unique(grn$source)))
message("  Targets: ", length(unique(grn$target)))
message("  CREs:    ", length(unique(grn$cre)))
message("  Per-cluster edge counts:")
ct_counts <- table(grn_celltype$cluster)
for (ct in sort(names(ct_counts))) {
  message("    ", ct, ": ", ct_counts[ct])
}