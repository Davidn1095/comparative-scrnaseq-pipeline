#!/usr/bin/env Rscript
################################################################################
# 05_composition.R — Differential abundance (propeller) across conditions
#
# Computes per-donor proportions of the 25 retained cell types, then runs
# speckle::propeller separately for three pairwise comparisons:
#   HC vs SLE, HC vs SjS, SLE vs SjS
# under the arcsin-sqrt transform. BH-corrects within each comparison.
#
# Inputs:
#   results/objects/atlas_integrated.rds  (Seurat object with meta.data
#     columns: cell_type, donor_id, condition)
#
# Outputs (results/composition/):
#   propeller_HC_vs_SLE.tsv
#   propeller_HC_vs_SjS.tsv
#   propeller_SLE_vs_SjS.tsv
#   propeller_summary.tsv   (one row per cell type x comparison)
################################################################################

# Optional extra R library searched before the default ones, for packages
# installed outside the container image. Set ATLAS_EXTRA_R_LIB to use it.
extra_lib <- Sys.getenv("ATLAS_EXTRA_R_LIB", "")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m)) return(normalizePath(dirname(m)))
  return(".")
}
source(file.path(get_script_dir(), "00_config.R"))
source(file.path(get_script_dir(), "00_utils.R"))

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject)
  library(speckle); library(limma)
  library(dplyr)
})

OUTDIR <- file.path(ATLAS_ROOT, "results", "composition")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

log_msg("Loading atlas...")
obj <- readRDS(file.path(OBJ_DIR, "atlas_integrated.rds"))
meta <- obj@meta.data
rm(obj); gc()
log_msg("Atlas meta loaded: ", nrow(meta), " cells")

stopifnot(all(c("cell_type", "donor_id", "condition") %in% colnames(meta)))

# Restrict to 25 common cell types (the same set used for DE and classifier)
meta <- meta[is_common_celltype(meta$cell_type), ]
log_msg("After 25-cell-type filter: ", nrow(meta), " cells in ",
        length(unique(meta$cell_type)), " cell types")

# Normalise condition labels: expect healthy / sle / sjs (lowercase)
meta$condition <- tolower(as.character(meta$condition))
cond_counts <- table(meta$condition)
log_msg("Condition counts: ", paste0(names(cond_counts), "=", cond_counts,
                                      collapse = ", "))

# Drop conditions not in the 3-class scheme (e.g. sicca, ra remnants)
keep_conds <- c("healthy", "sle", "sjs")
meta <- meta[meta$condition %in% keep_conds, ]
log_msg("After 3-condition filter: ", nrow(meta), " cells, ",
        length(unique(meta$donor_id)), " donors")

run_pair <- function(meta_sub, label_g1, label_g2, comparison_name) {
  # propeller takes clusters, sample (donor), group (condition)
  res <- propeller(clusters = meta_sub$cell_type,
                   sample   = meta_sub$donor_id,
                   group    = meta_sub$condition,
                   transform = "asin")
  # res has columns: BaselineProp.clusters, BaselineProp.Freq, PropMean.<g>, PropRatio,
  #   Tstatistic, P.Value, FDR

  # Build proportions per donor per cell type for explicit mean-proportion columns
  # (propeller returns mean proportions in PropMean.* columns; we also compute
  # fold changes from those for clarity).
  prop_g1_col <- paste0("PropMean.", label_g1)
  prop_g2_col <- paste0("PropMean.", label_g2)

  out <- data.frame(
    cell_type     = rownames(res),
    mean_prop_g1  = res[[prop_g1_col]],
    mean_prop_g2  = res[[prop_g2_col]],
    fold_change   = res[[prop_g2_col]] / res[[prop_g1_col]],
    log2_fc       = log2(res[[prop_g2_col]] / res[[prop_g1_col]]),
    p_value       = res$P.Value,
    fdr           = p.adjust(res$P.Value, method = "BH"),
    stringsAsFactors = FALSE
  )
  colnames(out)[2:3] <- c(paste0("mean_prop_", label_g1),
                          paste0("mean_prop_", label_g2))
  out$significant <- out$fdr < 0.05
  out$comparison  <- comparison_name
  out <- out[order(out$fdr), ]

  fname <- paste0("propeller_", comparison_name, ".tsv")
  write.table(out, file.path(OUTDIR, fname),
              sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("Wrote: ", fname, " (", nrow(out), " rows, ",
          sum(out$significant), " significant at FDR<0.05)")
  out
}

# Three pairwise comparisons
res_hc_sle <- run_pair(meta[meta$condition %in% c("healthy", "sle"), ],
                       "healthy", "sle", "HC_vs_SLE")
res_hc_sjs <- run_pair(meta[meta$condition %in% c("healthy", "sjs"), ],
                       "healthy", "sjs", "HC_vs_SjS")
res_sle_sjs <- run_pair(meta[meta$condition %in% c("sle", "sjs"), ],
                        "sle", "sjs", "SLE_vs_SjS")

# Combined summary: long table with one row per (cell_type, comparison)
to_long <- function(df, g1, g2) {
  data.frame(
    comparison         = df$comparison,
    cell_type          = df$cell_type,
    group1             = g1,
    group2             = g2,
    mean_prop_group1   = df[[paste0("mean_prop_", g1)]],
    mean_prop_group2   = df[[paste0("mean_prop_", g2)]],
    log2_fc            = df$log2_fc,
    p_value            = df$p_value,
    fdr                = df$fdr,
    significant        = df$significant,
    stringsAsFactors = FALSE
  )
}
summary_df <- rbind(
  to_long(res_hc_sle,  "healthy", "sle"),
  to_long(res_hc_sjs,  "healthy", "sjs"),
  to_long(res_sle_sjs, "sle",     "sjs")
)
write.table(summary_df, file.path(OUTDIR, "propeller_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
log_msg("Wrote: propeller_summary.tsv (", nrow(summary_df), " rows)")

log_msg("=== Done ===")
