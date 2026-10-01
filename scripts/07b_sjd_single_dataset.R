#!/usr/bin/env Rscript
################################################################################
# 07b_sjd_single_dataset.R — SjD against healthy controls within one dataset
#
# Robustness check for the SjD pathway results. The pooled SjD analysis (06, 07)
# compares SjD patients with the controls of both SjD datasets, and the healthy
# controls of GSE157278 carry about twice the dissociation-gene reads of every
# other group. GSE253568 has patients and controls at the same dissociation load,
# so SjD is compared with healthy controls within GSE253568 alone, with the
# pipeline's gene exclusions (dissociation genes included), and the result is set
# beside the pooled SjD and the SLE results in Supplementary Table 7 (script 11).
#
# Everything is as in 06 and 07, within the one dataset: pseudobulk counts from 03,
# the same MIN_DONORS / MIN_CELLS / MIN_CELLS_CLASS / MIN_CELLS_PER_DONOR floors with
# cell counts from the atlas, DESeq2 with design ~ condition (one dataset, so no
# dataset covariate), the gene exclusion applied after results(), and fgsea on the
# Wald statistic with the 50 Hallmark sets (minSize 15, maxSize 500, seed 42).
#
# Inputs:  data/SjS/<ROBUSTNESS_DATASET>/derived/pseudobulk/counts/ (03),
#          results/objects/atlas_integrated.rds (04a)
# Outputs (results/plots/differential_expression/):
#   de_sjs_<dataset>/deseq2_<cell_type>.tsv   per cell type, 06's layout
#   de_sjs_<dataset>/tests.tsv                donors and cells per group, tested or not
#   gsea_sjs_<dataset>.tsv                    fgsea, 07's layout
#
# Environment: ROBUSTNESS_DATASET (default GSE253568); the floors as in 06.
################################################################################
SCRIPT_DIR <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (SCRIPT_DIR == "." || SCRIPT_DIR == "") {
  SCRIPT_DIR <- tryCatch(
    dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])),
    error = function(e) getwd()
  )
  SCRIPT_DIR <- sub("^--file=", "", SCRIPT_DIR)
}
source(file.path(SCRIPT_DIR, "00_config.R"))
source(file.path(SCRIPT_DIR, "00_utils.R"))

suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(Matrix)
})

DATASET <- get_env("ROBUSTNESS_DATASET", "GSE253568")
PB_DIR  <- file.path(ATLAS_ROOT, "data", "SjS", DATASET, "derived", "pseudobulk", "counts")
DE_DIR  <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression")
OUTDIR  <- file.path(DE_DIR, paste0("de_sjs_", tolower(DATASET)))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
log_msg("SjD vs HC within ", DATASET, " | MIN_DONORS=", MIN_DONORS, " MIN_CELLS=", MIN_CELLS,
        " MIN_CELLS_CLASS=", MIN_CELLS_CLASS, " MIN_CELLS_PER_DONOR=", MIN_CELLS_PER_DONOR)

counts <- readRDS(file.path(PB_DIR, "pb_counts.rds"))
meta <- read.delim(file.path(PB_DIR, "pb_meta.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(meta$pb_id, colnames(counts)))
ctrl_labels <- c("normal", "hd", "healthy", "control")
meta$condition <- ifelse(tolower(meta$condition) %in% ctrl_labels, "control",
                  ifelse(tolower(meta$condition) %in% c("sjs", "case", "disease"), "case", NA))
meta <- meta[!is.na(meta$condition), ]

# Cells per donor and cell type, from the atlas as in 06
atlas_meta <- load_atlas()@meta.data
atlas_meta <- atlas_meta[atlas_meta$dataset_id == DATASET, ]
donor_ncells <- aggregate(list(n_cells = rep(1L, nrow(atlas_meta))),
                          by = list(cell_type = atlas_meta$cell_type, donor = atlas_meta$donor),
                          FUN = length)
hallmark <- load_msigdb_hallmark(Sys.getenv("MSIGDB_RDS", file.path(tools::R_user_dir("msigdbr", which = "data"), "msigdb.2025.1.Hs.rds")))

tests <- list(); gsea <- list()
for (ct in sort(unique(meta$cell_type[is_common_celltype(meta$cell_type)]))) {
  ct_meta <- meta[meta$cell_type == ct, ]
  ct_donor <- donor_ncells[donor_ncells$cell_type == ct, ]
  ct_meta <- ct_meta[ct_meta$donor_id %in% ct_donor$donor[ct_donor$n_cells >= MIN_CELLS_PER_DONOR], ]
  n_case <- sum(ct_meta$condition == "case"); n_ctrl <- sum(ct_meta$condition == "control")
  c_case <- sum(ct_donor$n_cells[ct_donor$donor %in% ct_meta$donor_id[ct_meta$condition == "case"]])
  c_ctrl <- sum(ct_donor$n_cells[ct_donor$donor %in% ct_meta$donor_id[ct_meta$condition == "control"]])
  status <- if (n_case < MIN_DONORS || n_ctrl < MIN_DONORS) "skipped: too few donors" else
            if (c_case + c_ctrl < MIN_CELLS) "skipped: too few total cells" else
            if (c_case < MIN_CELLS_CLASS || c_ctrl < MIN_CELLS_CLASS) "skipped: too few cells per class" else "tested"
  tests[[ct]] <- data.frame(cell_type = ct, n_donors_SjS = n_case, n_donors_healthy = n_ctrl,
                            n_cells_SjS = c_case, n_cells_healthy = c_ctrl, status = status)
  if (status != "tested") { log_msg("  ", ct, ": ", status, " (", n_case, " vs ", n_ctrl, ")"); next }

  y <- as.matrix(counts[, ct_meta$pb_id, drop = FALSE]); storage.mode(y) <- "integer"
  y <- y[rowSums(y) > 0, , drop = FALSE]
  ct_meta$condition <- factor(ct_meta$condition, levels = c("control", "case"))
  dds <- DESeq(DESeqDataSetFromMatrix(countData = y, colData = ct_meta, design = ~ condition), quiet = TRUE)
  res <- results(dds, contrast = c("condition", "case", "control"))
  df <- data.frame(gene = rownames(res), log2FoldChange = res$log2FoldChange, pvalue = res$pvalue,
                   padj = res$padj, baseMean = res$baseMean, lfcSE = res$lfcSE, stat = res$stat,
                   stringsAsFactors = FALSE)
  df <- df[!is_excluded_gene(df$gene), ]
  df <- df[!is.na(df$pvalue), ]
  df <- df[order(df$padj, df$pvalue), ]
  write.table(df, file.path(OUTDIR, paste0("deseq2_", sanitize_name(ct), ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  r <- df[!is.na(df$stat), ]
  if (nrow(r) >= 100) {
    ranks <- sort(setNames(r$stat, r$gene), decreasing = TRUE)
    set.seed(42)
    g <- fgsea(pathways = hallmark, stats = ranks, minSize = 15, maxSize = 500, nproc = N_CORES)
    g$leadingEdge <- vapply(g$leadingEdge, paste, character(1), collapse = ";")
    g$disease <- "SjS"; g$dataset <- DATASET; g$cell_type <- ct
    gsea[[ct]] <- as.data.frame(g)
  }
  log_msg("  ", ct, ": SjD ", n_case, " vs HC ", n_ctrl, " donors, ",
          sum(df$padj < PADJ_THRESH & abs(df$log2FoldChange) >= LFC_THRESH, na.rm = TRUE), " DEGs")
}
tests <- do.call(rbind, tests)
write.table(tests, file.path(OUTDIR, "tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
out <- do.call(rbind, gsea)
write.table(out, file.path(DE_DIR, paste0("gsea_sjs_", tolower(DATASET), ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg(sum(tests$status == "tested"), " of ", nrow(tests), " cell types tested; wrote ", OUTDIR,
        " and gsea_sjs_", tolower(DATASET), ".tsv")
