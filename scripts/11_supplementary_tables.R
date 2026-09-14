#!/usr/bin/env Rscript
################################################################################
# 11_supplementary_tables.R — Regenerate Supplementary Tables 1-6
#
# Purpose: assemble every supplementary table from current pipeline outputs.
# These six files previously had no generator in the repository and had drifted
# 60-76 days behind the results they summarise.
#
# Schemas below were read off the existing files and are reproduced exactly, so
# regenerating is a drop-in replacement. The one exception is S6, which gained
# three per-class columns and dropped the fold-SD column when Fig 4 moved to
# per-class SHAP.
#
#   S1 lineage mapping  <- LINEAGE_MAP (00_config.R), restricted to the 25
#                          common cell types
#   S2 composition      <- results/composition/propeller_summary.tsv (script 05)
#                          (byte-identical schema; verified identical to the
#                          previous supplementary file)
#   S3 DE SLE           <- results/plots/.../de_sle/deseq2_*.tsv (script 06), gzipped
#   S4 DE SjS           <- results/plots/.../de_sjs/deseq2_*.tsv (script 06), gzipped
#   S6 GSEA             <- results/plots/.../gsea_all_results.tsv (script 07)
#   S5 discordant       <- the two DE tables, restricted to genes significant in
#                          both diseases within a cell type with opposite sign
#   S7 SHAP             <- results/donor_classifiers_<method>/shap_glmnet/ (script 09):
#                          class-averaged mean |SHAP| (ordering) + per-class
#                          HC / SLE / SjS columns; no fold SDs
#
# Numbering follows first mention in main.tex, counting figure captions:
# the discordant table is first cited in Methods (Differential expression),
# after S3 and S4 and before the GSEA table, so it is S5.
#
# Output: manuscript/supplementary/
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
  library(dplyr)
})

# DENOISE_METHOD IS MANDATORY - no default.
#
# It used to default to "none" HERE and to "limma_modulescore" in 09, 10 and 11, so a
# run that did not set it had the PRODUCER writing results/donor_classifiers/ while
# every CONSUMER read results/donor_classifiers_limma_modulescore/. Nothing failed:
# 08 trained and wrote a complete tree, 09 then read the OTHER tree and reported
# success off stale models. Silent defaults that disagree across a pipeline are worse
# than no default, so every caller must now say which configuration it wants.
DENOISE_METHOD <- tolower(Sys.getenv("DENOISE_METHOD", ""))
if (!nzchar(DENOISE_METHOD)) {
  stop("DENOISE_METHOD is not set. It must be one of none / scvi / limma / ",
       "limma_modulescore, and it selects BOTH the analysis and the results tree. ",
       "The manuscript configuration is limma_modulescore. Refusing to guess, because ",
       "the producer and the consumers used to guess differently.")
}
DE_DIR  <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression")
CLF_DIR <- file.path(ATLAS_ROOT, "results", paste0("donor_classifiers_", DENOISE_METHOD))
OUTDIR  <- file.path(ATLAS_ROOT, "manuscript", "supplementary")
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

written <- character(0)
note <- function(f, n) {
  log_msg("  wrote ", basename(f), " (", n, " rows)")
  written <<- c(written, basename(f))
}

## ---- S1: lineage mapping ---------------------------------------------------
log_msg("Supplementary Table 1: lineage mapping")
lm <- data.frame(cell_type = names(LINEAGE_MAP), lineage = unname(LINEAGE_MAP),
                 stringsAsFactors = FALSE)
lm <- lm[is_common_celltype(lm$cell_type), , drop = FALSE]
# Report the space-normalised COMMON_CELL_TYPES spelling, matching the other
# supplementary tables and the figure labels. LINEAGE_MAP carries both the
# slash and space spellings of Th1/Th17, so de-duplicate after normalising.
lm$cell_type <- gsub("[-/]", " ", lm$cell_type)
lm <- lm[!duplicated(lm$cell_type), , drop = FALSE]
# Progenitors are mapped to "Other" in LINEAGE_MAP but presented as their own
# lineage in the 8-lineage grouping (mirrors 10_manuscript_figures.R).
lm$lineage[lm$lineage == "Other"] <- "Progenitor"
lm <- lm[order(lm$lineage, lm$cell_type), , drop = FALSE]
f1 <- file.path(OUTDIR, "supp_table_lineage_mapping.tsv")
write.table(lm, f1, sep = "\t", quote = FALSE, row.names = FALSE)
note(f1, nrow(lm))
if (nrow(lm) != length(COMMON_CELL_TYPES)) {
  log_msg("  WARNING: expected ", length(COMMON_CELL_TYPES), " cell types, got ", nrow(lm))
}

## ---- S2: composition -------------------------------------------------------
log_msg("Supplementary Table 2: composition (propeller)")
psum <- file.path(ATLAS_ROOT, "results", "composition", "propeller_summary.tsv")
if (file.exists(psum)) {
  d2 <- read.delim(psum, stringsAsFactors = FALSE)
  f2 <- file.path(OUTDIR, "supp_table_composition.tsv")
  write.table(d2, f2, sep = "\t", quote = FALSE, row.names = FALSE)
  note(f2, nrow(d2))
} else {
  log_msg("  SKIPPED: ", psum, " not found — run 05_composition.R")
}

## ---- S3 / S4: full DE tables ----------------------------------------------
collect_de <- function(disease) {
  dd <- file.path(DE_DIR, paste0("de_", disease))
  fs <- list.files(dd, pattern = "^deseq2_.*\\.tsv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(sort(fs), function(f) {
    ct <- gsub("_", " ", sub("^deseq2_", "", sub("\\.tsv$", "", basename(f))))
    df <- read.delim(f, stringsAsFactors = FALSE)
    if (!nrow(df)) return(NULL)
    df$cell_type <- ct
    # unfiltered: every tested gene, ordered by significance within cell type
    df <- df[order(df$padj, df$pvalue), , drop = FALSE]
    df[, c("cell_type", "gene", "baseMean", "log2FoldChange",
           "lfcSE", "stat", "pvalue", "padj"), drop = FALSE]
  }))
}
for (spec in list(c("sle", "supp_table_de_sle.tsv.gz", "3"),
                  c("sjs", "supp_table_de_sjs.tsv.gz", "4"))) {
  log_msg("Supplementary Table ", spec[3], ": full DE results (", toupper(spec[1]), ")")
  d <- collect_de(spec[1])
  if (is.null(d)) { log_msg("  SKIPPED: no DE tables in de_", spec[1]); next }
  f <- file.path(OUTDIR, spec[2])
  con <- gzfile(f, "w")
  write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  note(f, nrow(d))
}

## ---- S5: discordant DEGs ---------------------------------------------------
# Genes called significant in both diseases within the same cell type but with
# opposite sign. Cited from the Fig 2 caption, where the count alone is given.
log_msg("Supplementary Table 5: discordant DEGs")
sle <- collect_de("sle"); sjs <- collect_de("sjs")
if (is.null(sle) || is.null(sjs)) {
  log_msg("  SKIPPED: DE tables missing - run 06_differential_analysis.R")
} else {
  keep <- function(d) d[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1 &
                        !grepl(EXCLUDED_GENE_REGEX, d$gene, ignore.case = TRUE),
                        c("cell_type", "gene", "log2FoldChange", "padj")]
  a <- keep(sle); b <- keep(sjs)
  names(a)[3:4] <- c("log2FC_SLE", "padj_SLE")
  names(b)[3:4] <- c("log2FC_SjS", "padj_SjS")
  m <- merge(a, b, by = c("cell_type", "gene"))
  m <- m[sign(m$log2FC_SLE) != sign(m$log2FC_SjS), ]
  m <- m[order(m$cell_type, m$gene), ]
  m$log2FC_SLE <- round(m$log2FC_SLE, 3); m$log2FC_SjS <- round(m$log2FC_SjS, 3)
  m$padj_SLE <- signif(m$padj_SLE, 3);    m$padj_SjS <- signif(m$padj_SjS, 3)
  f7 <- file.path(OUTDIR, "supp_table_discordant_degs.tsv")
  write.table(m, f7, sep = "\t", quote = FALSE, row.names = FALSE)
  note(f7, nrow(m))
}

## ---- S6: GSEA --------------------------------------------------------------
log_msg("Supplementary Table 6: GSEA")
gf <- file.path(DE_DIR, "gsea_all_results.tsv")
if (file.exists(gf)) {
  g <- read.delim(gf, stringsAsFactors = FALSE)
  g$disease <- ifelse(g$disease %in% c("DE_SLE", "SLE"), "SLE",
               ifelse(g$disease %in% c("DE_SJS", "SjS"), "SjS", g$disease))
  names(g)[names(g) == "pval"] <- "pvalue"
  g <- g[, c("cell_type", "disease", "pathway", "NES", "pvalue", "padj", "leadingEdge"),
         drop = FALSE]
  g <- g[order(g$cell_type, g$disease, g$padj), , drop = FALSE]
  f5 <- file.path(OUTDIR, "supp_table_gsea.tsv")
  write.table(g, f5, sep = "\t", quote = FALSE, row.names = FALSE)
  note(f5, nrow(g))
} else {
  log_msg("  SKIPPED: ", gf, " not found — run 07_pathways.R")
}

## ---- S7: SHAP --------------------------------------------------------------
# Per-class mean |SHAP| (HC / SLE / SjS) alongside the class-averaged value that
# defines the ordering in Fig 4 and the text. No fold SDs: per-class fold SDs
# are large by construction for the 14-donor SjS class and are not reported.
log_msg("Supplementary Table 7: SHAP importance")
sh <- file.path(CLF_DIR, "shap_glmnet")
rd <- function(f) if (file.exists(file.path(sh, f))) read.delim(file.path(sh, f), stringsAsFactors = FALSE) else NULL
ct  <- rd("shap_importance_per_celltype.tsv")
pw  <- rd("shap_importance_per_pathway.tsv")
ft  <- rd("shap_importance_per_feature.tsv")
ctc <- rd("shap_importance_per_celltype_per_class.tsv")
pwc <- rd("shap_importance_per_pathway_per_class.tsv")
ftc <- rd("shap_importance_per_feature_per_class.tsv")

if (is.null(ct) || is.null(pw) || is.null(ft) || is.null(ctc) || is.null(pwc) || is.null(ftc)) {
  log_msg("  SKIPPED: SHAP outputs missing in ", sh, " — run 09_shap_importance.R")
} else {
  unsp <- function(x) gsub("_", " ", x)
  classes   <- c("healthy", "sle", "sjs")
  class_col <- c(healthy = "mean_abs_shap_HC", sle = "mean_abs_shap_SLE", sjs = "mean_abs_shap_SjS")
  # long (unit, class, value) -> wide, one column per class, rows keyed on `key`
  widen <- function(long, key, value) {
    keys <- unique(long[[key]])
    w <- data.frame(keys, stringsAsFactors = FALSE); colnames(w) <- key
    for (k in classes) {
      sub <- long[long$class == k, ]
      w[[class_col[[k]]]] <- sub[[value]][match(keys, sub[[key]])]
    }
    w
  }
  ct_w <- widen(ctc, "celltype", "mean_mean_abs_shap")
  pw_w <- widen(pwc, "pathway",  "mean_mean_abs_shap")
  ft_w <- widen(ftc, "feature",  "mean_abs_shap")

  s_ct <- data.frame(level = "cell_type",
                     cell_type = unsp(ct$celltype),
                     pathway = NA_character_,
                     mean_abs_shap = ct$mean_mean_abs_shap,
                     ct_w[match(ct$celltype, ct_w$celltype), class_col],
                     stringsAsFactors = FALSE)
  s_pw <- data.frame(level = "pathway",
                     cell_type = NA_character_,
                     pathway = pw$pathway,
                     mean_abs_shap = pw$mean_mean_abs_shap,
                     pw_w[match(pw$pathway, pw_w$pathway), class_col],
                     stringsAsFactors = FALSE)
  s_ft <- data.frame(level = "feature",
                     cell_type = unsp(ft$celltype),
                     pathway = ft$pathway,
                     mean_abs_shap = ft$mean_abs_shap,
                     ft_w[match(ft$feature, ft_w$feature), class_col],
                     stringsAsFactors = FALSE)
  stopifnot(!anyNA(s_ct[, class_col]), !anyNA(s_pw[, class_col]), !anyNA(s_ft[, class_col]))

  rank_within <- function(d) { d <- d[order(-d$mean_abs_shap), ]; d$rank <- seq_len(nrow(d)); d }
  s <- rbind(rank_within(s_ct), rank_within(s_pw), rank_within(s_ft))
  rownames(s) <- NULL
  f6 <- file.path(OUTDIR, "supp_table_shap.tsv")
  write.table(s, f6, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  note(f6, nrow(s))
}

log_msg("Done. ", length(written), " supplementary tables written to ", OUTDIR)
