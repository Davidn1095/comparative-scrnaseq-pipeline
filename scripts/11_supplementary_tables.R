#!/usr/bin/env Rscript
################################################################################
# 11_supplementary_tables.R — Regenerate Supplementary Tables 1-9
#
# Purpose: assemble every supplementary table from current pipeline outputs.
# These six files previously had no generator in the repository and had drifted
# 60-76 days behind the results they summarise.
#
# Schemas below were read off the existing files and are reproduced exactly, so
# regenerating is a drop-in replacement. The one exception is S6, which gained
# three per-class columns and dropped the fold-SD column when the classifier
# figure (now Fig 6 and Supplementary Fig 1) moved to per-class SHAP.
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
#   S7 SjD robustness   <- gsea_all_results.tsv (07) and gsea_sjs_gse253568.tsv (07b):
#                          per cell type and Hallmark pathway, NES and adjusted P for SLE,
#                          pooled SjD and SjD within GSE253568 alone
#   S8 Toro GO recovery <- results/manuscript1_figures/toro_comparison/toro_go_pathway_recovery.tsv
#                          (06b): per GO term reported as enriched (P < 0.01) by Toro-Dominguez
#                          et al. 2014, genes testable and recovered in SLE, SjD and both
#   S9 SHAP             <- results/donor_classifiers_<method>/shap_glmnet/ (script 09):
#                          class-averaged mean |SHAP| (ordering) + per-class
#                          HC / SLE / SjS columns; no fold SDs
#
# Numbering follows first mention in main.tex, counting figure captions:
# the discordant table is first cited in Methods (Differential expression),
# after S3 and S4 and before the GSEA table, so it is S5. The robustness table is first
# cited in Methods (Pathway enrichment analysis), after the GSEA table and before the SHAP
# table (Methods, Donor-level classification), so it is S7. The Toro-Dominguez table is first
# cited in Methods (External signature comparison), before the SHAP table, so it is S8 and
# SHAP is S9.
#
# Output: manuscript/supplementary/, one file per table, named by its number
# (Supplementary_Table_<n>_<content>.tsv[.gz]).
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

# ---- Display labels ----
# The manuscript names the disease Sjogren's disease (SjD). Input files,
# directories, filenames and internal keys keep "SjS"/"sjs" so upstream outputs
# still resolve; this one mapping relabels headers and text cells on output.
# DISPLAY_LABELS is matched as a substring ("HC_vs_SjS", "mean_abs_shap_SjS");
# DISPLAY_KEYS as a whole cell (the lowercase group keys of Supplementary Table 2).
# Cell types take their display names from CELLTYPE_DISPLAY (00_config.R), the
# map the figures use, matched as a whole cell.
DISPLAY_LABELS <- c(SjS = "SjD")
DISPLAY_KEYS   <- c(sjs = "sjd")
relabel_df <- function(d) {
  rl <- function(x) {
    for (k in names(DISPLAY_LABELS)) x <- gsub(k, DISPLAY_LABELS[[k]], x, fixed = TRUE)
    hit <- x %in% names(DISPLAY_KEYS); x[hit] <- DISPLAY_KEYS[x[hit]]
    display_celltype(x)
  }
  names(d) <- rl(names(d))
  for (j in seq_along(d)) if (is.character(d[[j]]) || is.factor(d[[j]])) d[[j]] <- rl(as.character(d[[j]]))
  d
}

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
f1 <- file.path(OUTDIR, "Supplementary_Table_1_lineage_mapping.tsv")
write.table(relabel_df(lm), f1, sep = "\t", quote = FALSE, row.names = FALSE)
note(f1, nrow(lm))
if (nrow(lm) != length(COMMON_CELL_TYPES)) {
  log_msg("  WARNING: expected ", length(COMMON_CELL_TYPES), " cell types, got ", nrow(lm))
}

## ---- S2: composition -------------------------------------------------------
log_msg("Supplementary Table 2: composition (propeller)")
psum <- file.path(ATLAS_ROOT, "results", "composition", "propeller_summary.tsv")
if (file.exists(psum)) {
  d2 <- read.delim(psum, stringsAsFactors = FALSE)
  f2 <- file.path(OUTDIR, "Supplementary_Table_2_composition.tsv")
  write.table(relabel_df(d2), f2, sep = "\t", quote = FALSE, row.names = FALSE)
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
for (spec in list(c("sle", "Supplementary_Table_3_DE_SLE.tsv.gz", "3"),
                  c("sjs", "Supplementary_Table_4_DE_SjD.tsv.gz", "4"))) {
  log_msg("Supplementary Table ", spec[3], ": full DE results (", toupper(spec[1]), ")")
  d <- collect_de(spec[1])
  if (is.null(d)) { log_msg("  SKIPPED: no DE tables in de_", spec[1]); next }
  f <- file.path(OUTDIR, spec[2])
  con <- gzfile(f, "w")
  write.table(relabel_df(d), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  note(f, nrow(d))
}

## ---- S5: discordant DEGs ---------------------------------------------------
# Genes called significant in both diseases within the same cell type but with
# opposite sign. Cited from the Fig 3 caption, where the count alone is given.
log_msg("Supplementary Table 5: discordant DEGs")
sle <- collect_de("sle"); sjs <- collect_de("sjs")
if (is.null(sle) || is.null(sjs)) {
  log_msg("  SKIPPED: DE tables missing - run 06_differential_analysis.R")
} else {
  keep <- function(d) d[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1 &
                        !is_excluded_gene(d$gene),
                        c("cell_type", "gene", "log2FoldChange", "padj")]
  a <- keep(sle); b <- keep(sjs)
  names(a)[3:4] <- c("log2FC_SLE", "padj_SLE")
  names(b)[3:4] <- c("log2FC_SjS", "padj_SjS")
  m <- merge(a, b, by = c("cell_type", "gene"))
  m <- m[sign(m$log2FC_SLE) != sign(m$log2FC_SjS), ]
  m <- m[order(m$cell_type, m$gene), ]
  m$log2FC_SLE <- round(m$log2FC_SLE, 3); m$log2FC_SjS <- round(m$log2FC_SjS, 3)
  m$padj_SLE <- signif(m$padj_SLE, 3);    m$padj_SjS <- signif(m$padj_SjS, 3)
  f7 <- file.path(OUTDIR, "Supplementary_Table_5_discordant_DEGs.tsv")
  write.table(relabel_df(m), f7, sep = "\t", quote = FALSE, row.names = FALSE)
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
  f5 <- file.path(OUTDIR, "Supplementary_Table_6_GSEA.tsv")
  write.table(relabel_df(g), f5, sep = "\t", quote = FALSE, row.names = FALSE)
  note(f5, nrow(g))
} else {
  log_msg("  SKIPPED: ", gf, " not found — run 07_pathways.R")
}

## ---- S7: SjD robustness (GSE253568 alone) ----------------------------------
# NES and adjusted P per cell type and Hallmark pathway for SLE (full cohort), SjD
# pooled over both SjD datasets and SjD within GSE253568, whose patients and controls
# carry the same dissociation-gene load (07b). NA where a cell type was not tested.
log_msg("Supplementary Table 7: SjD robustness, GSE253568 alone")
gf  <- file.path(DE_DIR, "gsea_all_results.tsv")
gf1 <- file.path(DE_DIR, "gsea_sjs_gse253568.tsv")
if (file.exists(gf) && file.exists(gf1)) {
  g  <- read.delim(gf, stringsAsFactors = FALSE)
  g1 <- read.delim(gf1, stringsAsFactors = FALSE)
  # 07b writes "Non-switched", "Non-Vd2" and "Th1/Th17" where 07 writes spaces; join on one form.
  key <- function(d) paste(gsub("[-/]", " ", d$cell_type), d$pathway, sep = "\r")
  sle <- g[g$disease == "SLE", ]; sjs <- g[g$disease == "SjS", ]
  keys <- unique(c(key(sle), key(sjs)))
  s7 <- data.frame(cell_type = sub("\r.*", "", keys), pathway = sub(".*\r", "", keys), stringsAsFactors = FALSE)
  pick <- function(d, col) d[[col]][match(keys, key(d))]
  s7$NES_SLE <- pick(sle, "NES");              s7$padj_SLE <- pick(sle, "padj")
  s7$NES_SjS_pooled <- pick(sjs, "NES");       s7$padj_SjS_pooled <- pick(sjs, "padj")
  s7$NES_SjS_GSE253568 <- pick(g1, "NES");     s7$padj_SjS_GSE253568 <- pick(g1, "padj")
  s7 <- s7[order(s7$cell_type, s7$pathway), ]
  f7 <- file.path(OUTDIR, "Supplementary_Table_7_SjD_robustness.tsv")
  write.table(relabel_df(s7), f7, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  note(f7, nrow(s7))
} else {
  log_msg("  SKIPPED: ", if (!file.exists(gf)) gf else gf1, " not found — run 07_pathways.R and 07b_sjd_single_dataset.R")
}

## ---- S8: Toro-Dominguez GO term recovery -----------------------------------
# Per GO biological-process term reported as enriched (P < 0.01) by Toro-Dominguez et al. 2014:
# its genes, those testable in both diseases, and those recovered as DEGs in the term's
# direction in SLE, in SjD and in both (06b).
log_msg("Supplementary Table 8: Toro-Dominguez GO term recovery")
tf <- file.path(ATLAS_ROOT, "results", "manuscript1_figures", "toro_comparison", "toro_go_pathway_recovery.tsv")
if (file.exists(tf)) {
  t8 <- read.delim(tf, stringsAsFactors = FALSE)
  s8 <- data.frame(direction = ifelse(t8$direction == "up", "over-expressed", "under-expressed"),
                   go_id = t8$go_id, term = sub(" \\(BP\\)$", "", t8$term),
                   p_value_original = t8$p_value, n_genes = t8$n_genes, n_testable = t8$n_testable,
                   n_recovered_SLE = t8$n_recovered_sle, n_recovered_SjD = t8$n_recovered_sjs,
                   n_recovered_both = t8$n_recovered_both, genes_recovered_both = t8$recovered_both,
                   stringsAsFactors = FALSE)
  f8 <- file.path(OUTDIR, "Supplementary_Table_8_Toro_GO_recovery.tsv")
  write.table(relabel_df(s8), f8, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  note(f8, nrow(s8))
} else {
  log_msg("  SKIPPED: ", tf, " not found — run 06b_toro_comparison.R")
}

## ---- S9: SHAP --------------------------------------------------------------
# Per-class mean |SHAP| (HC / SLE / SjS) alongside the class-averaged value that
# defines the ordering in Fig 6, Supplementary Fig 1 and the text. No fold SDs: per-class fold SDs
# are large by construction for the 14-donor SjS class and are not reported.
log_msg("Supplementary Table 9: SHAP importance")
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
  f9 <- file.path(OUTDIR, "Supplementary_Table_9_SHAP.tsv")
  write.table(relabel_df(s), f9, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  note(f9, nrow(s))
}

log_msg("Done. ", length(written), " supplementary tables written to ", OUTDIR)
