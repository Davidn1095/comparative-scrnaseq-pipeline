#!/usr/bin/env Rscript
################################################################################
# Script 06: Pseudobulk DESeq2 Differential Expression
#
# Purpose: Run pseudobulk DESeq2 differential expression comparing disease vs
#          healthy donors per cell type.
#
# Input:
#   - Pseudobulk counts from script 03
#
# Output:
#   - DESeq2 results: results/plots/differential_expression/de_{disease}/deseq2_{cell_type}.tsv
#   - Skipped cell-type audit: de_skipped_celltypes.tsv
#
# Usage:
#   Rscript 06_differential_analysis.R
#
# Environment variables:
#   ATLAS_ROOT    - Atlas root directory
#   DISEASE       - Disease to analyse (SLE, SjS) or "all" (default: all)
#   N_CORES       - Number of parallel cores
#   PADJ_THRESH   - Adjusted p-value threshold (default: 0.05)
#   LFC_THRESH    - Log2 fold change threshold (default: 0.25)
#   MIN_DONORS          - Minimum donors per condition after floor (default: 5)
#   MIN_CELLS           - Minimum total cells per cell type (default: 100)
#   MIN_CELLS_CLASS     - Minimum cells per condition (default: 10)
#   MIN_CELLS_PER_DONOR - Minimum cells per (donor, cell type) (default: 20)
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
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
})

DISEASE_FILTER <- Sys.getenv("DISEASE", "all")
# MIN_DONORS, MIN_CELLS, MIN_CELLS_CLASS and MIN_CELLS_PER_DONOR come
# from 00_config.R (env-overridable). The unified thresholds match the
# singleDeep prep filters in 09_prepare_singledeep.R.

cat("\n")
cat("=======================================================================\n")
cat("         PSEUDOBULK DESeq2 DIFFERENTIAL EXPRESSION\n")
cat("=======================================================================\n")
cat("\n")
cat("Configuration:\n")
cat("  Disease filter:", DISEASE_FILTER, "\n")
cat("  Cores:", N_CORES, "\n")
cat("  P-adj threshold:", PADJ_THRESH, "\n")
cat("  LFC threshold:", LFC_THRESH, "\n")
cat("\n")

################################################################################
# SHARED: Condition standardization
################################################################################

standardize_conditions <- function(conditions) {
  conditions_std <- tolower(conditions)
  case_when(
    conditions_std %in% c("healthy", "hd", "normal", "control") ~ "Healthy",
    conditions_std %in% c("sle", "lupus") ~ "SLE",
    conditions_std %in% c("sjs", "sjd", "sjogren") ~ "SjS",
    conditions_std %in% c("ra", "rheumatoid") ~ "RA",
    TRUE ~ NA_character_
  )
}

################################################################################
# PART A: DESeq2 PSEUDOBULK DIFFERENTIAL EXPRESSION
################################################################################

run_deseq2_analysis <- function(disease, celltype_ncells = NULL,
                                 celltype_donor_ncells = NULL) {
  cat("\n--- DESeq2 Analysis for", disease, "---\n")

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    cat("  DESeq2 package not available, skipping\n")
    return(NULL)
  }
  suppressPackageStartupMessages(library(DESeq2))

  DISEASE_DIR <- file.path(ATLAS_ROOT, "data", disease)
  DESEQ2_OUTDIR <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression", paste0("de_", tolower(disease)))

  if (!dir.exists(DISEASE_DIR)) {
    cat("  Disease directory not found:", DISEASE_DIR, "\n")
    return(NULL)
  }

  dir.create(DESEQ2_OUTDIR, recursive = TRUE, showWarnings = FALSE)

  # Load pseudobulk data from all accessions
  acc_dirs <- list.dirs(DISEASE_DIR, recursive = FALSE, full.names = TRUE)

  all_counts <- list()
  all_meta <- list()
  skipped_gap <- character(0)      # looks like a dataset, but pseudobulk is missing
  skipped_nondata <- character(0)  # not a dataset directory at all

  for (d in acc_dirs) {
    acc <- basename(d)
    counts_file <- file.path(d, "derived", "pseudobulk", "counts", "pb_counts.rds")
    meta_file <- file.path(d, "derived", "pseudobulk", "counts", "pb_meta.tsv")

    # A cohort must never drop out of the meta-analysis in silence. This used to be a
    # bare `next`: a dataset whose pseudobulk had failed to build simply vanished from
    # the DE, the run reported success, and the only trace was a smaller "Loaded N
    # datasets" that nobody had a reference value for. Distinguish a real pipeline gap
    # - an annotated object exists but its pseudobulk does not - from a directory that
    # was never a dataset, and shout about the former.
    if (!file.exists(counts_file) || !file.exists(meta_file)) {
      has_obj <- length(list.files(file.path(d, "derived", "objects"),
                                   pattern = "\\.rds$")) > 0
      if (has_obj) {
        skipped_gap <- c(skipped_gap, acc)
        cat("  !! WARNING: ", acc, " has an annotated object but NO pseudobulk -\n",
            "  !!   missing: ",
            if (!file.exists(counts_file)) basename(counts_file) else basename(meta_file),
            "\n  !!   EXCLUDED from the ", disease, " meta-analysis.\n", sep = "")
      } else {
        skipped_nondata <- c(skipped_nondata, acc)
      }
      next
    }

    counts <- readRDS(counts_file)
    meta <- read.delim(meta_file, stringsAsFactors = FALSE)

    if (!"dataset_id" %in% colnames(meta)) {
      meta$dataset_id <- acc
    }

    meta$pb_id <- paste0(acc, ":", meta$pb_id)
    colnames(counts) <- paste0(acc, ":", colnames(counts))

    all_counts[[acc]] <- counts
    all_meta[[acc]] <- meta
  }

  if (length(all_counts) == 0) {
    cat("  No pseudobulk data found for", disease, "\n")
    return(NULL)
  }

  cat("  Loaded", length(all_counts), "datasets:",
      paste(names(all_counts), collapse = ", "), "\n")
  if (length(skipped_nondata))
    cat("  (ignored", length(skipped_nondata), "non-dataset directories)\n")
  if (length(skipped_gap)) {
    cat("\n", strrep("!", 74), "\n", sep = "")
    cat("!! ", length(skipped_gap), " dataset(s) EXCLUDED from ", disease,
        " for want of pseudobulk:\n", sep = "")
    for (a in skipped_gap) cat("!!   ", a, "\n", sep = "")
    cat("!! These have annotated objects, so this is a pipeline gap, not a design\n")
    cat("!! choice. Run 03_pseudobulk.R for them, or the meta-analysis below is\n")
    cat("!! reporting on fewer cohorts than the study contains.\n")
    cat(strrep("!", 74), "\n\n", sep = "")
  }

  # Combine data
  common_genes <- Reduce(intersect, lapply(all_counts, rownames))
  all_counts <- lapply(all_counts, function(m) m[common_genes, , drop = FALSE])
  combined_counts <- do.call(cbind, all_counts)
  combined_meta <- do.call(rbind, all_meta)
  rownames(combined_meta) <- combined_meta$pb_id

  # Preserve the upstream clinical-condition strings ("healthy", "sle", "sjs",
  # "ra") under a clearer name; reuse the `condition` slot for the binary
  # case/control factor that DESeq2 contrasts on.
  combined_meta$clinical_condition <- combined_meta$condition

  # Standardize conditions
  disease_l <- tolower(disease)
  ctrl_labels <- c("normal", "hd", "healthy", "control")
  case_labels <- c(disease_l, "sle", "sjs", "ra", "case", "disease")

  combined_meta$condition <- NA_character_
  combined_meta$condition[tolower(combined_meta$clinical_condition) %in% ctrl_labels] <- "control"
  combined_meta$condition[tolower(combined_meta$clinical_condition) %in% case_labels & is.na(combined_meta$condition)] <- "case"

  combined_meta <- combined_meta[!is.na(combined_meta$condition), ]
  combined_counts <- combined_counts[, combined_meta$pb_id, drop = FALSE]

  cat("  Case:", sum(combined_meta$condition == "case"),
      "| Control:", sum(combined_meta$condition == "control"), "\n")

  # Run DESeq2 per cell type — restricted to the canonical 25 (config).
  cell_types <- sort(unique(combined_meta$cell_type))
  cell_types <- cell_types[is_common_celltype(cell_types)]

  run_deseq2_ct <- function(ct) {
    idx <- combined_meta$cell_type == ct
    ct_meta <- combined_meta[idx, ]

    # Per-donor cell-count floor (mirrors 09_prepare_singledeep.R): drop
    # donors whose cell count for this cell type is < MIN_CELLS_PER_DONOR.
    # Stabilises pseudobulk by removing donors whose count vectors are
    # dominated by sampling noise. Done before the donor count check so
    # MIN_DONORS is applied to the surviving donor count.
    if (!is.null(celltype_donor_ncells)) {
      ct_donor <- celltype_donor_ncells[celltype_donor_ncells$cell_type == ct, ]
      donors_above_floor <- ct_donor$donor[ct_donor$n_cells >= MIN_CELLS_PER_DONOR]
      ct_meta <- ct_meta[ct_meta$donor_id %in% donors_above_floor, ]
    }
    ct_counts <- combined_counts[, ct_meta$pb_id, drop = FALSE]

    n_case <- sum(ct_meta$condition == "case")
    n_ctrl <- sum(ct_meta$condition == "control")

    if (n_case < MIN_DONORS || n_ctrl < MIN_DONORS) {
      return(list(ct = ct, status = "skipped", n_case = n_case, n_ctrl = n_ctrl,
                  msg = paste0("too few donors after MIN_CELLS_PER_DONOR=",
                               MIN_CELLS_PER_DONOR, " floor (case: ",
                               n_case, ", ctrl: ", n_ctrl, ")")))
    }

    # Cell-count thresholds (belt-and-braces; trivially satisfied when
    # MIN_DONORS x MIN_CELLS_PER_DONOR >= MIN_CELLS, which holds for
    # the unified set 5 x 20 = 100).
    if (!is.null(celltype_donor_ncells)) {
      ct_donor <- celltype_donor_ncells[celltype_donor_ncells$cell_type == ct, ]
      n_cells_case <- sum(ct_donor$n_cells[ct_donor$donor %in% ct_meta$donor_id[ct_meta$condition == "case"]])
      n_cells_ctrl <- sum(ct_donor$n_cells[ct_donor$donor %in% ct_meta$donor_id[ct_meta$condition == "control"]])
      if ((n_cells_case + n_cells_ctrl) < MIN_CELLS) {
        return(list(ct = ct, status = "skipped", n_case = n_case, n_ctrl = n_ctrl,
                    msg = paste0("too few total cells (", n_cells_case + n_cells_ctrl, ")")))
      }
      if (n_cells_case < MIN_CELLS_CLASS || n_cells_ctrl < MIN_CELLS_CLASS) {
        return(list(ct = ct, status = "skipped", n_case = n_case, n_ctrl = n_ctrl,
                    msg = paste0("too few cells per class (case: ",
                                 n_cells_case, ", ctrl: ", n_cells_ctrl, ")")))
      }
    }

    ct_counts <- as.matrix(ct_counts)
    storage.mode(ct_counts) <- "integer"

    keep_genes <- rowSums(ct_counts) > 0
    ct_counts <- ct_counts[keep_genes, , drop = FALSE]

    ct_meta$condition <- factor(ct_meta$condition, levels = c("control", "case"))
    ct_meta$dataset_id <- factor(ct_meta$dataset_id)

    n_datasets <- nlevels(ct_meta$dataset_id)
    use_batch <- FALSE

    if (n_datasets > 1) {
      tab <- table(ct_meta$dataset_id, ct_meta$condition)
      confounded <- any(colSums(tab > 0) < 2)
      if (!confounded) {
        use_batch <- TRUE
      }
    }

    design_formula <- if (use_batch) ~ dataset_id + condition else ~ condition

    dds <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = ct_counts,
        colData = ct_meta,
        design = design_formula
      )
      dds <- DESeq(dds, quiet = TRUE)
      dds
    }, error = function(e) {
      return(list(ct = ct, status = "error", msg = conditionMessage(e)))
    })

    if (is.list(dds) && !is(dds, "DESeqDataSet")) return(dds)

    res <- tryCatch(
      results(dds, contrast = c("condition", "case", "control")),
      error = function(e) list(ct = ct, status = "error", msg = conditionMessage(e))
    )

    if (is.list(res) && !is(res, "DESeqResults")) return(res)

    df <- data.frame(
      gene = rownames(res),
      log2FoldChange = res$log2FoldChange,
      pvalue = res$pvalue,
      padj = res$padj,
      baseMean = res$baseMean,
      lfcSE = res$lfcSE,
      stat = res$stat,
      stringsAsFactors = FALSE
    )

    df <- df[!grepl(EXCLUDED_GENE_REGEX, df$gene, ignore.case = TRUE), ]
    df <- df[!is.na(df$pvalue), ]
    df <- df[order(df$padj, df$pvalue), ]

    tag <- sanitize_name(ct)
    out_file <- file.path(DESEQ2_OUTDIR, paste0("deseq2_", tag, ".tsv"))
    write.table(df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

    n_sig <- sum(df$padj < PADJ_THRESH & abs(df$log2FoldChange) > LFC_THRESH, na.rm = TRUE)
    list(ct = ct, status = "done", n_case = n_case, n_ctrl = n_ctrl,
         design = deparse(design_formula), n_sig = n_sig)
  }

  results_list <- parallel::mclapply(cell_types, run_deseq2_ct, mc.cores = N_CORES)

  for (r in results_list) {
    if (r$status == "skipped") {
      cat("    ", r$ct, "- SKIPPED:", r$msg, "\n")
    } else if (r$status == "error") {
      cat("    ", r$ct, "- ERROR:", r$msg, "\n")
    } else {
      cat("    ", r$ct, "- Sig DEGs:", r$n_sig, "\n")
    }
  }

  return(results_list)
}


################################################################################
# MAIN EXECUTION
################################################################################

# Determine diseases to analyze
if (DISEASE_FILTER == "all") {
  diseases_to_analyze <- DISEASES_ACTIVE
} else {
  diseases_to_analyze <- DISEASE_FILTER
}

# Load atlas to look up per-cell-type cell counts for the MIN_CELLS check
cat("\nLoading integrated atlas...\n")
atlas <- tryCatch({
  load_atlas()
}, error = function(e) {
  cat("  Error loading atlas:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(atlas)) {
  cat("  Cells:", ncol(atlas), "\n")
}

celltype_ncells <- NULL
celltype_donor_ncells <- NULL
if (!is.null(atlas)) {
  meta <- atlas@meta.data
  meta$clinical_condition <- meta$condition
  meta$condition_std <- standardize_conditions(meta$clinical_condition)
  meta_clean <- meta[!is.na(meta$condition_std), ]
  celltype_ncells <- meta_clean |>
    (\(d) aggregate(list(n_cells = rep(1L, nrow(d))),
                    by = list(cell_type = d$cell_type, condition_std = d$condition_std),
                    FUN = length))()
  # Per-donor, per-cell-type cell counts for the MIN_CELLS_PER_DONOR floor.
  celltype_donor_ncells <- meta_clean |>
    (\(d) aggregate(list(n_cells = rep(1L, nrow(d))),
                    by = list(cell_type = d$cell_type, donor = d$donor,
                              condition_std = d$condition_std),
                    FUN = length))()
}

# Run DESeq2 for each disease
deseq2_results_all <- list()
cat("\n")
cat("=======================================================================\n")
cat("  DESeq2 PSEUDOBULK DIFFERENTIAL EXPRESSION\n")
cat("=======================================================================\n")

for (disease in diseases_to_analyze) {
  result <- tryCatch({
    run_deseq2_analysis(disease, celltype_ncells, celltype_donor_ncells)
  }, error = function(e) {
    cat("  Error in DESeq2 for", disease, ":", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(result)) deseq2_results_all[[disease]] <- result
}

# Write skipped cell-type summary for supplementary materials
skip_records <- list()
if (length(deseq2_results_all) > 0) {
  for (disease in names(deseq2_results_all)) {
    for (r in deseq2_results_all[[disease]]) {
      if (!is.null(r$status) && r$status == "skipped") {
        skip_records[[paste0("DESeq2_", disease, "_", r$ct)]] <- data.frame(
          disease = disease, cell_type = r$ct,
          n_donors_case = r$n_case, n_donors_ctrl = r$n_ctrl,
          skip_reason = r$msg, stringsAsFactors = FALSE
        )
      }
    }
  }
}
# ALWAYS write this file, even when nothing was skipped. It is the record of which
# cell types the floors excluded, and 10_manuscript_figures.R stops without it. If it
# were written only when non-empty, "no skips" and "06 never ran" would be
# indistinguishable on disk - and a stale copy from an earlier run would be read as
# though it described this one. An empty file with headers says "06 ran, nothing was
# skipped"; an absent file says "06 did not run".
skip_out <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression", "de_skipped_summary.tsv")
dir.create(dirname(skip_out), recursive = TRUE, showWarnings = FALSE)
skip_df <- if (length(skip_records) > 0) do.call(rbind, skip_records) else
  data.frame(disease = character(0), cell_type = character(0),
             n_donors_case = integer(0), n_donors_ctrl = integer(0),
             skip_reason = character(0), stringsAsFactors = FALSE)
write.table(skip_df, skip_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSkipped cell types written to:", skip_out, "\n")
cat("Total skipped:", nrow(skip_df), "cell-type/disease combinations\n")

cat("\n")
cat("=======================================================================\n")
cat("              DIFFERENTIAL ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("\n")
cat("Output: results/plots/differential_expression/de_{disease}/deseq2_*.tsv\n")
cat("\n")
