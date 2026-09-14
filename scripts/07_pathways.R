#!/usr/bin/env Rscript
################################################################################
# Script 07: Signatures and Pathway Analysis
#
# Purpose: Generate disease signatures and run GSEA pathway enrichment
#          on the per-cell-type DESeq2 outputs.
#
# Input:
#   - DESeq2 results from script 06 (de_*/deseq2_*.tsv)
#
# Output:
#   - Disease signatures (gene lists)
#   - GSEA enrichment results
#   - Condition similarity matrices
#
# Usage:
#   Rscript 07_pathways.R
#
# Environment variables:
#   ATLAS_ROOT    - Atlas root directory
#   PADJ_THRESH   - Adjusted p-value threshold (default: 0.05)
#   LFC_THRESH    - Log2 fold change threshold (default: 1.0)
################################################################################

# Source shared configuration and utilities
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
  library(Matrix)
})

# Output directories (all under differential_expression)
DE_DIR <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression")
OUTDIR <- file.path(DE_DIR, "signatures")
GSEA_DIR <- DE_DIR
SIMILARITY_DIR <- file.path(DE_DIR, "similarity")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SIMILARITY_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Signatures and Pathway Analysis ===\n")
cat("Output directories:\n")
cat("  Signatures:", OUTDIR, "\n")
cat("  GSEA:", GSEA_DIR, "\n")
cat("  Similarity:", SIMILARITY_DIR, "\n\n")

################################################################################
# 1. Collect DESeq2 results from all datasets
################################################################################

cat("Collecting DE results...\n")

# Helper: canonical per-cell-type DE TSVs only.
# Skips legacy de_de_*/ duplicates and the *_annotated.tsv copies that the
# annotation step (below) writes alongside the canonical files; reading those
# back in would propagate "_annotated" into cell-type labels and fan out into
# duplicate similarity_*.tsv / signature_*.tsv files downstream.
canonical_de_files <- function() {
  de_root <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression")
  dirs <- list.dirs(de_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^de_(sle|sjs)$", basename(dirs))]
  unlist(lapply(dirs, function(d) {
    files <- list.files(d, pattern = "^deseq2_.*\\.tsv$", full.names = TRUE)
    files[!grepl("_annotated\\.tsv$", basename(files))]
  }), use.names = FALSE)
}

deseq2_files <- canonical_de_files()
if (length(deseq2_files) == 0) {
  stop("No DE results found. Run 04_differential_expression.R first.\n")
}

# Restrict to canonical 25 cell types. DE may write extras if stale outputs
# linger from prior runs; this guarantees signatures/GSEA/similarity stay
# aligned to the COMMON_CELL_TYPES set.
file_celltype <- gsub("_", " ", sub("\\.tsv$", "", sub("^deseq2_", "", basename(deseq2_files))))
deseq2_files <- deseq2_files[is_common_celltype(file_celltype)]

cat("  Found", length(deseq2_files), "DE result files (canonical only)\n")

# Load and combine results
deseq2_list <- lapply(deseq2_files, function(f) {
  # Extract disease from path (de_{disease}/deseq2_*.tsv)
  raw_label <- basename(dirname(f))
  disease_key <- sub("^de_", "", raw_label)

  # Extract cell type from filename
  ct <- sub("deseq2_", "", sub("\\.tsv$", "", basename(f)))
  ct <- gsub("_", " ", ct)

  res <- read.delim(f, stringsAsFactors = FALSE)
  # Map folder names to standard disease labels
  disease_map <- c("sle" = "SLE", "sjs" = "SjS")
  res$disease <- ifelse(disease_key %in% names(disease_map), disease_map[disease_key], toupper(disease_key))
  res$cell_type <- ct

  # Compute stat for GSEA ranking: signed -log10(pvalue)
  if (!"stat" %in% colnames(res) && "pvalue" %in% colnames(res)) {
    res$stat <- -log10(pmax(res$pvalue, 1e-300)) * sign(res$log2FoldChange)
  }

  res
})

deseq2_all <- bind_rows(deseq2_list)
cat("  Total DE tests:", nrow(deseq2_all), "\n")

################################################################################
# 2. Generate disease signatures per cell type
################################################################################

cat("\nGenerating disease signatures...\n")

# Filter significant genes
sig_genes <- deseq2_all %>%
  filter(
    !is.na(padj),
    padj < PADJ_THRESH,
    abs(log2FoldChange) > LFC_THRESH,
    !grepl(EXCLUDED_GENE_REGEX, gene, ignore.case = TRUE)
  )

cat("  Significant DE genes:", n_distinct(sig_genes$gene), "\n")

# Save signatures per disease per cell type — uncapped: writing every
# gene that passes the strict thresholds (padj < PADJ_THRESH,
# |log2FC| > LFC_THRESH, MT/RPL/RPS/HB excluded). Previously capped at
# TOP_N_GENES=50 per (disease, cell_type), but that ceiling clipped panel a
# of fig2 (DEG counts) to 50 for cell types with hundreds of true DEGs.
signatures <- sig_genes %>%
  arrange(disease, cell_type, padj)

# Save as combined file
write.table(
  signatures,
  file.path(OUTDIR, "disease_signatures_all.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Per-(disease, cell type) signature files and signature_summary.tsv were
# written here but had no consumer anywhere in the pipeline or manuscript;
# removed. disease_signatures_all.tsv above is the file script 10 reads.

cat("  Saved", nrow(signatures), "signature rows\n")

################################################################################
# 3. GSEA pathway enrichment
################################################################################

cat("\nRunning GSEA pathway enrichment...\n")

has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
has_fgsea <- requireNamespace("fgsea", quietly = TRUE)

# Fallback: use cached msigdb RDS if msigdbr is not installed
# MSigDB 2025.1 Hs cache in msigdbr's format; set MSIGDB_RDS to use another location.
cached_msigdb <- Sys.getenv("MSIGDB_RDS", file.path(tools::R_user_dir("msigdbr", which = "data"), "msigdb.2025.1.Hs.rds"))
use_cached <- !has_msigdbr && file.exists(cached_msigdb)

if (!has_fgsea || (!has_msigdbr && !use_cached)) {
  cat("  WARNING: msigdbr or fgsea not installed. Skipping GSEA.\n")
  cat("  Install with: BiocManager::install(c('msigdbr', 'fgsea'))\n")
} else {
  library(fgsea)

  if (use_cached) {
    cat("  Using cached msigdb: ", cached_msigdb, "\n")
    msig <- readRDS(cached_msigdb)
    # Cached RDS uses db_gene_symbol, gs_collection, gs_subcollection
    hallmark <- msig[msig$gs_collection == "H", ]
    hallmark_list <- split(hallmark$db_gene_symbol, hallmark$gs_name)
    reactome <- msig[msig$gs_collection == "C2" & msig$gs_subcollection == "CP:REACTOME", ]
  } else {
    library(msigdbr)
    hallmark <- msigdbr(species = "Homo sapiens", category = "H")
    hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
    reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  }

  all_pathways <- hallmark_list

  # Run GSEA for each disease/cell type combination
  gsea_results <- list()

  for (d in unique(deseq2_all$disease)) {
    for (ct in unique(deseq2_all$cell_type[deseq2_all$disease == d])) {

      de_subset <- deseq2_all %>%
        filter(disease == d, cell_type == ct, !is.na(stat))

      if (nrow(de_subset) < 100) next

      # Create ranked gene list
      ranks <- setNames(de_subset$stat, de_subset$gene)
      ranks <- sort(ranks, decreasing = TRUE)

      # Run fgsea
      set.seed(42)
      fgsea_res <- fgsea(
        pathways = all_pathways,
        stats = ranks,
        minSize = 15,
        maxSize = 500,
        nproc = N_CORES
      )

      fgsea_res$disease <- d
      fgsea_res$cell_type <- ct

      gsea_results[[paste(d, ct, sep = "__")]] <- fgsea_res
    }
  }

  # Combine and save
  if (length(gsea_results) > 0) {
    gsea_combined <- bind_rows(gsea_results)

    # Convert leadingEdge list column to semicolon-separated string
    if ("leadingEdge" %in% colnames(gsea_combined)) {
      gsea_combined$leadingEdge <- sapply(gsea_combined$leadingEdge,
                                          function(x) paste(x, collapse = ";"))
    }

    write.table(
      gsea_combined,
      file.path(GSEA_DIR, "gsea_all_results.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    # Save significant pathways
    gsea_sig <- gsea_combined %>%
      filter(padj < PADJ_THRESH) %>%
      arrange(disease, cell_type, padj)

    write.table(
      gsea_sig,
      file.path(GSEA_DIR, "gsea_significant.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    cat("  GSEA complete:", nrow(gsea_sig), "significant pathway enrichments\n")
  }
}

################################################################################
# 4. Condition similarity analysis
################################################################################

cat("\nCalculating condition similarity...\n")

# Same-sign Jaccard, aligned with the pathway-level definition used for Fig 3c:
# a gene counts as shared only when it is significant in both diseases AND moves
# in the same direction; the denominator is the union of genes significant in
# either. The magnitude criterion still differs between the two levels (genes
# require |log2FC| > LFC_THRESH, pathways have no |NES| threshold) and that is
# deliberate, so it is not reconciled here.
jaccard_same_sign <- function(g1, g2) {
  u <- union(names(g1), names(g2))
  if (length(u) == 0) return(0)
  shared <- intersect(names(g1), names(g2))
  same <- sum(sign(g1[shared]) == sign(g2[shared]))
  same / length(u)
}

# Calculate similarity for each cell type on FULL significant DEG lists
# (not top-50-clipped signatures, which would underestimate cross-disease overlap
# in cell types with >50 DEGs by clipping shared genes that rank outside the top 50).
# Every canonical cell type is emitted, including those with an empty union
# (Jaccard 0), so the gene- and pathway-level tables share one 25-cell-type basis.
cell_types <- sort(unique(deseq2_all$cell_type))
diseases <- sort(unique(deseq2_all$disease))

for (ct in cell_types) {

  # Named log2FC vectors per disease: names give set membership, sign gives direction
  ct_sigs <- lapply(setNames(diseases, diseases), function(d) {
    s <- sig_genes %>% filter(disease == d, cell_type == ct)
    setNames(s$log2FoldChange, s$gene)
  })

  conds <- names(ct_sigs)
  sim_mat <- matrix(0, nrow = length(conds), ncol = length(conds),
                    dimnames = list(conds, conds))

  for (i in seq_along(conds)) {
    for (j in seq_along(conds)) {
      sim_mat[i, j] <- if (i == j) {
        if (length(ct_sigs[[conds[i]]]) > 0) 1 else 0
      } else {
        jaccard_same_sign(ct_sigs[[conds[i]]], ct_sigs[[conds[j]]])
      }
    }
  }

  # Save
  ct_safe <- sanitize_name(ct)
  write.table(
    as.data.frame(sim_mat),
    file.path(SIMILARITY_DIR, paste0("similarity_", ct_safe, ".tsv")),
    sep = "\t", quote = FALSE
  )
}

cat("  Similarity matrices saved\n")

cat("\n=== Signatures and Pathway Analysis Complete ===\n")
cat("Output files:\n")
cat("  Signatures:\n")
cat("    - disease_signatures_all.tsv\n")
cat("  GSEA:\n")
cat("    - gsea_all_results.tsv\n")
cat("    - gsea_significant.tsv\n")
cat("  Similarity:\n")
cat("    - similarity_*.tsv (per cell type)\n")
