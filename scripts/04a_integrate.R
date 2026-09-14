#!/usr/bin/env Rscript
################################################################################
# Script 04: Integrate Atlas
#
# Purpose: Merge all annotated datasets into a unified atlas. No integration is
#          performed here: PCA/Harmony/UMAP were stripped as unused. The
#          batch-corrected embedding used by Fig 1 comes from scVI (04b).
#
# Input:
#   - Annotated Seurat objects from script 02 (*_qc_preproc_singler_monaco.rds)
#
# Output:
#   - Integrated Seurat object (atlas_integrated.rds)
#   - Integration summary statistics
#
# Usage:
#   Rscript 04a_integrate.R
#
# Environment variables:
#   ATLAS_ROOT - Atlas root directory
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
  library(SeuratObject)
  library(dplyr)
})

# Output paths
# ATLAS_OVERVIEW_DIR: atlas-level overview plots/tables (results/plots/atlas_overview)
ATLAS_OVERVIEW_DIR <- file.path(ATLAS_ROOT, "results", "plots", "atlas_overview")
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ATLAS_OVERVIEW_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg("=== Atlas Integration ===")
log_msg("Output directory: ", ATLAS_OVERVIEW_DIR)
log_msg("Diseases excluded: ", DISEASES_EXCLUDE)
log_msg("Conditions excluded: ", CONDITIONS_EXCLUDE)

################################################################################
# 1. Find and load all annotated objects
#    Active datasets (4): GSE157278, GSE253568 (SjS) +
#                         CXG_436154da-bcf1-4130-9c8b-120ff9a888f2, GSE162577 (SLE)
#    Anything under data/excluded/ is skipped by definition.
################################################################################

ACTIVE_DATASETS <- c(
  "CXG_436154da-bcf1-4130-9c8b-120ff9a888f2",
  "GSE162577",
  "GSE157278",
  "GSE253568"
)

log_msg("Searching for annotated Seurat objects...")

data_root <- file.path(ATLAS_ROOT, "data")
disease_dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)
disease_dirs <- disease_dirs[basename(disease_dirs) != "excluded"]

# Exclude specified diseases (e.g., RA)
if (DISEASES_EXCLUDE != "") {
  diseases_to_exclude <- strsplit(DISEASES_EXCLUDE, ",")[[1]]
  disease_dirs <- disease_dirs[!basename(disease_dirs) %in% diseases_to_exclude]
  log_msg("  Excluding diseases: ", paste(diseases_to_exclude, collapse = ", "))
}
dataset_dirs <- unlist(lapply(disease_dirs, function(d) {
  list.dirs(d, recursive = FALSE, full.names = TRUE)
}), use.names = FALSE)

# Sanity check: filesystem discovery must match the active-dataset whitelist.
found <- basename(dataset_dirs)
if (!setequal(found, ACTIVE_DATASETS)) {
  stop("Dataset discovery mismatch.\n  Expected: ",
       paste(ACTIVE_DATASETS, collapse = ", "),
       "\n  Found:    ", paste(found, collapse = ", "),
       "\n  Move stray datasets to data/excluded/ before integrating.")
}
log_msg("  Active datasets (", length(ACTIVE_DATASETS), "): ",
        paste(found, collapse = ", "))

# Load function for a per-dataset annotated Seurat object
read_obj <- function(d) {
  acc <- basename(d)
  disease <- basename(dirname(d))

  # Try different possible filenames
  rds_candidates <- c(
    file.path(d, "derived", "objects", paste0(acc, "_qc_preproc_singler_monaco.rds")),
    file.path(d, "derived", "seurat", "seurat_annotated.rds")
  )

  rds <- rds_candidates[file.exists(rds_candidates)][1]
  if (is.na(rds)) return(NULL)

  obj <- readRDS(rds)

  # Standardize metadata (always use folder-based disease label)
  obj$dataset_id <- acc
  obj$disease <- disease
  obj$study <- acc

  # Standardize cell type column
  if ("cell_type_pruned" %in% colnames(obj@meta.data)) {
    obj$cell_type <- obj$cell_type_pruned
  } else if ("cell_type_singler" %in% colnames(obj@meta.data)) {
    obj$cell_type <- obj$cell_type_singler
  }

  # Standardize donor column
  if (!"donor" %in% colnames(obj@meta.data)) {
    if ("donor_id" %in% colnames(obj@meta.data)) {
      obj$donor <- obj$donor_id
    } else if ("subject_id" %in% colnames(obj@meta.data)) {
      obj$donor <- obj$subject_id
    } else {
      obj$donor <- obj$sample_id
    }
  }

  log_msg("  Loading: ", disease, "/", acc, " (", ncol(obj), " cells)")

  obj
}

# Load all objects
objs <- lapply(dataset_dirs, read_obj)
objs <- objs[!vapply(objs, is.null, logical(1))]

if (length(objs) == 0) {
  stop("No annotated Seurat objects found.\n",
       "Please run scripts 01-02 first to generate annotated objects.\n")
}

log_msg("Loaded ", length(objs), " datasets")

################################################################################
# 2. Merge objects
################################################################################

log_msg("Merging datasets...")

if (length(objs) == 1) {
  seu <- objs[[1]]
} else {
  seu <- merge(objs[[1]], y = objs[-1])
}

log_msg("  Total cells: ", ncol(seu))
log_msg("  Total genes: ", nrow(seu))
log_msg("  Datasets: ", n_distinct(seu$dataset_id))
log_msg("  Diseases: ", paste(unique(seu$disease), collapse = ", "))

# Filter excluded conditions (e.g., sicca)
if (CONDITIONS_EXCLUDE != "" && "condition" %in% colnames(seu@meta.data)) {
  conds_to_exclude <- strsplit(CONDITIONS_EXCLUDE, ",")[[1]]
  cells_before <- ncol(seu)
  keep <- !seu$condition %in% conds_to_exclude
  if (sum(!keep) > 0) {
    seu <- subset(seu, cells = colnames(seu)[keep])
    log_msg("  Excluded conditions: ", paste(conds_to_exclude, collapse = ", "))
    log_msg("  Cells removed: ", cells_before - ncol(seu))
    log_msg("  Cells remaining: ", ncol(seu))
  }
}

# Note: Per-condition subsampling already done per-dataset during loading
# This keeps memory manageable during merge
log_msg("Cells per condition after merge:")
print(table(seu$condition))

# Join layers if Seurat v5
seu <- safe_join_layers(seu)

################################################################################
# 3. Standard Seurat workflow
################################################################################

log_msg("Running Seurat workflow (normalisation only; PCA/Harmony/UMAP stripped — not consumed by any downstream script; Fig 1 UMAP comes from scVI latent in results/integration/scvi_full_gpu/)...")

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)

################################################################################
# 4. Generate summary statistics
################################################################################

log_msg("Generating summary statistics...")

# Overall summary
summary_overall <- data.frame(
  metric = c("Total cells", "Total genes", "Datasets", "Diseases",
             "Cell types", "Donors"),
  value = c(
    ncol(seu),
    nrow(seu),
    n_distinct(seu$dataset_id),
    n_distinct(seu$disease),
    n_distinct(seu$cell_type),
    n_distinct(seu$donor)
  )
)

write.table(
  summary_overall,
  file.path(ATLAS_OVERVIEW_DIR, "atlas_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Per-dataset summary
summary_dataset <- seu@meta.data %>%
  group_by(dataset_id, disease) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(donor),
    n_celltypes = n_distinct(cell_type),
    median_genes = median(nFeature_RNA, na.rm = TRUE),
    median_umis = median(nCount_RNA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.table(
  summary_dataset,
  file.path(ATLAS_OVERVIEW_DIR, "atlas_summary_by_dataset.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Per-disease summary
summary_disease <- seu@meta.data %>%
  group_by(disease) %>%
  summarise(
    n_cells = n(),
    n_donors = n_distinct(donor),
    n_datasets = n_distinct(dataset_id),
    n_celltypes = n_distinct(cell_type),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.table(
  summary_disease,
  file.path(ATLAS_OVERVIEW_DIR, "atlas_summary_by_disease.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Per-celltype summary
summary_celltype <- seu@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    n_cells = n(),
    pct_of_total = round((n() / ncol(seu)) * 100, 2),
    n_donors = n_distinct(donor),
    n_datasets = n_distinct(dataset_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write.table(
  summary_celltype,
  file.path(ATLAS_OVERVIEW_DIR, "atlas_summary_by_celltype.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 7. Save integrated object
################################################################################

log_msg("Saving integrated atlas...")
saveRDS(seu, file.path(OBJ_DIR, "atlas_integrated.rds"))

################################################################################
# 8. Optional: extract counts + metadata for scVI training (visualisation only)
#    Triggered by env var SCVI_EXTRACT=1. Writes a sparse Matrix Market file +
#    per-cell metadata + per-gene table to results/integration/
#    scvi_visualisation/extracted/. Restricted to the 25 common cell types
#    (matches Fig 1 panels b-e) and the top 2000 HVGs. scripts/04b_train_scvi.py
#    consumes these files; the resulting latent UMAP is used for Fig 1 only,
#    not for any downstream DE/GSEA/composition/classifier step.
################################################################################

if (as.integer(Sys.getenv("SCVI_EXTRACT", "0")) == 1L) {
  suppressPackageStartupMessages(library(Matrix))
  SCVI_OUT <- file.path(ATLAS_ROOT, "results", "integration",
                        "scvi_visualisation", "extracted")
  dir.create(SCVI_OUT, recursive = TRUE, showWarnings = FALSE)
  log_msg("scVI extract requested (SCVI_EXTRACT=1)")

  scvi_seu <- subset(seu, cells = colnames(seu)[
    is_common_celltype(as.character(seu$cell_type))])
  log_msg("  Restricted to 25 common cell types: ", ncol(scvi_seu), " cells")

  hvgs <- VariableFeatures(scvi_seu)
  if (length(hvgs) == 0) hvgs <- VariableFeatures(
    FindVariableFeatures(scvi_seu, nfeatures = 2000, verbose = FALSE))
  stopifnot(length(hvgs) >= 1500)
  hvgs <- head(hvgs, 2000)
  log_msg("  HVGs: ", length(hvgs))

  counts <- LayerData(scvi_seu, assay = "RNA", layer = "counts")[hvgs, , drop = FALSE]
  Matrix::writeMM(counts, file.path(SCVI_OUT, "counts.mtx"))
  writeLines(colnames(counts), file.path(SCVI_OUT, "barcodes.tsv"))
  writeLines(rownames(counts), file.path(SCVI_OUT, "features.tsv"))
  write.table(
    data.frame(
      barcode    = colnames(counts),
      dataset_id = as.character(scvi_seu$dataset_id),
      condition  = as.character(scvi_seu$condition),
      cell_type  = as.character(scvi_seu$cell_type),
      donor_id   = as.character(scvi_seu$donor_id),
      stringsAsFactors = FALSE),
    file.path(SCVI_OUT, "metadata.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  Wrote scVI extract to: ", SCVI_OUT)
}

log_msg("=== Atlas Integration Complete ===")
log_msg("Output files:")
log_msg("  - atlas_integrated.rds")
log_msg("  - atlas_summary.tsv")
log_msg("  - atlas_summary_by_dataset.tsv")
log_msg("  - atlas_summary_by_disease.tsv")
log_msg("  - atlas_summary_by_celltype.tsv")
