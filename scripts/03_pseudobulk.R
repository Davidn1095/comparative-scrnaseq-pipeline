#!/usr/bin/env Rscript
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
  library(Matrix)
})

ACC <- trimws(ACCESSION)
BASE <- ACC_DIR
RDS <- get_env("RDS", file.path(BASE, "derived", "objects", paste0(ACC, "_qc_preproc_singler_monaco.rds")))
# OUTDIR: per-dataset pseudobulk counts output (BASE/derived/pseudobulk/counts)
OUTDIR <- get_env("OUTDIR", file.path(BASE, "derived", "pseudobulk", "counts"))
OUT_COUNTS <- get_env("OUT_COUNTS", file.path(OUTDIR, "pb_counts.rds"))
OUT_META <- get_env("OUT_META", file.path(OUTDIR, "pb_meta.tsv"))

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(RDS)
DefaultAssay(x) <- "RNA"

x <- safe_join_layers(x)
if (!"donor_id" %in% colnames(x@meta.data) && "subject_id" %in% colnames(x@meta.data)) {
  x$donor_id <- x$subject_id
}

need <- c("donor_id", "condition", "cell_type_pruned")
stopifnot(all(need %in% colnames(x@meta.data)))

keep <- rownames(x@meta.data)[!is.na(x@meta.data$cell_type_pruned)]
x <- subset(x, cells = keep)

x$pb_id <- paste(x$donor_id, x$condition, x$cell_type_pruned, sep = "__")

# Seurat v5: layers already joined by safe_join_layers above; use LayerData.
counts <- LayerData(x, assay = "RNA", layer = "counts")

pb <- factor(x$pb_id)
M <- Matrix::sparse.model.matrix(~ 0 + pb)
colnames(M) <- sub("^pb", "", colnames(M))
pb_counts <- counts %*% M
colnames(pb_counts) <- levels(pb)

parts <- do.call(rbind, strsplit(colnames(pb_counts), "__", fixed = TRUE))
meta <- data.frame(
  pb_id = colnames(pb_counts),
  donor_id = parts[, 1],
  condition = parts[, 2],
  cell_type = parts[, 3],
  dataset_id = ACC,
  stringsAsFactors = FALSE
)

saveRDS(pb_counts, OUT_COUNTS)
write.table(meta, OUT_META, sep = "\t", quote = FALSE, row.names = FALSE)
