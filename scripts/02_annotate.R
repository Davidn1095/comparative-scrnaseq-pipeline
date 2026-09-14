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
  library(SeuratObject)
  library(SingleR)
  library(SingleCellExperiment)
})

ACC <- trimws(ACCESSION)
BASE <- ACC_DIR
RDS_IN <- get_env("RDS_IN", file.path(BASE, "derived", "objects", paste0(ACC, "_qc_preproc.rds")))
RDS_OUT <- get_env("RDS_OUT", file.path(BASE, "derived", "objects", paste0(ACC, "_qc_preproc_singler_monaco.rds")))

# --- SingleR reference ---------------------------------------------------------
# The Monaco reference is a cache, not an input. If the .rds is missing it is
# fetched with celldex::MonacoImmuneData() and saved to the same path, so a fresh
# clone builds it on first run. The celldex and ExperimentHub versions used are
# recorded in the object's metadata and in a sidecar file beside it; a cache that
# carries neither was built outside this script and its provenance is unknown.
# REF_RDS overrides the path, but must point at an existing file: an explicit
# path that is missing is an error, not something to fetch into.
monaco_sidecar <- function(path) sub("\\.rds$", ".provenance.tsv", path)

load_monaco_reference <- function(path) {
  sidecar <- monaco_sidecar(path)
  if (file.exists(path)) {
    ref <- readRDS(path)
    prov <- S4Vectors::metadata(ref)$provenance
    if (is.null(prov) && file.exists(sidecar)) {
      kv <- utils::read.delim(sidecar, header = FALSE, col.names = c("key", "value"),
                              colClasses = "character", quote = "")
      prov <- as.list(stats::setNames(kv$value, kv$key))
    }
    if (is.null(prov)) {
      log_msg("Reference provenance UNKNOWN: ", path, " carries no build record and has no ",
              basename(sidecar), "; it was not built by this script, so the celldex version ",
              "that produced it cannot be established")
    } else {
      log_msg("Reference built ", prov$fetched, " by ", prov$call, " with celldex ", prov$celldex,
              ", ExperimentHub ", prov$ExperimentHub, " (hub snapshot ", prov$hub_snapshot,
              ", Bioconductor ", prov$Bioconductor, ")")
    }
    return(ref)
  }

  if (!requireNamespace("celldex", quietly = TRUE)) {
    stop("SingleR reference not found at ", path, " and the celldex package is not installed, ",
         "so it cannot be fetched. Install celldex (Bioconductor) or place the Monaco reference ",
         ".rds at that path.", call. = FALSE)
  }
  log_msg("Reference not found at ", path, "; fetching celldex::MonacoImmuneData() from ExperimentHub")
  ref <- tryCatch(
    celldex::MonacoImmuneData(ensembl = FALSE, cell.ont = "all"),
    error = function(e) stop(
      "celldex::MonacoImmuneData() failed, so the SingleR reference could not be built at ", path,
      ". The fetch downloads from ExperimentHub and needs internet access; run the script once ",
      "on a machine that has it to build the cache. Underlying error: ",
      conditionMessage(e), call. = FALSE))
  if (!inherits(ref, "SummarizedExperiment") ||
      !"label.fine" %in% colnames(SummarizedExperiment::colData(ref))) {
    stop("celldex::MonacoImmuneData() returned an object without a label.fine column; ",
         "nothing was saved to ", path, call. = FALSE)
  }

  prov <- list(
    fetched       = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    call          = 'celldex::MonacoImmuneData(ensembl = FALSE, cell.ont = "all")',
    celldex       = as.character(utils::packageVersion("celldex")),
    ExperimentHub = as.character(utils::packageVersion("ExperimentHub")),
    hub_snapshot  = tryCatch(as.character(AnnotationHub::snapshotDate(ExperimentHub::ExperimentHub())),
                             error = function(e) "unrecorded"),
    Bioconductor  = tryCatch(as.character(BiocManager::version()), error = function(e) "unrecorded"),
    R             = paste(R.version$major, R.version$minor, sep = ".")
  )
  S4Vectors::metadata(ref)$provenance <- prov

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  ok <- tryCatch({ saveRDS(ref, tmp); file.rename(tmp, path) }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    unlink(tmp)
    stop("Fetched the Monaco reference but could not save it to ", path,
         "; check that the directory is writable.", call. = FALSE)
  }
  writeLines(paste(names(prov), unlist(prov), sep = "\t"), sidecar)
  log_msg("Saved reference to ", path, " (celldex ", prov$celldex, ", ExperimentHub ",
          prov$ExperimentHub, "); provenance written to ", sidecar)
  ref
}

ref_override <- get_env("REF_RDS", "")
if (ref_override != "" && !file.exists(ref_override)) {
  stop("REF_RDS is set to ", ref_override, " but no file exists there. Unset REF_RDS to use, ",
       "or build, the cached reference at ", SINGLER_REF, call. = FALSE)
}
ref_rds <- if (ref_override != "") ref_override else SINGLER_REF
log_msg("Using reference: ", ref_rds)
ref <- load_monaco_reference(ref_rds)

if (!file.exists(RDS_IN)) stop("Missing: ", RDS_IN)

x <- readRDS(RDS_IN)
md <- x@meta.data
need <- c("sample_id", "condition")
miss <- setdiff(need, colnames(md))
if (length(miss) > 0) stop("Missing meta columns: ", paste(miss, collapse = ","))

# Join layers BEFORE normalization for Seurat v5
x <- safe_join_layers(x)

x <- NormalizeData(x, verbose = FALSE)

# Join again after normalization just in case
x <- safe_join_layers(x)

# Determine label column based on reference type
# Load SummarizedExperiment for colData if needed
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("SummarizedExperiment package is required")
}
library(SummarizedExperiment)

lab_col <- NULL
if (inherits(ref, c("SingleCellExperiment", "SummarizedExperiment"))) {
  # SCE/SE object - use colData
  if ("label.fine" %in% colnames(colData(ref))) lab_col <- "label.fine"
  if (is.null(lab_col) && "label.main" %in% colnames(colData(ref))) lab_col <- "label.main"
  if (is.null(lab_col)) lab_col <- colnames(colData(ref))[1]
} else if (inherits(ref, "Seurat")) {
  # Seurat object - use metadata
  if ("label.fine" %in% colnames(ref@meta.data)) lab_col <- "label.fine"
  if (is.null(lab_col) && "label.main" %in% colnames(ref@meta.data)) lab_col <- "label.main"
  if (is.null(lab_col)) lab_col <- colnames(ref@meta.data)[1]
} else {
  stop("Reference must be SingleCellExperiment, SummarizedExperiment, or Seurat object")
}

genes <- rownames(x)
ref_genes <- rownames(ref)
if (length(intersect(genes, ref_genes)) == 0) {
  g <- sub("\\..*$", "", genes)
  if (all(grepl("^ENSG", g))) {
    suppressPackageStartupMessages({
      library(org.Hs.eg.db)
      library(AnnotationDbi)
    })
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(g), columns = "SYMBOL", keytype = "ENSEMBL")
    map <- map[!is.na(map$SYMBOL), , drop = FALSE]
    sym <- map$SYMBOL[match(g, map$ENSEMBL)]
    keep <- !is.na(sym)
    x <- x[keep, ]
    rownames(x) <- make.unique(sym[keep])
  }
}

if (length(intersect(rownames(x), ref_genes)) == 0) {
  x$cell_type_singler <- "unknown"
  x$cell_type_pruned <- "unknown"
  x$singler_ref <- "monaco"
  x$accession <- ACC
  saveRDS(x, RDS_OUT)
  quit(save = "no", status = 0)
}

# Get labels from reference
if (inherits(ref, c("SingleCellExperiment", "SummarizedExperiment"))) {
  ref_labels <- colData(ref)[[lab_col]]
  # Keep as SCE/SE for SingleR - it knows how to handle these
  ref_data <- ref
} else if (inherits(ref, "Seurat")) {
  ref_labels <- ref@meta.data[[lab_col]]
  # Extract expression data from Seurat reference
  if (inherits(ref[["RNA"]], "Assay5")) {
    ref_data <- SeuratObject::LayerData(ref, assay = "RNA", layer = "data")
  } else {
    ref_data <- LayerData(ref, layer = "data", assay = "RNA")
  }
} else {
  stop("Reference must be SingleCellExperiment, SummarizedExperiment, or Seurat object")
}

# For Seurat v5, extract counts/data matrix directly instead of converting to SCE
# This avoids the layer conversion issues
message("Running SingleR annotation...")

# Get normalized counts matrix
if (inherits(x[["RNA"]], "Assay5")) {
  # Seurat v5: extract data layer
  test_mat <- SeuratObject::LayerData(x, assay = "RNA", layer = "data")
} else {
  # Seurat v3: use LayerData (works in v5 too; replaces defunct GetAssayData(slot=))
  test_mat <- LayerData(x, layer = "data", assay = "RNA")
}

# Ensure matrix format
if (!inherits(test_mat, "dgCMatrix") && !inherits(test_mat, "matrix")) {
  test_mat <- as.matrix(test_mat)
}

message(sprintf("  Test data: %d genes x %d cells", nrow(test_mat), ncol(test_mat)))

# Cell-level SingleR annotation for all datasets (uniform; no cluster-level fallback)
message("  Using cell-level annotation")
pred <- SingleR(
  test = test_mat,
  ref = ref_data,
  labels = ref_labels
)
# cell_type_singler keeps the raw SingleR calls; cell_type_pruned carries the
# pruned labels, where NA marks a low-confidence assignment. Cells with NA are
# filtered downstream (03_pseudobulk.R) and so excluded from all analyses.
x$cell_type_singler <- pred$labels
x$cell_type_pruned <- pred$pruned.labels

x$singler_ref <- "monaco"
x$accession <- ACC

message(sprintf("  Annotated %d cells", ncol(x)))
message("  Unique cell types: ", paste(unique(x$cell_type_singler), collapse = ", "))

saveRDS(x, RDS_OUT)
message("Saved annotated object to: ", RDS_OUT)
