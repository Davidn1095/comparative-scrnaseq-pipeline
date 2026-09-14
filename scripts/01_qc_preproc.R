#!/usr/bin/env Rscript
################################################################################
# Script 01: QC Filtering and Preprocessing
#
# Purpose: Load raw scRNA-seq data (10X MTX, H5AD, H5, or pre-built RDS),
#          perform quality control, normalization, and feature selection.
#
# QC Pipeline:
#   1. Basic threshold filtering (nFeature, nCount, %MT, %ribo)
#   2. Doublet detection (scDblFinder)
#   3. Normalization and variable feature selection
#
# Input:
#   - Raw data in multiple formats (auto-detected)
#   - Metadata from meta/samples.tsv
#
# Output:
#   - {ACCESSION}_qc_preproc.rds (Seurat object)
#
# Usage:
#   export DISEASE=SLE ACCESSION=GSE162577
#   export ACC_DIR="$ATLAS_ROOT/data/$DISEASE/$ACCESSION"
#   Rscript scripts/01_qc_preproc.R
################################################################################

# Source shared configuration and utilities
force_qc <- tolower(Sys.getenv("FORCE_QC", "0")) %in% c("1", "true", "yes")
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

# Basilisk env is baked into the container (BASILISK_USE_SYSTEM_DIR=1).
# No need to manually set RETICULATE_PYTHON.

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  if (requireNamespace("hdf5r", quietly = TRUE)) library(hdf5r)
})

# --- Derived paths ---
BASE <- ACC_DIR
RAW_DIR <- file.path(BASE, "raw")
EX_DIR <- file.path(RAW_DIR, "extracted")
EX_FALLBACK <- file.path(BASE, "extracted")
META_DIR <- file.path(BASE, "meta")
LOG_DIR <- file.path(ATLAS_ROOT, "logs")
# OUTDIR: per-dataset QC/preproc Seurat output (BASE/derived/objects)
OUTDIR <- file.path(BASE, "derived", "objects")
CELL_BATCH <- file.path(META_DIR, "cell_batch.tsv")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Locate samples metadata ---
samples_path <- file.path(META_DIR, "samples.tsv")
if (!file.exists(samples_path)) {
  candidates <- list.files(META_DIR, pattern = "samples.*\\.tsv$", full.names = TRUE)
  samples_path <- if (length(candidates) > 0) candidates[1] else ""
}
if (samples_path == "") stop("Missing samples.tsv under ", META_DIR)

rds_out <- file.path(OUTDIR, paste0(ACCESSION, "_qc_preproc.rds"))

################################################################################
# Helper: finalize object (QC + normalize + save)
################################################################################

finalize_object <- function(x) {
  x <- add_qc_metrics(x)
  x <- apply_qc_filters(x)

  # Normalize before doublet detection (scDblFinder works on normalized data).
  # Returning x unchanged on error left the object with no data layer, which
  # scDblFinder and every downstream stage then consumed as though it had been
  # normalised. Stop instead: an unnormalised object is not a degraded result, it is
  # a different thing wearing the same filename.
  x <- tryCatch({
    NormalizeData(x, verbose = FALSE)
  }, error = function(e) {
    stop("NormalizeData failed for ", ACCESSION, ". Refusing to continue with an ",
         "unnormalised object, because doublet detection and all downstream stages ",
         "assume the data layer exists. Underlying error: ", conditionMessage(e))
  })

  # Doublet removal, per cohort design. Cohorts the original authors already
  # doublet-filtered are listed in DOUBLETS_PREFILTERED (00_config.R, with the full
  # rationale) and are skipped here: they have had doublets removed once, by the
  # method their design calls for, and a second scDblFinder pass would remove
  # singlets rather than doublets. Everything else is multiplex-free and gets
  # scDblFinder, run per sample to account for sample-specific doublet rates.
  if (ACCESSION %in% DOUBLETS_PREFILTERED) {
    log_msg("Skipping scDblFinder for ", ACCESSION, ": the deposited object is ",
            "already doublet-filtered by the original authors (demuxlet on a pooled ",
            "design, median 16 donors per 10x lane). Cross-donor doublets are ",
            "already removed; a second pass would subtract singlets. See ",
            "DOUBLETS_PREFILTERED in 00_config.R before re-adding this.")
  } else {
    x <- detect_doublets(x, sample_col = "sample_id", remove_doublets = TRUE)
  }

  # --- Infer missing sex from gene expression ---
  # Check if samples.tsv has missing sex values and infer from expression
  if (file.exists(samples_path)) {
    samples <- read.table(samples_path, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)

    # Check if sex column is missing or has NA/empty values
    needs_sex_inference <- FALSE
    if (!"sex" %in% colnames(samples)) {
      samples$sex <- NA
      needs_sex_inference <- TRUE
    } else if (any(is.na(samples$sex) | samples$sex == "")) {
      needs_sex_inference <- TRUE
    }

    if (needs_sex_inference) {
      log_msg("Inferring sex from gene expression for samples with missing annotation...")

      # Infer sex at cell level
      x <- tryCatch({
        infer_sex_from_expression(x)
      }, error = function(e) {
        log_msg("Sex inference failed: ", conditionMessage(e))
        x
      })

      # Aggregate to donor level if cell-level inference succeeded
      if ("inferred_sex" %in% colnames(x@meta.data)) {
        donor_col <- if ("donor_id" %in% colnames(x@meta.data)) "donor_id" else "sample_id"
        donor_sex <- tryCatch({
          infer_donor_sex(x, donor_col = donor_col)
        }, error = function(e) {
          log_msg("Donor-level sex aggregation failed: ", conditionMessage(e))
          NULL
        })

        if (!is.null(donor_sex)) {
          # Add sex_source column if needed
          if (!"sex_source" %in% colnames(samples)) {
            samples$sex_source <- NA
          }

          # Match by donor_id or sample_id
          match_col <- if ("donor_id" %in% colnames(samples)) "donor_id" else "sample_id"

          n_inferred <- 0
          for (i in seq_len(nrow(samples))) {
            donor <- samples[[match_col]][i]
            match_idx <- which(donor_sex$donor_id == donor)

            if (length(match_idx) > 0) {
              inferred <- donor_sex$inferred_sex[match_idx[1]]
              # Only update if currently missing and inference is confident
              if (is.na(samples$sex[i]) || samples$sex[i] == "") {
                if (!grepl("ambiguous|low_conf", inferred)) {
                  samples$sex[i] <- inferred
                  samples$sex_source[i] <- "inferred_expression"
                  n_inferred <- n_inferred + 1
                }
              }
            }
          }

          # Save updated samples.tsv
          write.table(samples, samples_path, sep = "\t", quote = FALSE, row.names = FALSE)
          log_msg("Inferred sex for ", n_inferred, " samples, updated: ", samples_path)
        }
      }
    }
  }

  # Find variable features after all QC steps. Returning x unchanged on error saved an
  # object with no variable features, which 04a then merged and the scVI extract read
  # for its HVG selection - silently producing a different feature set than intended.
  x <- tryCatch({
    FindVariableFeatures(x, nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)
  }, error = function(e) {
    stop("FindVariableFeatures failed for ", ACCESSION, ". Refusing to save an object ",
         "with no variable features, because downstream HVG selection would silently ",
         "use a different feature set. Underlying error: ", conditionMessage(e))
  })

  saveRDS(x, rds_out)
  log_msg("Saved: ", rds_out, " (", ncol(x), " cells, ", nrow(x), " genes)")
}

################################################################################
# Pre-built raw objects are refused
################################################################################
# This script once preferred an <accession>_seurat_raw*.rds, when one existed,
# over the raw data. No script in the pipeline writes such a file, so one found
# here was built elsewhere and would silently bypass QC from the raw files.
# Stop and name it instead.

raw_candidates <- c(
  list.files(OUTDIR, pattern = paste0("^", ACCESSION, "_seurat_raw.*\\.rds$"), full.names = TRUE),
  list.files(file.path(BASE, "objects"), pattern = paste0("^", ACCESSION, "_seurat_raw.*\\.rds$"), full.names = TRUE)
)
raw_candidates <- raw_candidates[file.exists(raw_candidates)]
if (length(raw_candidates) > 0) {
  stop("Found pre-built object(s) that this pipeline does not create: ",
       paste(raw_candidates, collapse = ", "), ". 01_qc_preproc.R builds each dataset from ",
       "the files under ", RAW_DIR, "; move or delete the object(s) and rerun.", call. = FALSE)
}

################################################################################
# Path 2a: Pre-converted 10X from h5ad (avoids reticulate memory issues)
################################################################################

converted_10x_dir <- file.path(RAW_DIR, "converted_10x")
if (dir.exists(converted_10x_dir) &&
    file.exists(file.path(converted_10x_dir, "matrix.mtx.gz"))) {

  # Skip if already processed, UNLESS forced.
  #
  # This skip is an operational trap and has bitten once. The wrapper runs 01 -> 02 ->
  # 03 in sequence, so when 01 exits early the run does NOT stop: it goes straight on
  # to ANNOTATE THE STALE OBJECT and rebuild pseudobulk from it, reporting success
  # throughout. That is how a rerun intended to fix the Perez ingestion silently
  # re-annotated the very object it was meant to replace. The skip is still the right
  # default for resuming an interrupted batch, but it now announces itself loudly and
  # can be overridden.
  #
  #   FORCE_QC=1   recompute and overwrite, ignoring any existing rds
  force_qc <- tolower(Sys.getenv("FORCE_QC", "0")) %in% c("1", "true", "yes")
  if (file.exists(rds_out) && !force_qc) {
    log_msg(strrep("!", 74))
    log_msg("!! SKIPPING QC: output already exists")
    log_msg("!!   ", rds_out)
    log_msg("!!   modified ", format(file.info(rds_out)$mtime, "%Y-%m-%d %H:%M:%S"))
    log_msg("!! Downstream steps in this wrapper will now run against THAT object,")
    log_msg("!! not a fresh one. If you changed QC, 00_utils.R, the container or the")
    log_msg("!! raw input, this is NOT what you want.")
    log_msg("!! Re-run with FORCE_QC=1 to recompute, or delete the file above.")
    log_msg(strrep("!", 74))
    quit(save = "no", status = 0)
  }
  if (file.exists(rds_out) && force_qc)
    log_msg("FORCE_QC=1: recomputing and overwriting ", rds_out)

  log_msg("Found pre-converted 10X data at ", converted_10x_dir)
  x <- tryCatch({
    counts <- Seurat::ReadMtx(
      mtx = file.path(converted_10x_dir, "matrix.mtx.gz"),
      features = file.path(converted_10x_dir, "features.tsv.gz"),
      cells = file.path(converted_10x_dir, "barcodes.tsv.gz")
    )
    obj <- Seurat::CreateSeuratObject(counts = counts, project = ACCESSION)
    log_msg("  Loaded: ", ncol(obj), " cells x ", nrow(obj), " genes")

    # Add metadata from CSV
    meta_csv <- file.path(converted_10x_dir, "metadata.csv")
    if (file.exists(meta_csv)) {
      obs_df <- read.csv(meta_csv, row.names = 1, check.names = FALSE)
      obs_df <- obs_df[match(colnames(obj), rownames(obs_df)), , drop = FALSE]
      obj <- AddMetaData(obj, metadata = obs_df)
      log_msg("  Added metadata: ", ncol(obs_df), " columns")
    }
    obj
  }, error = function(e) {
    log_msg("  Pre-converted 10X read failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(x)) {
    # Map conditions
    if ("mapped_condition" %in% colnames(x@meta.data) && all(!is.na(x$mapped_condition))) {
      x$condition <- x$mapped_condition
      if (!"sample_id" %in% colnames(x@meta.data)) x$sample_id <- "sample_1"
    } else {
      ann <- readr::read_tsv(samples_path, show_col_types = FALSE)
      stopifnot(all(c("sample_id", "gsm_id", "condition") %in% colnames(ann)))
      if (!"sample_id" %in% colnames(x@meta.data)) {
        if (all(c("donor_id", "batch") %in% colnames(x@meta.data))) {
          x$sample_id <- paste(x$donor_id, x$batch, sep = "__")
        } else if (all(c("donor_id", "batch_id") %in% colnames(x@meta.data))) {
          x$sample_id <- paste(x$donor_id, x$batch_id, sep = "__")
        } else {
          x$sample_id <- "sample_1"
        }
      }
      m_condition <- setNames(as.character(ann$condition), as.character(ann$sample_id))
      m_gsm <- setNames(as.character(ann$gsm_id), as.character(ann$sample_id))
      x$condition <- unname(m_condition[x$sample_id])
      x$gsm_id <- unname(m_gsm[x$sample_id])
    }

    # CXG fallback: the CELLxGENE-converted h5ad path leaves sample_id="sample_1"
    # for all cells when no batch/batch_id is present, which makes the samples.tsv
    # mapping above fail silently (sample_id key mismatch with the author-format
    # donor IDs in samples.tsv → condition all-NA). Derive condition from the
    # CELLxGENE `disease` ontology field as a last resort. This is the case for
    # CXG_436154da (perez2022) — see also Methods §4.1.
    if (all(is.na(x$condition))) {
      if ("disease" %in% colnames(x@meta.data)) {
        dx <- tolower(as.character(x$disease))
        x$condition <- ifelse(dx == "normal", "healthy",
                       ifelse(grepl("lupus", dx), "sle",
                       ifelse(grepl("sj.gren", dx), "sjs", NA_character_)))
        log_msg("  Derived condition from `disease` (CXG-style fallback). ",
                "Non-NA cells: ", sum(!is.na(x$condition)), " of ", ncol(x))
      } else if ("disease_state" %in% colnames(x@meta.data)) {
        x$condition <- as.character(x$disease_state)
        log_msg("  Derived condition from `disease_state` fallback.")
      }
    }

    x <- ensure_metadata(x)
    finalize_object(x)
    quit(save = "no", status = 0)
  }
}

################################################################################
# Path 2b: H5AD files (via reticulate/anndata)
################################################################################

h5ad_candidates <- list.files(RAW_DIR, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE)
if (length(h5ad_candidates) > 0) {
  set.seed(1)
  h5ad_in <- h5ad_candidates[1]
  x <- NULL

  # Try baked-in basilisk Python if RETICULATE_PYTHON is not explicitly set
  if (Sys.getenv("RETICULATE_PYTHON", "") == "") {
    baked_py <- file.path(find.package("zellkonverter", quiet = TRUE),
      "basilisk", "zellkonverterAnnDataEnv-0.10.2", "bin", "python3")
    if (length(baked_py) == 1 && file.exists(baked_py)) Sys.setenv(RETICULATE_PYTHON = baked_py)
  }

  if (Sys.getenv("RETICULATE_PYTHON", "") != "") {
    x <- tryCatch({
      ad <- reticulate::import("anndata")

      adata <- ad$read_h5ad(h5ad_in)
      log_msg("h5ad has ", adata$n_obs, " cells total")

      if (!is.null(adata$raw)) {
        counts_py <- adata$raw$X
        gene_names <- adata$raw$var_names$tolist()
      } else {
        counts_py <- adata$X
        gene_names <- adata$var_names$tolist()
      }
      cell_names <- adata$obs_names$tolist()
      obs_csv <- tempfile(fileext = ".csv")
      adata$obs$to_csv(obs_csv)
      obs_df <- read.csv(obs_csv, row.names = 1, check.names = FALSE)
      unlink(obs_csv)

      # Convert sparse matrix to R (transpose: cells x genes -> genes x cells)
      counts <- reticulate::py_to_r(counts_py)
      counts <- Matrix::t(counts)
      rownames(counts) <- gene_names
      colnames(counts) <- cell_names

      # Convert Ensembl IDs to gene symbols if needed
      if (mean(grepl("^ENSG[0-9]+", gene_names)) > 0.5) {
        log_msg("Detected Ensembl IDs, converting to gene symbols...")
        mapping <- AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = gene_names, keytype = "ENSEMBL",
          column = "SYMBOL", multiVals = "first"
        )
        new_names <- ifelse(is.na(mapping[gene_names]), gene_names, mapping[gene_names])
        rownames(counts) <- make.unique(new_names)
      }

      obj <- Seurat::CreateSeuratObject(counts = counts, project = ACCESSION)
      obs_df <- obs_df[match(colnames(obj), rownames(obs_df)), , drop = FALSE]
      AddMetaData(obj, metadata = obs_df)
    }, error = function(e) {
      log_msg("reticulate/anndata read failed: ", conditionMessage(e))
      NULL
    })
  }

  if (is.null(x)) stop("Failed to read h5ad. Ensure RETICULATE_PYTHON is set.")

  # Map conditions
  if ("mapped_condition" %in% colnames(x@meta.data) && all(!is.na(x$mapped_condition))) {
    x$condition <- x$mapped_condition
    if (!"sample_id" %in% colnames(x@meta.data)) x$sample_id <- "sample_1"
  } else {
    ann <- read_tsv(samples_path, show_col_types = FALSE)
    stopifnot(all(c("sample_id", "gsm_id", "condition") %in% colnames(ann)))

    if (!"sample_id" %in% colnames(x@meta.data)) {
      if (all(c("donor_id", "batch") %in% colnames(x@meta.data))) {
        x$sample_id <- paste(x$donor_id, x$batch, sep = "__")
      } else if (all(c("donor_id", "batch_id") %in% colnames(x@meta.data))) {
        x$sample_id <- paste(x$donor_id, x$batch_id, sep = "__")
      } else {
        x$sample_id <- "sample_1"
      }
    }

    m_condition <- setNames(as.character(ann$condition), as.character(ann$sample_id))
    m_gsm <- setNames(as.character(ann$gsm_id), as.character(ann$sample_id))
    x$condition <- unname(m_condition[x$sample_id])
    x$gsm_id <- unname(m_gsm[x$sample_id])
  }

  x <- ensure_metadata(x)
  finalize_object(x)
  quit(save = "no", status = 0)
}

################################################################################
# Path 3: Extracted 10X files (MTX/barcodes or H5)
################################################################################

ex_path <- if (dir.exists(EX_DIR)) EX_DIR else if (dir.exists(EX_FALLBACK)) EX_FALLBACK else ""
if (ex_path == "") stop("Missing extracted matrices under ", EX_DIR, " or ", EX_FALLBACK)

# SECOND skip-on-exists, on the extracted-10X path used by the GEO cohorts. This one
# was entirely SILENT - no log line at all - so a rerun of any GEO dataset would patch
# metadata on the existing object and quit, while the wrapper carried on to annotate
# and pseudobulk it as though QC had run. Honours the same FORCE_QC as the converted
# 10X path above.
if (file.exists(rds_out) && !force_qc) {
  log_msg(strrep("!", 74))
  log_msg("!! SKIPPING QC (extracted 10X path): output already exists")
  log_msg("!!   ", rds_out)
  log_msg("!!   modified ", format(file.info(rds_out)$mtime, "%Y-%m-%d %H:%M:%S"))
  log_msg("!! Only metadata is patched below. Downstream steps will run against THAT")
  log_msg("!! object. Re-run with FORCE_QC=1 to recompute, or delete the file above.")
  log_msg(strrep("!", 74))
  x <- readRDS(rds_out)
  if (!"condition" %in% colnames(x@meta.data) && "disease_state" %in% colnames(x@meta.data)) {
    x$condition <- x$disease_state
  }
  if (!"donor_id" %in% colnames(x@meta.data) && "subject_id" %in% colnames(x@meta.data)) {
    x$donor_id <- x$subject_id
  }
  saveRDS(x, rds_out)
  quit(save = "no", status = 0)
}

ann <- read.delim(samples_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
need <- c("dataset_id", "batch_id", "sample_id", "gsm_id", "donor_id", "condition", "age_group")
stopifnot(all(need %in% colnames(ann)))

bc_files <- sort(list.files(ex_path, pattern = "_barcodes\\.tsv\\.gz$", full.names = TRUE))
h5_files <- sort(list.files(ex_path, pattern = "\\.h5$", full.names = TRUE))

# --- Path 3a: Combined MTX with cell_batch.tsv ---
if (file.exists(CELL_BATCH) && length(bc_files) == 1) {
  cb <- read.delim(CELL_BATCH, stringsAsFactors = FALSE)
  if (all(c("Cell", "batch") %in% colnames(cb))) {
    prefix <- sub("_barcodes\\.tsv\\.gz$", "", basename(bc_files[1]))
    if (!prefix %in% ann$sample_id) {
      log_msg("Detected combined MTX with cell_batch.tsv, splitting by sample...")
      mtx_file <- file.path(ex_path, paste0(prefix, "_matrix.mtx.gz"))
      genes_file <- file.path(ex_path, paste0(prefix, "_features.tsv.gz"))
      if (!file.exists(genes_file)) genes_file <- file.path(ex_path, paste0(prefix, "_genes.tsv.gz"))

      counts <- ReadMtx(mtx = mtx_file, features = genes_file, cells = bc_files[1])
      obj <- CreateSeuratObject(counts = counts, project = ACCESSION)

      obj$sample_id  <- unname(setNames(cb$batch, cb$Cell)[colnames(obj)])
      valid <- !is.na(obj$sample_id)
      if (sum(!valid) > 0) {
        log_msg("Dropping ", sum(!valid), " cells not found in cell_batch.tsv")
        obj <- subset(obj, cells = colnames(obj)[valid])
      }

      obj$condition  <- unname(setNames(ann$condition, ann$sample_id)[obj$sample_id])
      obj$gsm_id     <- unname(setNames(ann$gsm_id, ann$sample_id)[obj$sample_id])
      obj$donor_id   <- unname(setNames(ann$donor_id, ann$sample_id)[obj$sample_id])
      obj$batch_id   <- unname(setNames(ann$batch_id, ann$sample_id)[obj$sample_id])
      obj$dataset_id <- ACCESSION

      finalize_object(obj)
      quit(save = "no", status = 0)
    }
  }
}

# --- Helper: read 10X H5 file ---
read_10x_h5 <- function(path) {
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    f <- hdf5r::H5File$new(path, mode = "r")
    on.exit(f$close_all(), add = TRUE)
    data <- f[["matrix/data"]][]
    indices <- f[["matrix/indices"]][]
    indptr <- f[["matrix/indptr"]][]
    shape <- f[["matrix/shape"]][]
    features <- make.unique(as.character(f[["matrix/features/name"]][]))
    barcodes <- f[["matrix/barcodes"]][]
    return(Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data,
                                dims = shape, dimnames = list(features, barcodes)))
  }
  # Fallback: Python h5py
  log_msg("hdf5r unavailable, using Python h5py for: ", basename(path))
  py_bin <- Sys.getenv("RETICULATE_PYTHON", "python3")
  out_dir <- tempfile("h5read_")
  dir.create(out_dir)
  writeLines(c(
    "import h5py, numpy as np, scipy.sparse as sp, os",
    sprintf("f = h5py.File('%s', 'r')", path),
    "g = f['matrix']",
    sprintf("out_dir = '%s'", out_dir),
    "mat = sp.csc_matrix((g['data'][:].astype(np.float64), g['indices'][:], g['indptr'][:]), shape=g['shape'][:])",
    "mat.sort_indices()",
    "np.savetxt(os.path.join(out_dir, 'data.txt'), mat.data, fmt='%.6g')",
    "np.savetxt(os.path.join(out_dir, 'indices.txt'), mat.indices, fmt='%d')",
    "np.savetxt(os.path.join(out_dir, 'indptr.txt'), mat.indptr, fmt='%d')",
    "np.savetxt(os.path.join(out_dir, 'shape.txt'), np.array(g['shape'][:]), fmt='%d')",
    "features = [x.decode() if isinstance(x, bytes) else x for x in g['features/name'][:]]",
    "barcodes = [x.decode() if isinstance(x, bytes) else x for x in g['barcodes'][:]]",
    "np.savetxt(os.path.join(out_dir, 'features.txt'), features, fmt='%s')",
    "np.savetxt(os.path.join(out_dir, 'barcodes.txt'), barcodes, fmt='%s')",
    "f.close()"
  ), file.path(out_dir, "read_h5.py"))
  system2(py_bin, file.path(out_dir, "read_h5.py"), stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(out_dir, "data.txt"))) stop("Python h5py read failed for: ", path)
  mat_r <- new("dgCMatrix",
    i = scan(file.path(out_dir, "indices.txt"), what = integer(), quiet = TRUE),
    p = scan(file.path(out_dir, "indptr.txt"), what = integer(), quiet = TRUE),
    x = scan(file.path(out_dir, "data.txt"), what = numeric(), quiet = TRUE),
    Dim = scan(file.path(out_dir, "shape.txt"), what = integer(), quiet = TRUE)
  )
  rownames(mat_r) <- make.unique(readLines(file.path(out_dir, "features.txt")))
  colnames(mat_r) <- readLines(file.path(out_dir, "barcodes.txt"))
  unlink(out_dir, recursive = TRUE)
  mat_r
}

# --- Path 3b: Per-sample H5 files ---
if (length(bc_files) == 0 && length(h5_files) > 0) {
  sample_id <- sub("\\.h5$", "", basename(h5_files))
  stopifnot(all(sample_id %in% ann$sample_id))

  map_dataset <- setNames(ann$dataset_id, ann$sample_id)
  map_batch   <- setNames(ann$batch_id, ann$sample_id)
  map_gsm     <- setNames(ann$gsm_id, ann$sample_id)
  map_donor   <- setNames(ann$donor_id, ann$sample_id)
  map_cond    <- setNames(ann$condition, ann$sample_id)
  map_age     <- setNames(ann$age_group, ann$sample_id)

  objs <- vector("list", length(h5_files))
  for (i in seq_along(sample_id)) {
    p <- sample_id[i]
    counts <- read_10x_h5(h5_files[i])
    if (is.list(counts)) counts <- counts[[1]]
    obj <- CreateSeuratObject(counts = counts, project = ACCESSION)
    obj <- RenameCells(obj, add.cell.id = p)
    obj$dataset_id <- unname(map_dataset[p])
    obj$batch_id   <- unname(map_batch[p])
    obj$sample_id  <- p
    obj$gsm_id     <- unname(map_gsm[p])
    obj$donor_id   <- unname(map_donor[p])
    obj$condition  <- unname(map_cond[p])
    obj$age_group  <- unname(map_age[p])
    objs[[i]] <- obj
  }
  x <- merge(objs[[1]], y = objs[-1], project = ACCESSION)

} else {
  # --- Path 3c: Per-sample MTX files ---
  stopifnot(length(bc_files) > 0)
  sample_id <- sub("_barcodes\\.tsv\\.gz$", "", basename(bc_files))
  stopifnot(all(sample_id %in% ann$sample_id))

  map_dataset <- setNames(ann$dataset_id, ann$sample_id)
  map_batch   <- setNames(ann$batch_id, ann$sample_id)
  map_gsm     <- setNames(ann$gsm_id, ann$sample_id)
  map_donor   <- setNames(ann$donor_id, ann$sample_id)
  map_cond    <- setNames(ann$condition, ann$sample_id)
  map_age     <- setNames(ann$age_group, ann$sample_id)

  pick_features <- function(prefix) {
    f1 <- file.path(ex_path, paste0(prefix, "_features.tsv.gz"))
    f2 <- file.path(ex_path, paste0(prefix, "_genes.tsv.gz"))
    if (file.exists(f1)) return(f1)
    if (file.exists(f2)) return(f2)
    stop("Missing features/genes file for ", prefix)
  }

  read_counts <- function(mtx, genes, cells) {
    gene_len <- length(readLines(gzfile(genes), n = 1e6))
    cell_len <- length(readLines(gzfile(cells), n = 1e6))
    hdr <- readLines(gzfile(mtx), n = 10)
    dims_line <- hdr[!grepl("^%", hdr)][1]
    dims <- as.integer(strsplit(trimws(dims_line), "\\s+")[[1]][1:2])
    if (dims[1] == gene_len && dims[2] == cell_len) {
      return(ReadMtx(mtx = mtx, features = genes, cells = cells))
    }
    if (dims[1] == cell_len && dims[2] == gene_len) {
      m <- Matrix::t(Matrix::readMM(mtx))
      rownames(m) <- readLines(gzfile(genes), n = gene_len)
      colnames(m) <- readLines(gzfile(cells), n = cell_len)
      return(m)
    }
    stop("Matrix dimensions do not match genes/barcodes")
  }

  objs <- vector("list", length(bc_files))
  for (i in seq_along(sample_id)) {
    p <- sample_id[i]
    mtx   <- file.path(ex_path, paste0(p, "_matrix.mtx.gz"))
    genes <- pick_features(p)
    cells <- file.path(ex_path, paste0(p, "_barcodes.tsv.gz"))
    stopifnot(file.exists(mtx), file.exists(genes), file.exists(cells))

    counts <- read_counts(mtx = mtx, genes = genes, cells = cells)
    obj <- CreateSeuratObject(counts = counts, project = ACCESSION)
    obj <- RenameCells(obj, add.cell.id = p)
    obj$dataset_id <- unname(map_dataset[p])
    obj$batch_id   <- unname(map_batch[p])
    obj$sample_id  <- p
    obj$gsm_id     <- unname(map_gsm[p])
    obj$donor_id   <- unname(map_donor[p])
    obj$condition  <- unname(map_cond[p])
    obj$age_group  <- unname(map_age[p])
    objs[[i]] <- obj
  }
  x <- merge(objs[[1]], y = objs[-1], project = ACCESSION)
}

# --- Final processing for MTX/H5 paths ---
finalize_object(x)
