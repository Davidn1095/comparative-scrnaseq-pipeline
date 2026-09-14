################################################################################
# 00_utils.R - Shared Utility Functions for the scTSig Atlas Pipeline
#
# This file provides reusable helper functions used across multiple scripts.
# Each R script sources this file after 00_config.R.
################################################################################

# --- Logging ---

log_msg <- function(...) {
  msg <- paste0(...)
  cat(format(Sys.time(), "[%H:%M:%S] "), msg, "\n", sep = "")
}

# --- String utilities ---

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

# --- Seurat v5 layer handling ---

safe_join_layers <- function(seu, assay = "RNA") {

  if (!assay %in% Seurat::Assays(seu)) return(seu)
  if (inherits(seu[[assay]], "Assay5")) {
    lyr <- SeuratObject::Layers(seu[[assay]])
    if (length(lyr) > 1) {
      log_msg("Joining ", length(lyr), " layers in assay ", assay)
      seu <- SeuratObject::JoinLayers(seu, assay = assay)
    }
  }
  seu
}

# --- QC Filtering ---

iqr_bounds <- function(v, mult = IQR_MULT) {
  q <- quantile(v, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 7)
  iqr <- q[2] - q[1]
  c(q1 = q[1], q3 = q[2], iqr = iqr, lo = q[1] - mult * iqr, hi = q[2] + mult * iqr)
}

apply_qc_filters <- function(seu,
                              min_features = MIN_FEATURES,
                              max_mt_pct = MAX_MT_PCT,
                              max_ribo_pct = MAX_RIBO_PCT,
                              iqr_mult = IQR_MULT) {
  md <- seu@meta.data

  nf_b <- iqr_bounds(md$nFeature_RNA, iqr_mult)
  nc_b <- iqr_bounds(md$nCount_RNA, iqr_mult)
  mt_b <- iqr_bounds(md$percent.mt, iqr_mult)

  nf_lo <- max(min_features, nf_b["lo"])
  nf_hi <- nf_b["hi"]
  nc_hi <- nc_b["hi"]
  mt_hi <- min(max_mt_pct, mt_b["hi"])

  use_ribo <- "percent.ribo" %in% colnames(md) && any(is.finite(md$percent.ribo))
  if (use_ribo) {
    rb_b <- iqr_bounds(md$percent.ribo, iqr_mult)
    rb_hi <- min(max_ribo_pct, rb_b["hi"])
  }

  n_before <- ncol(seu)
  seu_pre <- seu

  seu <- tryCatch({
    if (use_ribo) {
      subset(seu,
        subset = nFeature_RNA >= nf_lo & nFeature_RNA <= nf_hi &
          nCount_RNA <= nc_hi & percent.mt <= mt_hi & percent.ribo <= rb_hi)
    } else {
      subset(seu,
        subset = nFeature_RNA >= nf_lo & nFeature_RNA <= nf_hi &
          nCount_RNA <= nc_hi & percent.mt <= mt_hi)
    }
  }, error = function(e) {
    # Returning seu_pre here handed back the UNFILTERED object and carried on, so a
    # broken subset expression produced a dataset in which every cell passed QC and
    # the run reported success. The log line below would then report
    # "N -> N cells", which reads as "nothing was filtered" rather than "filtering
    # failed". A QC step that cannot filter must stop, not pass everything.
    stop("apply_qc_filters: subsetting failed, so no cell was filtered. Refusing to ",
         "return the unfiltered object, because that is indistinguishable downstream ",
         "from data that passed QC. Underlying error: ", conditionMessage(e))
  })

  if (ncol(seu) == 0) {
    stop("apply_qc_filters: the thresholds removed every cell (nF ", round(nf_lo), "-",
         round(nf_hi), ", MT<", round(mt_hi, 1), "%). Refusing to fall back to the ",
         "unfiltered object. Check the QC metrics on this dataset.")
  }

  log_msg("QC filter: ", n_before, " -> ", ncol(seu), " cells",
          " (nF: ", round(nf_lo), "-", round(nf_hi),
          ", MT<", round(mt_hi, 1), "%",
          if (use_ribo) paste0(", Ribo<", round(rb_hi, 1), "%") else "",
          ")")
  seu
}

# --- Atlas loading ---

load_atlas <- function(atlas_path = NULL) {
  if (is.null(atlas_path)) {
    atlas_path <- ATLAS_PATH
  }
  if (!file.exists(atlas_path)) {
    stop("Integrated atlas not found at ", atlas_path,
         ". Run 04a_integrate.R first.")
  }
  seu <- readRDS(atlas_path)
  log_msg("Loaded atlas: ", ncol(seu), " cells, ",
          length(unique(seu$disease)), " diseases")
  seu
}

# --- Condition mapping for plots ---

map_plot_condition <- function(seu) {
  seu@meta.data$plot_condition <- factor(
    unname(CONDITION_MAP[as.character(seu$condition)]),
    levels = CONDITION_ORDER
  )
  seu
}

# --- Standard metadata defaults ---

ensure_metadata <- function(seu, accession = ACCESSION) {
  if (!"disease_state" %in% colnames(seu@meta.data)) seu$disease_state <- seu$condition
  if (!"donor_id" %in% colnames(seu@meta.data)) seu$donor_id <- seu$sample_id
  if (!"batch_id" %in% colnames(seu@meta.data)) seu$batch_id <- seu$sample_id
  if (!"dataset_id" %in% colnames(seu@meta.data)) seu$dataset_id <- accession
  seu
}

# --- Percent feature set ---

add_qc_metrics <- function(seu) {
  # GUARD: a QC fence whose pattern matches NOTHING silently becomes a no-op. That is
  # what happened to CXG_436154da (perez2022): the h5ad -> 10x conversion wrote Ensembl
  # IDs into the gene-symbol column, so "^MT-" and "^RPL|^RPS" matched zero features,
  # percent.mt and percent.ribo were identically zero for all 1.26M cells, the IQR
  # fences collapsed to "<= 0" and passed everything. The log recorded it as
  # "MT<0%, Ribo<0%" and nothing else complained. Recomputed on the finished object the
  # true values were 3.5% mitochondrial and 38.2% ribosomal, so the data was always
  # there and only the identifiers were wrong. Stop rather than pass everything.
  n_mt   <- length(grep("^MT-", rownames(seu)))
  n_ribo <- length(grep("^RPL|^RPS", rownames(seu)))
  if (n_mt == 0 || n_ribo == 0) {
    stop(sprintf(paste0(
      "add_qc_metrics: QC pattern matched no features (MT- %d, RPL/RPS %d) over %d genes ",
      "starting %s. The matrix is almost certainly keyed by something other than gene ",
      "symbols - check that the feature file carries symbols, not Ensembl IDs. Refusing ",
      "to continue, because a zero-match fence passes every cell."),
      n_mt, n_ribo, nrow(seu), paste(utils::head(rownames(seu), 3), collapse = ", ")))
  }
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RPL|^RPS")
  seu
}

# --- Doublet Detection using scDblFinder ---

#' Detect and remove doublets using scDblFinder
#'
#' @param seu Seurat object (should be normalized)
#' @param sample_col Column name for sample/batch grouping (doublets detected per sample)
#' @param remove_doublets If TRUE, remove doublets; if FALSE, just annotate
#' @return Seurat object with doublet_score and doublet_class columns
detect_doublets <- function(seu, sample_col = "sample_id", remove_doublets = TRUE) {
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("scDblFinder not installed - doublet detection is required. ",
         "Install with: BiocManager::install('scDblFinder')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment not installed - required for doublet detection. ",
         "Install with: BiocManager::install('SingleCellExperiment')")
  }

  n_before <- ncol(seu)
  log_msg("Running scDblFinder doublet detection on ", n_before, " cells...")

  # Join layers if Seurat v5 assay has multiple layers (multi-sample merge)
  # Required because as.SingleCellExperiment fails with multi-layer v5 assays
  if (inherits(seu[["RNA"]], "Assay5") &&
      length(SeuratObject::Layers(seu[["RNA"]])) > 1) {
    log_msg("  Joining layers in v5 assay for doublet detection...")
    seu <- SeuratObject::JoinLayers(seu)
  }

  # Convert to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(seu)

  # Run scDblFinder per sample/batch. The grouping column MATTERS: scDblFinder's
  # expected doublet rate scales with the number of cells it is given, so running one
  # pass over a whole multi-donor dataset infers a rate appropriate to a single huge
  # sample. CXG_436154da (perez2022) hit exactly this - the CELLxGENE import leaves
  # sample_id = "sample_1" for every cell when the source h5ad has no batch field, so
  # the guard below failed, scDblFinder ran once over 1.21M cells and removed 39.8% of
  # them, against 7.2-8.5% for every dataset with real per-sample IDs. Fall back to
  # donor_id, gsm_id or batch_id before giving up.
  #
  # That fallback is a safety net for cohorts whose per-sample ID failed to import, not
  # a way to make a pooled multiplexed design tractable: grouping by a biological donor
  # hides cross-donor doublets from the simulation. CXG_436154da is now exempt from
  # doublet detection entirely (DOUBLETS_PREFILTERED in 00_config.R) and no longer
  # reaches this function.
  grp <- NULL
  for (cand in c(sample_col, "donor_id", "gsm_id", "batch_id")) {
    if (cand %in% colnames(seu@meta.data) &&
        length(unique(seu@meta.data[[cand]])) > 1) { grp <- cand; break }
  }
  if (!is.null(grp)) {
    if (!identical(grp, sample_col))
      log_msg("  ", sample_col, " has a single value; grouping doublet detection by ", grp,
              " instead (", length(unique(seu@meta.data[[grp]])), " groups)")
    sce <- scDblFinder::scDblFinder(sce, samples = grp)
  } else {
    log_msg("  WARNING: no column with >1 unique value found; running scDblFinder ",
            "ungrouped over all ", ncol(seu), " cells. Expect an inflated doublet rate.")
    sce <- scDblFinder::scDblFinder(sce)
  }

  # Transfer annotations back to Seurat
  seu$doublet_score <- sce$scDblFinder.score
  seu$doublet_class <- as.character(sce$scDblFinder.class)

  n_doublets <- sum(seu$doublet_class == "doublet", na.rm = TRUE)
  pct_doublets <- round(100 * n_doublets / n_before, 1)
  log_msg("  Detected ", n_doublets, " doublets (", pct_doublets, "%)")

  if (remove_doublets && n_doublets > 0) {
    seu <- subset(seu, subset = doublet_class == "singlet")
    log_msg("  Removed doublets: ", n_before, " -> ", ncol(seu), " cells")
  }

  seu
}

# =============================================================================
# Sex Inference from Gene Expression
# =============================================================================

#' Infer sex from gene expression using XIST and Y-chromosome genes
#'
#' Uses XIST expression (X-linked, female-specific) and Y-chromosome genes
#' to predict biological sex from single-cell or bulk RNA-seq data.
#'
#' @param seu Seurat object with normalized expression data
#' @param xist_thresh XIST expression threshold (log-normalized) for female call
#' @param y_genes Character vector of Y-chromosome genes to use
#' @param min_y_expr Minimum mean Y-gene expression for male call
#' @param assay Assay to use (default: "RNA")
#' @return Seurat object with inferred_sex column added to metadata
#'
#' @details
#' Classification logic:
#' - Female: High XIST (>thresh) AND low Y-gene expression
#' - Male: Low XIST AND high Y-gene expression (>min_y_expr)
#' - Ambiguous: Both high or both low (flagged for manual review)
#'
#' Common Y-chromosome genes used: RPS4Y1, EIF1AY, DDX3Y, KDM5D, UTY, USP9Y
#'
#' @examples
#' seu <- infer_sex_from_expression(seu)
#' table(seu$inferred_sex)
infer_sex_from_expression <- function(seu,
                                       xist_thresh = 0.5,
                                       y_genes = c("RPS4Y1", "EIF1AY", "DDX3Y",
                                                   "KDM5D", "UTY", "USP9Y", "ZFY"),
                                       min_y_expr = 0.3,
                                       assay = "RNA") {

  log_msg("Inferring sex from gene expression...")

  # Handle Seurat v5 multiple layers
  seu <- safe_join_layers(seu, assay = assay)

  # Get expression matrix
  expr <- Seurat::GetAssayData(seu, layer = "data", assay = assay)
  genes_available <- rownames(expr)

  # Check for XIST
  has_xist <- "XIST" %in% genes_available
  if (!has_xist) {
    log_msg("  WARNING: XIST not found in expression matrix")
  }


  # Check for Y-chromosome genes
  y_genes_found <- intersect(y_genes, genes_available)
  if (length(y_genes_found) == 0) {
    log_msg("  WARNING: No Y-chromosome genes found")
  } else {
    log_msg("  Found ", length(y_genes_found), "/", length(y_genes), " Y-chromosome genes: ",
            paste(y_genes_found, collapse = ", "))
  }

  # Calculate per-cell expression
  n_cells <- ncol(expr)

  # XIST expression
  if (has_xist) {
    xist_expr <- as.numeric(expr["XIST", ])
  } else {
    xist_expr <- rep(0, n_cells)
  }

  # Mean Y-gene expression
  if (length(y_genes_found) > 0) {
    y_expr_mat <- expr[y_genes_found, , drop = FALSE]
    y_mean_expr <- Matrix::colMeans(y_expr_mat)
  } else {
    y_mean_expr <- rep(0, n_cells)
  }

  # Classify
  inferred_sex <- rep("ambiguous", n_cells)

  # Female: high XIST, low Y
  female_mask <- (xist_expr > xist_thresh) & (y_mean_expr < min_y_expr)
  inferred_sex[female_mask] <- "female"

  # Male: low XIST, high Y
  male_mask <- (xist_expr <= xist_thresh) & (y_mean_expr >= min_y_expr)
  inferred_sex[male_mask] <- "male"

  # Add to metadata
  seu$inferred_sex <- inferred_sex
  seu$xist_expr <- xist_expr
  seu$y_gene_expr <- y_mean_expr

  # Summary
  sex_table <- table(inferred_sex)
  log_msg("  Sex inference results:")
  for (s in names(sex_table)) {
    pct <- round(100 * sex_table[s] / n_cells, 1)
    log_msg("    ", s, ": ", sex_table[s], " cells (", pct, "%)")
  }

  seu
}

#' Infer sex at donor/sample level from single-cell data
#'
#' Aggregates cell-level sex inference to donor level using majority voting.
#'
#' @param seu Seurat object with inferred_sex column (from infer_sex_from_expression)
#' @param donor_col Column name for donor/sample grouping
#' @param min_agreement Minimum proportion of cells that must agree (default: 0.8)
#' @return Data frame with donor_id, inferred_sex, n_cells, agreement columns
#'
#' @examples
#' seu <- infer_sex_from_expression(seu)
#' donor_sex <- infer_donor_sex(seu, donor_col = "donor_id")
infer_donor_sex <- function(seu, donor_col = "donor_id", min_agreement = 0.8) {

  if (!"inferred_sex" %in% colnames(seu@meta.data)) {
    stop("Run infer_sex_from_expression() first")
  }

  meta <- seu@meta.data

  # Aggregate by donor
  donors <- unique(meta[[donor_col]])
  rows <- list()

  for (d in donors) {
    cells <- meta[meta[[donor_col]] == d, ]
    n <- nrow(cells)

    n_female <- sum(cells$inferred_sex == "female")
    n_male <- sum(cells$inferred_sex == "male")
    n_ambiguous <- sum(cells$inferred_sex == "ambiguous")

    # Majority vote (excluding ambiguous)
    if (n_female > n_male) {
      sex <- "female"
      agreement <- n_female / (n_female + n_male + 0.001)
    } else if (n_male > n_female) {
      sex <- "male"
      agreement <- n_male / (n_female + n_male + 0.001)
    } else {
      sex <- "ambiguous"
      agreement <- 0
    }

    # Flag low agreement
    if (agreement < min_agreement && sex != "ambiguous") {
      sex <- paste0(sex, "_low_conf")
    }

    rows[[length(rows) + 1]] <- data.frame(
      donor_id = d,
      inferred_sex = sex,
      n_cells = n,
      n_female = n_female,
      n_male = n_male,
      n_ambiguous = n_ambiguous,
      agreement = round(agreement, 3),
      stringsAsFactors = FALSE
    )
  }

  results <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(
    donor_id = character(), inferred_sex = character(),
    n_cells = integer(), n_female = integer(), n_male = integer(),
    n_ambiguous = integer(), agreement = numeric(),
    stringsAsFactors = FALSE
  )

  log_msg("Donor-level sex inference:")
  log_msg("  ", sum(grepl("^female", results$inferred_sex)), " female donors")
  log_msg("  ", sum(grepl("^male", results$inferred_sex)), " male donors")
  log_msg("  ", sum(results$inferred_sex == "ambiguous"), " ambiguous donors")

  results
}
