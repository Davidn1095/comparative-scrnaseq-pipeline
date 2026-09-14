#!/usr/bin/env Rscript
################################################################################
# 08_donor_classifiers.R — Donor-level Hallmark-pathway classifiers
#
# Builds a 291-donor × 1250-feature matrix where each feature is the mean
# Seurat AddModuleScore for one MSigDB Hallmark pathway aggregated over a
# donor's cells of a given cell type. Trains three classifiers (ranger RF,
# XGBoost, glmnet multinomial elastic net) under stratified 5×5 CV (25 folds).
#
# Outputs (results/donor_classifiers/):
#   feature_matrix.rds          291 × 1250 numeric matrix + donor metadata
#   cv_predictions.tsv          per-fold per-donor true + predicted + prob
#   benchmark.tsv               mean ± sd metrics per classifier
#   feature_importance.tsv      top-30 features per classifier
#   confusion_matrices.pdf      3-panel CM (concatenated CV predictions)
#
# Preprocessing leak protection: median imputation + z-score normalisation
# fit inside each fold's training partition only, then applied to its test.
################################################################################

# Optional extra R library searched before the default ones, for packages
# installed outside the container image. Set ATLAS_EXTRA_R_LIB to use it.
extra_lib <- Sys.getenv("ATLAS_EXTRA_R_LIB", "")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m)) return(normalizePath(dirname(m)))
  return(".")
}
source(file.path(get_script_dir(), "00_config.R"))
source(file.path(get_script_dir(), "00_utils.R"))

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat)
  library(dplyr); library(tidyr)
  library(ranger); library(xgboost); library(glmnet)
  library(pROC); library(ggplot2); library(cowplot)
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
stopifnot(DENOISE_METHOD %in% c("none", "scvi", "limma", "limma_modulescore"))
if (DENOISE_METHOD == "scvi") {
  stop("DENOISE_METHOD=scvi builds features from scVI-denoised expression in ",
       "data/scvi_input/denoised.npy, but no script in this pipeline writes that file ",
       "(04b_train_scvi.py saves only the latent space). The option is disabled rather than ",
       "left to read a file of unknown origin. The reported classifier uses limma_modulescore.",
       call. = FALSE)
}
log_msg("Denoising method: ", DENOISE_METHOD)

OUTDIR <- switch(DENOISE_METHOD,
  scvi              = file.path(ATLAS_ROOT, "results", "donor_classifiers_scvi"),
  limma             = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma"),
  limma_modulescore = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma_modulescore"),
  file.path(ATLAS_ROOT, "results", "donor_classifiers"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# MSigDB 2025.1 Hs cache, built on first use by load_msigdb_hallmark() in 00_utils.R;
# set MSIGDB_RDS to use another location.
MSIGDB_RDS <- Sys.getenv("MSIGDB_RDS", file.path(tools::R_user_dir("msigdbr", which = "data"), "msigdb.2025.1.Hs.rds"))
FEATURE_MATRIX_PATH <- file.path(OUTDIR, "feature_matrix.rds")
SCVI_INPUT_DIR <- file.path(ATLAS_ROOT, "data", "scvi_input")

set.seed(42)

# =========================================================================
# Phase 1: Build feature matrix (291 × 1250) — cached on disk if available
# =========================================================================

# The skip below is an operational trap, and it has bitten. The feature matrix is
# derived from atlas_integrated.rds, so it goes stale the moment the atlas is rebuilt -
# but this branch reuses it silently, the classifiers train on the OLD features, and
# 09 then reports SHAP off those models. A rerun after the Perez re-ingestion looked
# like it had worked and returned performance identical to four decimal places,
# because it was the cache. Announce it, and allow an override.
#
#   FORCE_FEATURES=1   rebuild the feature matrix, ignoring any cached copy
# The banner alone was not enough. On 2026-08-12 it printed exactly as designed -
# matrix built 09:59:33, atlas 12:05:40 - and the run carried on and trained every
# classifier on the stale features anyway, because warning and proceeding is still
# proceeding. A cache must invalidate itself: reuse only while it is NEWER than the
# atlas it was derived from, and rebuild otherwise. FORCE_FEATURES=1 additionally
# forces a rebuild even when the cache is still valid.
force_features <- tolower(Sys.getenv("FORCE_FEATURES", "0")) %in% c("1", "true", "yes")
ATLAS_RDS  <- file.path(OBJ_DIR, "atlas_integrated.rds")
fm_mtime    <- if (file.exists(FEATURE_MATRIX_PATH)) file.info(FEATURE_MATRIX_PATH)$mtime else NA
atlas_mtime <- if (file.exists(ATLAS_RDS)) file.info(ATLAS_RDS)$mtime else NA
stamp <- function(x) if (is.na(x)) "missing" else format(x, "%Y-%m-%d %H:%M:%S")
cache_stale <- !is.na(fm_mtime) && !is.na(atlas_mtime) && atlas_mtime > fm_mtime

if (file.exists(FEATURE_MATRIX_PATH) && force_features)
  log_msg("FORCE_FEATURES=1: rebuilding, ignoring ", FEATURE_MATRIX_PATH)
if (cache_stale && !force_features) {
  log_msg(strrep("!", 74))
  log_msg("!! CACHED FEATURE MATRIX IS STALE - rebuilding from the current atlas")
  log_msg("!!   matrix built ", stamp(fm_mtime))
  log_msg("!!   atlas  built ", stamp(atlas_mtime), "   <-- newer")
  log_msg("!! The matrix is derived from the atlas, so it cannot be reused. Training on")
  log_msg("!! it would make every classifier metric, and all SHAP downstream of it,")
  log_msg("!! silently describe a superseded atlas.")
  log_msg(strrep("!", 74))
}
if (file.exists(FEATURE_MATRIX_PATH) && !force_features && !cache_stale) {
  log_msg("Reusing cached feature matrix (matrix ", stamp(fm_mtime),
          " is newer than atlas ", stamp(atlas_mtime), "): ", FEATURE_MATRIX_PATH)
  fm <- readRDS(FEATURE_MATRIX_PATH)
} else {
  log_msg("Building feature matrix from scratch.")

  if (DENOISE_METHOD %in% c("limma", "limma_modulescore")) {
    USE_MODULESCORE <- (DENOISE_METHOD == "limma_modulescore")
    # In-memory regeneration: load atlas, per cell type subset → apply
    # limma::removeBatchEffect → Hallmark score per cell → donor mean.
    # When USE_MODULESCORE is TRUE the per-cell Hallmark score is
    # Seurat::AddModuleScore (control-gene-normalised); otherwise it is the
    # plain colMeans over set genes. All in memory; no intermediates.
    log_msg("Loading atlas for in-memory limma regeneration (~2.5 min)...")
    if (USE_MODULESCORE) log_msg("Per-cell scoring: AddModuleScore (nbin=25, ctrl=100, seed=42)")
    else                 log_msg("Per-cell scoring: colMeans (simple mean)")
    obj <- readRDS(file.path(OBJ_DIR, "atlas_integrated.rds"))

    log_msg("Loading Hallmark gene sets...")
    hallmark_list <- load_msigdb_hallmark(MSIGDB_RDS)

    # Gene pool selection.
    #   limma:             HVGs (signal-carrying genes for the simple-mean path).
    #   limma_modulescore: union of all 50 Hallmark pathway genes. Aligns the
    #                      gene pool with the feature definition and gives
    #                      AddModuleScore enough genes for nbin=25 binning
    #                      (~4000 / 25 ≈ 160 per bin, well above ctrl=100).
    if (USE_MODULESCORE) {
      hvg <- unique(unlist(hallmark_list))
      hvg <- hvg[!grepl(EXCLUDED_GENE_REGEX, hvg)]
      hvg <- intersect(hvg, rownames(obj))
      log_msg("  Hallmark gene union: ", length(hvg),
              " genes (post-MT/RPL/RPS/HB/sex-chromosome filter, intersected with atlas)")
    } else {
      hvg <- Seurat::VariableFeatures(obj)
      hvg <- hvg[!grepl(EXCLUDED_GENE_REGEX, hvg)]
      hvg <- intersect(hvg, rownames(obj))
      log_msg("  HVG set: ", length(hvg), " genes (post-MT/RPL/RPS/HB/sex-chromosome filter)")
    }

    # Pull log-normalised expression for the chosen gene pool once, keep sparse.
    expr_hvg <- LayerData(obj, layer = "data")[hvg, , drop = FALSE]
    meta_all <- obj@meta.data
    rm(obj); gc()
    log_msg("  Pool expression: ", nrow(expr_hvg), " genes × ", ncol(expr_hvg), " cells")

    # Restrict to canonical 25 cell types
    meta_all <- meta_all[is_common_celltype(meta_all$cell_type), ]
    log_msg("  Atlas cells in 25 common cell types: ", nrow(meta_all))

    cts_in_meta <- unique(meta_all$cell_type)
    log_msg("  Processing ", length(cts_in_meta), " cell types")

    donor_ct_long <- list()
    for (ct in cts_in_meta) {
      t0 <- Sys.time()
      cl_cells <- rownames(meta_all)[meta_all$cell_type == ct]
      cl_meta <- meta_all[cl_cells, , drop = FALSE]
      # Subset + dense for limma
      cl_expr <- as.matrix(expr_hvg[, cl_cells, drop = FALSE])

      batch_fac <- factor(cl_meta$dataset_id)
      if (nlevels(batch_fac) > 1) {
        cond_fac <- factor(cl_meta$condition)
        design <- model.matrix(~cond_fac)
        cl_expr <- limma::removeBatchEffect(cl_expr, batch = batch_fac,
                                            design = design)
      }
      # cl_expr is genes × cells, limma-corrected (in-memory only)

      donor_ids <- cl_meta$donor_id
      donors <- unique(donor_ids)
      df <- data.frame(donor_id = donors,
                       cell_type = ct, stringsAsFactors = FALSE)

      if (USE_MODULESCORE) {
        # Wrap the limma-corrected matrix in a minimal Seurat object and run
        # AddModuleScore for per-cell control-gene-normalised pathway scores.
        # Counts slot is required by CreateSeuratObject; we pass a clamped
        # copy of the corrected matrix (non-negatives), then immediately
        # overwrite the data slot with the actual residuals.
        cl_expr_sparse <- as(Matrix::Matrix(cl_expr, sparse = TRUE), "dgCMatrix")
        counts_placeholder <- cl_expr_sparse
        counts_placeholder@x[counts_placeholder@x < 0] <- 0
        tmp_obj <- CreateSeuratObject(counts = counts_placeholder,
                                      assay  = "RNA")
        LayerData(tmp_obj, layer = "data", assay = "RNA") <- cl_expr_sparse
        hs_score <- lapply(hallmark_list, function(g) intersect(g, rownames(tmp_obj)))
        keep <- sapply(hs_score, length) > 0
        hs_score <- hs_score[keep]
        suppressWarnings(suppressMessages({
          tmp_obj <- AddModuleScore(tmp_obj, features = hs_score,
                                    name = "HM_", nbin = 25, ctrl = 100,
                                    seed = 42)
        }))
        score_cols <- paste0("HM_", seq_along(hs_score))
        mod_scores <- tmp_obj@meta.data[, score_cols, drop = FALSE]
        colnames(mod_scores) <- names(hs_score)
        for (h_name in names(hallmark_list)) {
          if (h_name %in% colnames(mod_scores)) {
            donor_mean <- tapply(mod_scores[[h_name]], donor_ids, mean)
            df[[h_name]] <- donor_mean[donors]
          } else {
            df[[h_name]] <- NA_real_
          }
        }
        rm(tmp_obj, mod_scores, cl_expr_sparse, counts_placeholder); gc()
      } else {
        for (h_name in names(hallmark_list)) {
          fg <- intersect(hallmark_list[[h_name]], rownames(cl_expr))
          if (length(fg) == 0) {
            df[[h_name]] <- NA_real_
          } else {
            cell_score <- colMeans(cl_expr[fg, , drop = FALSE])
            donor_mean <- tapply(cell_score, donor_ids, mean)
            df[[h_name]] <- donor_mean[donors]
          }
        }
      }
      donor_ct_long[[ct]] <- df

      rm(cl_expr); gc()
      dt_sec <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
      log_msg("    ", ct, ": ", length(donors), " donors, ", length(cl_cells),
              " cells (", dt_sec, "s)")
    }

    rm(expr_hvg); gc()

    donor_ct_means <- do.call(rbind, donor_ct_long)
    log_msg("  donor × cell_type rows: ", nrow(donor_ct_means))

    long <- donor_ct_means %>%
      pivot_longer(cols = all_of(names(hallmark_list)),
                   names_to = "pathway", values_to = "score") %>%
      mutate(feature = paste0(pathway, "__",
                              gsub("[-/ ]", "_", cell_type)))
    wide <- long %>%
      select(donor_id, feature, score) %>%
      pivot_wider(names_from = feature, values_from = score)
    donor_ids_all <- wide$donor_id
    X <- as.matrix(wide[, -1]); rownames(X) <- donor_ids_all

    log_msg("  Pre-impute: ", nrow(X), " donors × ", ncol(X), " features")
    log_msg("  NA count before imputation: ", sum(is.na(X)),
            " (", round(100 * sum(is.na(X)) / length(X), 2), "%)")

    # Donor metadata: pheno (donor_id × condition) from atlas meta_all
    donor_first <- meta_all[!duplicated(meta_all$donor_id),
                            c("donor_id", "condition")]
    rownames(donor_first) <- donor_first$donor_id
    pheno <- data.frame(
      donor_id = donor_ids_all,
      condition = donor_first[donor_ids_all, "condition"],
      stringsAsFactors = FALSE
    )
    rownames(pheno) <- donor_ids_all

    fm <- list(X_raw = X, pheno = pheno,
               feature_names = colnames(X),
               hallmark_set_sizes = sapply(hallmark_list, length))
    saveRDS(fm, FEATURE_MATRIX_PATH)
    log_msg("  Saved: ", FEATURE_MATRIX_PATH)
    skip_atlas_block <- TRUE
  } else {
    log_msg("Loading atlas (~2.5 min)...")
    obj <- readRDS(file.path(OBJ_DIR, "atlas_integrated.rds"))
    skip_atlas_block <- FALSE
  }

  if (!skip_atlas_block) {
  log_msg("Loading Hallmark gene sets...")
  hallmark_list <- load_msigdb_hallmark(MSIGDB_RDS)
  log_msg("  ", length(hallmark_list), " Hallmark sets, ",
          length(unique(unlist(hallmark_list))), " unique genes")

  # Restrict each set to genes present in the atlas
  atlas_genes <- rownames(obj)
  hallmark_list <- lapply(hallmark_list,
                          function(g) intersect(g, atlas_genes))
  set_sizes <- sapply(hallmark_list, length)
  log_msg("  Per-Hallmark set size after intersecting atlas: median=",
          median(set_sizes), " min=", min(set_sizes), " max=", max(set_sizes))
  stopifnot(all(set_sizes >= 5))

  # The block below is only for DENOISE_METHOD == "none" (Seurat AddModuleScore
  # on raw log-normalised expression).
  if (DENOISE_METHOD == "none") {
  log_msg("AddModuleScore on 50 Hallmark sets (15-25 min)...")
  obj <- AddModuleScore(obj, features = hallmark_list,
                        name = "HM_", nbin = 25, ctrl = 100, seed = 42)

  # AddModuleScore writes columns HM_1 ... HM_50 (in the order of the list).
  # Rename to the actual hallmark names.
  score_cols <- paste0("HM_", seq_along(hallmark_list))
  stopifnot(all(score_cols %in% colnames(obj@meta.data)))
  colnames(obj@meta.data)[match(score_cols, colnames(obj@meta.data))] <-
    names(hallmark_list)

  log_msg("Aggregating per donor × cell type (mean score)...")
  meta <- obj@meta.data
  rm(obj); gc()

  # Restrict to the canonical 25 cell types
  meta <- meta[is_common_celltype(meta$cell_type), ]
  log_msg("  Cells in 25 common cell types: ", nrow(meta))

  # Aggregate
  donor_ct_means <- meta %>%
    group_by(donor_id, cell_type) %>%
    summarise(across(all_of(names(hallmark_list)), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop")
  log_msg("  donor × cell_type rows: ", nrow(donor_ct_means))

  # Pivot to donor × (pathway × cell_type)
  long <- donor_ct_means %>%
    pivot_longer(cols = all_of(names(hallmark_list)),
                 names_to = "pathway", values_to = "score") %>%
    mutate(feature = paste0(pathway, "__",
                            gsub("[-/ ]", "_", cell_type)))
  # 291 donors × 25 cell types × 50 pathways = 363,750 rows (some donors
  # have no cells of a given cell type → those rows absent → NAs after pivot)

  wide <- long %>%
    select(donor_id, feature, score) %>%
    pivot_wider(names_from = feature, values_from = score)

  donor_ids <- wide$donor_id
  X <- as.matrix(wide[, -1])
  rownames(X) <- donor_ids

  log_msg("  Pre-impute: ", nrow(X), " donors × ", ncol(X), " features")
  log_msg("  NA count before imputation: ", sum(is.na(X)),
          " (", round(100 * sum(is.na(X)) / length(X), 2), "%)")

  # Donor labels — one row per unique donor from meta
  pheno <- meta %>%
    group_by(donor_id) %>%
    summarise(condition = first(condition),
              dataset_id = first(dataset_id), .groups = "drop")
  pheno <- as.data.frame(pheno)
  rownames(pheno) <- pheno$donor_id
  pheno <- pheno[donor_ids, , drop = FALSE]

  fm <- list(X_raw = X, pheno = pheno,
             feature_names = colnames(X),
             hallmark_set_sizes = set_sizes)
  saveRDS(fm, FEATURE_MATRIX_PATH)
  log_msg("  Saved: ", FEATURE_MATRIX_PATH)
  }  # close: if (DENOISE_METHOD == "none")
  }  # close: if (!skip_atlas_block)
}

X_raw <- fm$X_raw
pheno <- fm$pheno
stopifnot(identical(rownames(X_raw), rownames(pheno)))
log_msg("Feature matrix: ", nrow(X_raw), " donors × ", ncol(X_raw),
        " features. Conditions: ",
        paste(names(table(pheno$condition)), table(pheno$condition),
              sep="=", collapse=", "))
log_msg("Range of values: [",
        round(min(X_raw, na.rm=TRUE), 4), ", ",
        round(max(X_raw, na.rm=TRUE), 4), "]")
log_msg("NA count: ", sum(is.na(X_raw)),
        " (", round(100*sum(is.na(X_raw))/length(X_raw), 2), "%)")

# =========================================================================
# Phase 2: Stratified 5×5 CV with three classifiers
# =========================================================================

# --- helpers --------------------------------------------------------------
stratified_folds <- function(y, k = 5, repeats = 5, seed = 42) {
  set.seed(seed)
  folds <- list()
  for (r in seq_len(repeats)) {
    by_class <- split(seq_along(y), y)
    fold_id <- integer(length(y))
    for (cls in names(by_class)) {
      idx <- by_class[[cls]]
      shuf <- sample(idx)
      fold_id[shuf] <- rep_len(seq_len(k), length(shuf))
    }
    for (f in seq_len(k)) {
      folds[[length(folds) + 1]] <- list(
        repeat_id = r, fold_id = f,
        test_idx = which(fold_id == f),
        train_idx = which(fold_id != f)
      )
    }
  }
  folds
}

# Multi-class Matthews correlation (Gorodkin 2004)
mcc_multiclass <- function(C) {
  N <- sum(C)
  if (N == 0) return(NA_real_)
  t_k <- rowSums(C); p_k <- colSums(C); c_k <- sum(diag(C))
  num <- c_k * N - sum(t_k * p_k)
  den <- sqrt((N^2 - sum(p_k^2)) * (N^2 - sum(t_k^2)))
  if (den == 0) return(0)
  num / den
}

per_class_sens_spec <- function(C) {
  k <- nrow(C)
  out <- list(sens = setNames(numeric(k), rownames(C)),
              spec = setNames(numeric(k), rownames(C)))
  for (i in seq_len(k)) {
    TP <- C[i, i]; FN <- sum(C[i, -i]); FP <- sum(C[-i, i])
    TN <- sum(C) - TP - FN - FP
    out$sens[i] <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    out$spec[i] <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  }
  out
}

macro_roc_auc <- function(y_true, prob_mat) {
  # one-vs-rest macro-averaged AUC; prob_mat: rows donors, cols classes
  classes <- colnames(prob_mat)
  aucs <- numeric(length(classes))
  for (i in seq_along(classes)) {
    cls <- classes[i]
    bin_true <- as.integer(y_true == cls)
    aucs[i] <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(bin_true, prob_mat[, cls],
                                     quiet = TRUE, direction = "<"))),
      error = function(e) NA_real_)
  }
  mean(aucs, na.rm = TRUE)
}

# --- preprocessing per fold ----------------------------------------------
fit_preproc <- function(X_train) {
  med <- apply(X_train, 2, median, na.rm = TRUE)
  X_imp <- X_train
  for (j in seq_len(ncol(X_imp))) {
    nas <- is.na(X_imp[, j]); if (any(nas)) X_imp[nas, j] <- med[j]
  }
  mn <- colMeans(X_imp); sd <- apply(X_imp, 2, sd)
  sd[sd == 0] <- 1
  list(median = med, mean = mn, sd = sd)
}
apply_preproc <- function(X, pp) {
  X_imp <- X
  for (j in seq_len(ncol(X_imp))) {
    nas <- is.na(X_imp[, j]); if (any(nas)) X_imp[nas, j] <- pp$median[j]
  }
  sweep(sweep(X_imp, 2, pp$mean, "-"), 2, pp$sd, "/")
}

# --- classifier wrappers --------------------------------------------------
fit_rf <- function(X, y, class_weights) {
  ranger(y = y, x = X, num.trees = 1000, classification = TRUE,
         probability = TRUE,
         class.weights = class_weights,
         importance = "permutation",
         seed = 42, num.threads = 4)
}
predict_rf <- function(model, X) predict(model, data = X)$predictions

fit_xgb <- function(X, y, class_weights) {
  classes <- levels(y)
  y_int <- as.integer(y) - 1L
  # No per-sample weights — XGBoost's gradient boosting handles the 113/164/14
  # class imbalance natively. Earlier attempts with inverse-frequency weights
  # (full or sqrt-scaled) destabilised training and yielded MCC ~0 on test.
  # Drop colsample_bytree (use all 1250 features), reduce eta, add subsample,
  # cap depth at 6, and rely on regularisation + sufficient rounds.
  d <- xgb.DMatrix(data = X, label = y_int)
  params <- list(objective = "multi:softprob", num_class = length(classes),
                 eta = 0.05, max_depth = 6, subsample = 0.8,
                 min_child_weight = 1, gamma = 0,
                 eval_metric = "mlogloss",
                 nthread = 4)
  set.seed(42)
  model <- xgb.train(params = params, data = d, nrounds = 300, verbose = 0)
  attr(model, "classes") <- classes
  model
}
predict_xgb <- function(model, X) {
  classes <- attr(model, "classes")
  d <- xgb.DMatrix(data = X)
  probs <- predict(model, d)
  # xgboost 3.x returns an n × k matrix directly for multi:softprob;
  # older versions returned a flat n*k vector. Handle both.
  if (!is.matrix(probs)) {
    probs <- matrix(probs, ncol = length(classes), byrow = TRUE)
  }
  colnames(probs) <- classes
  probs
}

fit_glmnet <- function(X, y, class_weights) {
  w <- as.numeric(class_weights[as.character(y)])
  set.seed(42)
  cv <- cv.glmnet(x = X, y = y, family = "multinomial", alpha = 0.5,
                  weights = w, nfolds = 5, parallel = FALSE)
  cv
}
predict_glmnet <- function(model, X) {
  probs <- predict(model, newx = X, s = "lambda.min", type = "response")
  probs <- probs[, , 1]  # drop the singleton lambda dim
  probs
}

# --- run CV --------------------------------------------------------------
y <- factor(pheno$condition, levels = c("healthy", "sle", "sjs"))
inv_freq <- 1 / (table(y) / length(y))
inv_freq <- inv_freq / sum(inv_freq)
log_msg("Inverse-frequency class weights: ",
        paste(names(inv_freq), round(inv_freq, 3), sep="=", collapse=", "))

folds <- stratified_folds(y, k = 5, repeats = 5, seed = 42)
log_msg("CV folds: ", length(folds))

FOLD_MODELS_DIR <- file.path(OUTDIR, "fold_models")
dir.create(FOLD_MODELS_DIR, recursive = TRUE, showWarnings = FALSE)

models <- list(rf = list(fit = fit_rf, pred = predict_rf),
               xgb = list(fit = fit_xgb, pred = predict_xgb),
               glmnet = list(fit = fit_glmnet, pred = predict_glmnet))

all_preds <- list()
all_metrics <- list()
fold_importances <- list()

for (mname in names(models)) {
  log_msg("=== ", mname, " ===")
  t0 <- Sys.time()
  per_fold_metrics <- list()
  for (f_idx in seq_along(folds)) {
    f <- folds[[f_idx]]
    X_tr <- X_raw[f$train_idx, , drop = FALSE]
    X_te <- X_raw[f$test_idx,  , drop = FALSE]
    y_tr <- y[f$train_idx]; y_te <- y[f$test_idx]

    pp <- fit_preproc(X_tr)
    X_tr <- apply_preproc(X_tr, pp)
    X_te <- apply_preproc(X_te, pp)

    model <- models[[mname]]$fit(X_tr, y_tr, inv_freq)
    saveRDS(
      list(model = model, pp = pp,
           train_idx = f$train_idx, test_idx = f$test_idx,
           donor_ids_test = rownames(X_te),
           y_test = as.character(y_te)),
      file.path(FOLD_MODELS_DIR,
                sprintf("%s_rep%d_fold%d.rds", mname, f$repeat_id, f$fold_id)))
    probs <- models[[mname]]$pred(model, X_te)
    pred  <- factor(colnames(probs)[max.col(probs, ties.method = "first")],
                    levels = levels(y))

    # capture per-fold predictions
    all_preds[[length(all_preds) + 1]] <- data.frame(
      classifier = mname,
      repeat_id = f$repeat_id, fold_id = f$fold_id,
      donor_id = rownames(X_te),
      true = as.character(y_te),
      pred = as.character(pred),
      prob_healthy = probs[, "healthy"],
      prob_sle     = probs[, "sle"],
      prob_sjs     = probs[, "sjs"],
      stringsAsFactors = FALSE)

    C <- table(factor(y_te, levels = levels(y)),
               factor(pred, levels = levels(y)))
    sens_spec <- per_class_sens_spec(C)
    per_fold_metrics[[f_idx]] <- data.frame(
      classifier = mname, repeat_id = f$repeat_id, fold_id = f$fold_id,
      MCC = mcc_multiclass(C),
      AUC = macro_roc_auc(as.character(y_te), probs),
      sens_macro = mean(sens_spec$sens, na.rm = TRUE),
      spec_macro = mean(sens_spec$spec, na.rm = TRUE),
      sens_healthy = sens_spec$sens["healthy"],
      sens_sle     = sens_spec$sens["sle"],
      sens_sjs     = sens_spec$sens["sjs"],
      spec_healthy = sens_spec$spec["healthy"],
      spec_sle     = sens_spec$spec["sle"],
      spec_sjs     = sens_spec$spec["sjs"]
    )

    # importance: collect on full-train fits later (saves time per-fold);
    # here keep importance per-fold for averaging
    imp <- switch(mname,
      rf = importance(model),
      xgb = {
        m <- xgb.importance(model = model)
        setNames(m$Gain, m$Feature)
      },
      glmnet = {
        coefs <- coef(model, s = "lambda.min")
        # multi-class glmnet returns list of sparse vectors per class
        abs_sum <- Reduce("+", lapply(coefs, function(b) abs(as.matrix(b)[-1, 1])))
        setNames(abs_sum, names(abs_sum))
      })
    fold_importances[[length(fold_importances) + 1]] <-
      data.frame(classifier = mname, feature = names(imp), value = unname(imp))
  }
  all_metrics[[mname]] <- do.call(rbind, per_fold_metrics)
  log_msg("  done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1),
          " min")
}

cv_preds <- do.call(rbind, all_preds)
cv_metrics <- do.call(rbind, all_metrics)
write.table(cv_preds, file.path(OUTDIR, "cv_predictions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Headline benchmark: macro-averaged metrics only
summary_metrics <- cv_metrics %>%
  group_by(classifier) %>%
  summarise(across(c(MCC, AUC, sens_macro, spec_macro),
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"), .groups = "drop")
write.table(summary_metrics, file.path(OUTDIR, "benchmark.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Supplementary: per-class sensitivity and specificity (mean ± sd)
per_class_metrics <- cv_metrics %>%
  group_by(classifier) %>%
  summarise(across(c(sens_healthy, sens_sle, sens_sjs,
                     spec_healthy, spec_sle, spec_sjs),
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"), .groups = "drop")
write.table(per_class_metrics,
            file.path(OUTDIR, "benchmark_per_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Feature importance — average across folds, top 30 per classifier
imp_df <- do.call(rbind, fold_importances) %>%
  group_by(classifier, feature) %>%
  summarise(value_mean = mean(value, na.rm = TRUE),
            value_sd   = sd(value, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(classifier) %>%
  slice_max(value_mean, n = 30) %>%
  arrange(classifier, desc(value_mean))
write.table(imp_df, file.path(OUTDIR, "feature_importance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Confusion matrices from concatenated predictions
cm_list <- lapply(unique(cv_preds$classifier), function(cls) {
  d <- cv_preds[cv_preds$classifier == cls, ]
  C <- table(true = factor(d$true, levels = levels(y)),
             pred = factor(d$pred, levels = levels(y)))
  Cn <- sweep(C, 1, rowSums(C), "/")  # row-normalised
  df <- as.data.frame(Cn)
  colnames(df) <- c("true", "pred", "frac")
  df$classifier <- cls
  df
})
cm_df <- do.call(rbind, cm_list)

plt <- ggplot(cm_df, aes(x = pred, y = true, fill = frac)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
  facet_wrap(~ classifier, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#3182bd",
                      limits = c(0, 1), name = "row-norm") +
  scale_y_discrete(limits = rev) +
  labs(x = "predicted", y = "true",
       title = "Confusion matrices (concatenated 25-fold predictions, row-normalised)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "confusion_matrices.pdf"),
       plt, width = 11, height = 4)

log_msg("=== Done ===")
log_msg("Outputs in: ", OUTDIR)
