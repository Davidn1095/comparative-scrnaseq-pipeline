#!/usr/bin/env Rscript
################################################################################
# 08_donor_classifiers.R — Donor-level Hallmark-pathway classifiers
#
# Each donor is described by 1,250 features: the mean Seurat AddModuleScore of one
# MSigDB Hallmark pathway over the donor's cells of one of 25 cell types, on
# expression batch-corrected per cell type. Three classifiers (ranger RF, XGBoost,
# glmnet multinomial elastic net) are trained under 5×5 CV (25 folds) stratified by
# class and dataset jointly.
#
# Outputs (results/donor_classifiers_<method>/):
#   donor_expression_sums.rds   per cell type, gene × donor expression sums (cache)
#   feature_matrix.rds          X_folds, the 25 per-fold feature matrices (291 × 1250),
#                               donor metadata (pheno), the folds and the cases left
#                               uncorrected
#   cv_predictions.tsv          per-fold per-donor true + predicted + prob
#   benchmark.tsv               mean ± sd metrics per classifier
#   feature_importance.tsv      top-30 features per classifier
#   confusion_matrices.pdf      3-panel CM (concatenated CV predictions)
#   fold_models/                per-fold models, preprocessing and held-out features
#
# Leak protection: the batch correction is fitted inside each fold on the training
# donors and applied to the held-out donors without their labels (Phase 1), and
# median imputation and z-score normalisation are fitted on the training partition
# only, then applied to its test.
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
# Phase 1: Features, from per-donor expression sums (cached on disk)
# =========================================================================
# limma / limma_modulescore. The dataset effect is removed per cell type as
# limma::removeBatchEffect removes it: least squares with condition (treatment
# contrasts) and dataset (sum-to-zero contrasts) in the design, the dataset
# coefficients then subtracted. Each Hallmark score is averaged over a donor's cells
# of the cell type. The fit uses the donors' labels, so it is made inside every
# cross-validation fold from the training donors only and applied to all donors, the
# held-out ones without their labels (Phase 2). Fitting it once on all 291 donors
# before cross-validation let the held-out donors' own labels into their features.
#
# Both the correction and a donor's mean score are linear in the cells' expression,
# so each fold's features are computed exactly from per-donor expression sums per cell
# type, built once here: least squares on the cells equals least squares on the donor
# means weighted by cells per donor. With limma_modulescore, AddModuleScore's
# control-gene draw (Seurat 5.4.0, AddModuleScore.Assay; nbin 25, ctrl 100, seed 42) is
# replayed from the genes' mean corrected expression over all cells of the cell type.
# Genes whose means tie exactly are ordered by floating-point noise, so a few control
# sets can differ from a per-cell AddModuleScore run.
#
# none. AddModuleScore on the uncorrected expression; no fit, so one matrix serves
# every fold.
#
# The cache (the donor sums, or for none the uncorrected matrix) is derived from the
# atlas, and reusing it silently has bitten twice. After the Perez re-ingestion a rerun
# returned performance identical to four decimal places because it was the cache; on
# 2026-08-12 a stale-cache banner printed exactly as designed and the run trained every
# classifier on the stale features anyway, because warning and proceeding is still
# proceeding. The cache is reused only while it is NEWER than the atlas, and rebuilt
# otherwise. FORCE_FEATURES=1 also forces a rebuild.
force_features <- tolower(Sys.getenv("FORCE_FEATURES", "0")) %in% c("1", "true", "yes")
ATLAS_RDS  <- file.path(OBJ_DIR, "atlas_integrated.rds")
CACHE_PATH <- file.path(OUTDIR, if (DENOISE_METHOD == "none") "features_uncorrected.rds"
                                else "donor_expression_sums.rds")
cache_mtime <- if (file.exists(CACHE_PATH)) file.info(CACHE_PATH)$mtime else NA
atlas_mtime <- if (file.exists(ATLAS_RDS)) file.info(ATLAS_RDS)$mtime else NA
stamp <- function(x) if (is.na(x)) "missing" else format(x, "%Y-%m-%d %H:%M:%S")
cache_stale <- !is.na(cache_mtime) && !is.na(atlas_mtime) && atlas_mtime > cache_mtime

if (file.exists(CACHE_PATH) && force_features)
  log_msg("FORCE_FEATURES=1: rebuilding, ignoring ", CACHE_PATH)
if (cache_stale && !force_features) {
  log_msg(strrep("!", 74))
  log_msg("!! CACHED FEATURE INPUT IS STALE - rebuilding from the current atlas")
  log_msg("!!   cache built ", stamp(cache_mtime))
  log_msg("!!   atlas built ", stamp(atlas_mtime), "   <-- newer")
  log_msg("!! The cache is derived from the atlas, so it cannot be reused. Training on")
  log_msg("!! it would make every classifier metric, and all SHAP downstream of it,")
  log_msg("!! silently describe a superseded atlas.")
  log_msg(strrep("!", 74))
}
# The sums also depend on the gene exclusion (the gene pool drops excluded genes), which
# the atlas date cannot see: the cache records the exclusion it was built with and is
# rebuilt when EXCLUDED_GENES changes.
EXCLUSION_KEY <- paste(sort(EXCLUDED_GENES), collapse = ",")
cache <- NULL
if (file.exists(CACHE_PATH) && !force_features && !cache_stale) {
  cache <- readRDS(CACHE_PATH)
  if (identical(cache$exclusion, EXCLUSION_KEY)) {
    log_msg("Reusing cached feature input (", stamp(cache_mtime), ", newer than atlas ",
            stamp(atlas_mtime), ", same gene exclusion): ", CACHE_PATH)
  } else {
    log_msg("Cached feature input was built with a different gene exclusion - rebuilding")
    cache <- NULL
  }
}
if (is.null(cache)) {
  log_msg("Loading atlas (~2.5 min)...")
  obj <- readRDS(ATLAS_RDS)
  log_msg("Loading Hallmark gene sets...")
  hallmark_list <- load_msigdb_hallmark(MSIGDB_RDS)

  if (DENOISE_METHOD == "none") {
    hl <- lapply(hallmark_list, function(g) intersect(g, rownames(obj)))
    stopifnot(all(lengths(hl) >= 5))
    log_msg("AddModuleScore on 50 Hallmark sets (15-25 min)...")
    obj <- AddModuleScore(obj, features = hl, name = "HM_", nbin = 25, ctrl = 100, seed = 42)
    score_cols <- paste0("HM_", seq_along(hl))
    stopifnot(all(score_cols %in% colnames(obj@meta.data)))
    meta <- obj@meta.data
    rm(obj); gc()
    colnames(meta)[match(score_cols, colnames(meta))] <- names(hl)
    meta <- meta[is_common_celltype(meta$cell_type), ]
    agg <- meta %>%
      group_by(donor_id, cell_type) %>%
      summarise(across(all_of(names(hl)), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
    first <- meta[!duplicated(meta$donor_id), ]
    cache <- list(scores = as.data.frame(agg), hallmark_list = hl,
                  donor = data.frame(donor_id = first$donor_id, condition = first$condition,
                                     dataset_id = first$dataset_id, stringsAsFactors = FALSE))
  } else {
    USE_MODULESCORE <- (DENOISE_METHOD == "limma_modulescore")
    # Gene pool.
    #   limma:             HVGs (signal-carrying genes for the simple-mean path).
    #   limma_modulescore: union of all 50 Hallmark pathway genes. Aligns the gene pool
    #                      with the feature definition and gives AddModuleScore enough
    #                      genes for nbin=25 binning (~4000 / 25 ≈ 160 per bin, well
    #                      above ctrl=100).
    genes <- if (USE_MODULESCORE) unique(unlist(hallmark_list)) else Seurat::VariableFeatures(obj)
    genes <- genes[!is_excluded_gene(genes)]
    genes <- intersect(genes, rownames(obj))
    stopifnot(!any(grepl("_", genes)))   # a per-cell Seurat object would rename them
    log_msg("  ", if (USE_MODULESCORE) "Hallmark gene union: " else "HVG set: ", length(genes),
            " genes (after the gene exclusion, intersected with atlas)")
    expr <- LayerData(obj, layer = "data")[genes, , drop = FALSE]
    meta <- obj@meta.data
    rm(obj); gc()
    meta <- meta[is_common_celltype(meta$cell_type), ]
    log_msg("  Atlas cells in 25 common cell types: ", nrow(meta))
    col_idx <- match(rownames(meta), colnames(expr))
    stopifnot(!anyNA(col_idx))
    sums <- list()
    for (ct in unique(as.character(meta$cell_type))) {
      sel <- which(meta$cell_type == ct)
      d <- factor(meta$donor_id[sel], levels = unique(meta$donor_id[sel]))
      D <- Matrix::sparseMatrix(i = seq_along(d), j = as.integer(d), x = 1,
                                dims = c(length(d), nlevels(d)))
      S <- as.matrix(expr[, col_idx[sel], drop = FALSE] %*% D)
      dimnames(S) <- list(genes, levels(d))
      first <- sel[match(levels(d), meta$donor_id[sel])]
      sums[[ct]] <- list(S = S, n = tabulate(as.integer(d), nlevels(d)),
                         condition = as.character(meta$condition[first]),
                         dataset = as.character(meta$dataset_id[first]))
      log_msg("    ", ct, ": ", nlevels(d), " donors, ", length(sel), " cells")
    }
    cache <- list(genes = genes, hallmark_list = hallmark_list, sums = sums,
                  score = if (USE_MODULESCORE) "modulescore" else "mean")
    rm(expr, meta); gc()
  }
  cache$exclusion <- EXCLUSION_KEY
  saveRDS(cache, CACHE_PATH)
  log_msg("Saved: ", CACHE_PATH)
}

pathways <- names(cache$hallmark_list)
if (DENOISE_METHOD == "none") {
  long <- cache$scores %>%
    pivot_longer(cols = all_of(pathways), names_to = "pathway", values_to = "score") %>%
    mutate(feature = paste0(pathway, "__", gsub("[-/ ]", "_", cell_type)))
  wide <- long %>% select(donor_id, feature, score) %>%
    pivot_wider(names_from = feature, values_from = score)
  X_fixed <- as.matrix(wide[, -1]); rownames(X_fixed) <- wide$donor_id
  donors <- sort(rownames(X_fixed)); feats <- colnames(X_fixed)
  X_fixed <- X_fixed[donors, , drop = FALSE]
  pheno <- cache$donor[match(donors, cache$donor$donor_id), ]
  n_genes <- NA_integer_
} else {
  # Donor-level feature matrix for one fold: dataset effect fitted on train_donors
  gene_sets <- lapply(cache$hallmark_list, function(g) intersect(g, cache$genes))
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  cts <- names(cache$sums)
  feats <- as.vector(t(outer(cts, pathways, function(ct, pw) paste0(pw, "__", gsub("[-/ ]", "_", ct)))))
  donors <- sort(unique(unlist(lapply(cache$sums, function(s) colnames(s$S)))))
  pheno <- do.call(rbind, lapply(cache$sums, function(s)
    data.frame(donor_id = colnames(s$S), condition = s$condition, dataset_id = s$dataset,
               stringsAsFactors = FALSE)))
  pheno <- pheno[match(donors, pheno$donor_id), ]
  n_genes <- length(cache$genes)
}
rownames(pheno) <- pheno$donor_id

# Dataset effect per gene (genes × datasets), fitted on the donors in `use` and 0 for a
# dataset without any of them, whose effect cannot be estimated.
dataset_effect <- function(s, use) {
  datasets <- sort(unique(s$dataset))
  eff <- matrix(0, nrow(s$S), length(datasets), dimnames = list(rownames(s$S), datasets))
  batch <- factor(s$dataset[use])
  if (nlevels(batch) < 2) return(eff)
  contrasts(batch) <- contr.sum(levels(batch))
  X_batch <- model.matrix(~batch)[, -1, drop = FALSE]
  cond <- factor(s$condition[use])
  design <- if (nlevels(cond) > 1) model.matrix(~cond) else matrix(1, length(cond), 1)
  fit <- lm.wfit(cbind(design, X_batch), t(s$S[, use, drop = FALSE]) / s$n[use], w = s$n[use])
  beta <- fit$coefficients[-seq_len(ncol(design)), , drop = FALSE]
  beta[is.na(beta)] <- 0
  eff[, levels(batch)] <- t(beta) %*% t(contr.sum(levels(batch)))
  eff
}

# AddModuleScore's control genes, replayed from the genes' mean expression
control_genes <- function(avg, gene_sets, nbin = 25, ctrl = 100, seed = 42) {
  set.seed(seed)
  avg <- avg[order(avg)]
  cut <- ggplot2::cut_number(x = avg + rnorm(n = length(avg)) / 1e+30, n = nbin,
                             labels = FALSE, right = FALSE)
  names(cut) <- names(avg)
  lapply(gene_sets, function(f) unique(unlist(lapply(f, function(g)
    names(sample(x = cut[which(cut == cut[g])], size = ctrl, replace = FALSE))))))
}

# Donor means of the per-cell scores of one cell type after subtracting `eff`
donor_scores <- function(s, eff) {
  shift <- eff[, s$dataset, drop = FALSE]
  M <- sweep(s$S, 2, s$n, "/") - shift
  if (cache$score == "modulescore") {
    avg <- (rowSums(s$S) - drop(shift %*% s$n)) / sum(s$n)
    ctrl <- control_genes(avg, gene_sets)
  }
  sc <- vapply(seq_along(gene_sets), function(i) {
    v <- colMeans(M[gene_sets[[i]], , drop = FALSE])
    if (cache$score == "modulescore") v - colMeans(M[ctrl[[i]], , drop = FALSE]) else v
  }, numeric(ncol(M)))
  if (!is.matrix(sc)) sc <- matrix(sc, nrow = 1)
  dimnames(sc) <- list(colnames(s$S), names(gene_sets))
  sc
}

fold_features <- function(train_donors) {
  if (DENOISE_METHOD == "none") return(structure(X_fixed, uncorrected = character(0)))
  X <- matrix(NA_real_, length(donors), length(feats), dimnames = list(donors, feats))
  uncorrected <- character(0)
  for (ct in names(cache$sums)) {
    s <- cache$sums[[ct]]
    use <- colnames(s$S) %in% train_donors
    gone <- setdiff(unique(s$dataset), s$dataset[use])
    if (length(gone)) uncorrected <- c(uncorrected, paste0(ct, ": ", paste(gone, collapse = ", ")))
    sc <- donor_scores(s, dataset_effect(s, use))
    X[rownames(sc), paste0(colnames(sc), "__", gsub("[-/ ]", "_", ct))] <- sc
  }
  attr(X, "uncorrected") <- uncorrected
  X
}

log_msg("Features: ", length(donors), " donors × ", length(feats), " features. Conditions: ",
        paste(names(table(pheno$condition)), table(pheno$condition), sep = "=", collapse = ", "))

# =========================================================================
# Phase 2: 5×5 CV, stratified by class and dataset, with three classifiers
# =========================================================================

# --- helpers --------------------------------------------------------------
# Folds stratified by class and dataset jointly: donors are ordered by dataset and
# class, shuffled within each class-by-dataset stratum and dealt to the k folds in
# turn from a random starting fold. A dataset's donors thus fall into as many folds as
# it has donors (up to k), so no fold holds out every donor of a dataset with two or
# more, and every fold keeps the class balance.
stratified_folds <- function(y, group, k = 5, repeats = 5, seed = 42) {
  set.seed(seed)
  key <- paste(group, y, sep = "\r")
  strata <- split(seq_along(y), factor(key, levels = sort(unique(key))))
  folds <- list()
  for (r in seq_len(repeats)) {
    ord <- unlist(lapply(strata, function(i) i[sample.int(length(i))]), use.names = FALSE)
    fold_id <- integer(length(y))
    fold_id[ord] <- (seq_along(ord) + sample.int(k, 1) - 2) %% k + 1
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

folds <- stratified_folds(y, pheno$dataset_id, k = 5, repeats = 5, seed = 42)
log_msg("CV folds: ", length(folds))
no_train <- unlist(lapply(folds, function(f)
  setdiff(unique(pheno$dataset_id), pheno$dataset_id[f$train_idx])))
if (length(no_train)) stop("A fold has no training donor of dataset(s) ",
                           paste(unique(no_train), collapse = ", "), call. = FALSE)
log_msg("Every fold has training donors of all ", length(unique(pheno$dataset_id)), " datasets")

# Per-fold features: batch correction fitted on each fold's training donors
log_msg("Building the 25 per-fold feature matrices...")
X_folds <- lapply(folds, function(f) fold_features(donors[f$train_idx]))
names(X_folds) <- vapply(folds, function(f) sprintf("rep%d_fold%d", f$repeat_id, f$fold_id), "")
uncorrected <- do.call(rbind, lapply(names(X_folds), function(nm) {
  u <- attr(X_folds[[nm]], "uncorrected")
  if (length(u)) data.frame(fold = nm, cell_type_dataset = u) }))
if (is.null(uncorrected)) {
  log_msg("Every cell type had training donors of every dataset in every fold")
} else {
  log_msg("Cell type × dataset pairs without training donors in a fold (left uncorrected): ",
          nrow(uncorrected), " - ", paste(uncorrected$fold, uncorrected$cell_type_dataset, collapse = "; "))
}
X_folds <- lapply(X_folds, function(X) { attr(X, "uncorrected") <- NULL; X })
X1 <- X_folds[[1]]
log_msg("NA count (donors without cells of a cell type): ", sum(is.na(X1)),
        " (", round(100 * sum(is.na(X1)) / length(X1), 2), "%)")
fm <- list(X_folds = X_folds, pheno = pheno, feature_names = feats, folds = folds,
           uncorrected = uncorrected, n_genes = n_genes)
saveRDS(fm, FEATURE_MATRIX_PATH)
log_msg("Saved: ", FEATURE_MATRIX_PATH)

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
    X_tr <- X_folds[[f_idx]][f$train_idx, , drop = FALSE]
    X_te <- X_folds[[f_idx]][f$test_idx,  , drop = FALSE]
    X_test_raw <- X_te
    y_tr <- y[f$train_idx]; y_te <- y[f$test_idx]

    pp <- fit_preproc(X_tr)
    X_tr <- apply_preproc(X_tr, pp)
    X_te <- apply_preproc(X_te, pp)

    model <- models[[mname]]$fit(X_tr, y_tr, inv_freq)
    saveRDS(
      list(model = model, pp = pp,
           train_idx = f$train_idx, test_idx = f$test_idx,
           donor_ids_test = rownames(X_te),
           y_test = as.character(y_te), X_test_raw = X_test_raw),
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
