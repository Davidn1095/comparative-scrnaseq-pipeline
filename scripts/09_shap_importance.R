#!/usr/bin/env Rscript
################################################################################
# 09_shap_importance.R — Out-of-fold SHAP importance for a donor classifier
#
# Parameterised by the CLASSIFIER env var:
#   CLASSIFIER=xgb     exact TreeSHAP via predict(..., predcontrib = TRUE)
#   CLASSIFIER=glmnet  exact linear SHAP from elastic-net coefficients
#                      (phi_j(x) = beta_j * x_j; features are already z-scored
#                      to training-partition mean within each fold's pp block,
#                      so E[X_j] = 0 within-fold).
#
# Outputs (results/donor_classifiers_<method>/shap_<classifier>/):
#   shap_values_oof.rds                    291 × 1250 × 3 array (named dims)
#   shap_importance_per_feature.tsv        1250 rows
#   shap_importance_per_celltype.tsv       25 rows
#   shap_importance_per_pathway.tsv        50 rows
#   shap_importance_per_celltype_per_fold.tsv  25 cell types × 25 folds
#   shap_importance_per_pathway_per_fold.tsv   50 pathways × 25 folds
#   shap_importance_per_feature_per_class.tsv   1250 features × 3 classes (long)
#   shap_importance_per_celltype_per_class.tsv  25 cell types × 3 classes (long)
#   shap_importance_per_pathway_per_class.tsv   50 pathways × 3 classes (long)
#
# Inputs (results/donor_classifiers_<method>/):
#   feature_matrix.rds                     X_raw + donor metadata
#   fold_models/<classifier>_rep*_fold*.rds  per-fold trained models + pp
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

# CLASSIFIER IS MANDATORY - no default.
#
# This script computes SHAP for ONE classifier per invocation and writes to
# shap_<CLASSIFIER>/. It used to default to "glmnet", so a run that did not set it
# refreshed shap_glmnet/ and left shap_xgb/ untouched - producing a results tree that
# was half current and half months old, with nothing to distinguish the two. A
# comparison against shap_xgb/ then compared a stale file to itself and returned a
# Kendall tau of exactly +1.000, which is how this was caught. Both classifiers must
# be run explicitly, one invocation each.
CLASSIFIER <- tolower(Sys.getenv("CLASSIFIER", ""))
if (!nzchar(CLASSIFIER)) {
  stop("CLASSIFIER is not set. It must be xgb or glmnet, and it selects BOTH the SHAP ",
       "method (TreeSHAP vs linear) and the output tree shap_<CLASSIFIER>/. Run this ",
       "script once per classifier. Refusing to guess, because guessing leaves the ",
       "other tree stale and indistinguishable from a fresh one.")
}
stopifnot(CLASSIFIER %in% c("xgb", "glmnet"))

suppressPackageStartupMessages({
  if (CLASSIFIER == "xgb")    library(xgboost)
  if (CLASSIFIER == "glmnet") { library(glmnet); library(Matrix) }
})

OUTDIR <- switch(DENOISE_METHOD,
  scvi              = file.path(ATLAS_ROOT, "results", "donor_classifiers_scvi"),
  limma             = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma"),
  limma_modulescore = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma_modulescore"),
  file.path(ATLAS_ROOT, "results", "donor_classifiers"))

FOLD_MODELS_DIR <- file.path(OUTDIR, "fold_models")
SHAP_OUTDIR     <- file.path(OUTDIR, paste0("shap_", CLASSIFIER))
dir.create(SHAP_OUTDIR, recursive = TRUE, showWarnings = FALSE)
log_msg("Denoising method:  ", DENOISE_METHOD)
log_msg("Classifier:        ", CLASSIFIER,
        " (",
        if (CLASSIFIER == "xgb") "TreeSHAP" else "linear SHAP",
        ")")
log_msg("Reading fold models from: ", FOLD_MODELS_DIR)
log_msg("Writing SHAP to:          ", SHAP_OUTDIR)

fm <- readRDS(file.path(OUTDIR, "feature_matrix.rds"))
X_raw  <- fm$X_raw
donors <- rownames(X_raw)
feats  <- fm$feature_names
classes <- c("healthy", "sle", "sjs")
log_msg("Feature matrix: ", length(donors), " donors × ", length(feats), " features")

apply_preproc <- function(X, pp) {
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j]); if (any(nas)) X[nas, j] <- pp$median[j]
  }
  sweep(sweep(X, 2, pp$mean, "-"), 2, pp$sd, "/")
}

bundle_pattern <- paste0("^", CLASSIFIER, "_rep[0-9]+_fold[0-9]+\\.rds$")
bundles <- list.files(FOLD_MODELS_DIR, pattern = bundle_pattern, full.names = TRUE)
log_msg("Found ", length(bundles), " ", CLASSIFIER, " fold bundles")
stopifnot(length(bundles) == 25L)

# Each donor is held out exactly once per repeat × 5 repeats = 5 times.
shap_sum <- array(0,
                  dim = c(length(donors), length(feats), length(classes)),
                  dimnames = list(donors, feats, classes))
oof_count <- integer(length(donors)); names(oof_count) <- donors

feat_pathway  <- vapply(strsplit(feats, "__", fixed = TRUE), `[`, character(1), 1)
feat_celltype <- vapply(strsplit(feats, "__", fixed = TRUE), `[`, character(1), 2)

per_fold_ct_rows <- list()
per_fold_pw_rows <- list()

fold_shap <- function(b, classifier) {
  X_te <- apply_preproc(X_raw[b$test_idx, , drop = FALSE], b$pp)
  out <- array(0, dim = c(nrow(X_te), length(feats), length(classes)),
               dimnames = list(b$donor_ids_test, feats, classes))
  if (classifier == "xgb") {
    d_te <- xgb.DMatrix(data = X_te)
    contribs <- predict(b$model, d_te, predcontrib = TRUE)
    # xgboost 3.x returns n_test × n_classes × (n_features + 1) array.
    stopifnot(is.array(contribs), length(dim(contribs)) == 3L,
              dim(contribs)[2] == length(classes),
              dim(contribs)[3] == length(feats) + 1L)
    for (k in seq_along(classes)) {
      cm <- contribs[, k, seq_len(length(feats))]
      if (!is.matrix(cm)) cm <- matrix(cm, nrow = 1L)
      colnames(cm) <- feats
      out[, , k] <- cm
    }
  } else if (classifier == "glmnet") {
    # multinomial cv.glmnet at s = lambda.min: list of K sparse matrices
    cf <- coef(b$model, s = "lambda.min")
    stopifnot(is.list(cf), length(cf) == length(classes))
    for (k in seq_along(classes)) {
      bk <- as.numeric(cf[[classes[k]]])
      stopifnot(length(bk) == length(feats) + 1L)
      beta_k <- bk[-1L]
      cm <- sweep(X_te, 2, beta_k, "*")
      colnames(cm) <- feats
      out[, , k] <- cm
    }
  }
  out
}

for (b_path in bundles) {
  b <- readRDS(b_path)
  cm_arr <- fold_shap(b, CLASSIFIER)   # n_test × n_features × n_classes

  fold_abs_per_feat <- numeric(length(feats))
  names(fold_abs_per_feat) <- feats
  for (k in seq_along(classes)) {
    shap_sum[b$donor_ids_test, , k] <-
      shap_sum[b$donor_ids_test, , k] + cm_arr[, , k]
    fold_abs_per_feat <- fold_abs_per_feat + colMeans(abs(cm_arr[, , k]))
  }
  fold_abs_per_feat <- fold_abs_per_feat / length(classes)
  oof_count[b$donor_ids_test] <- oof_count[b$donor_ids_test] + 1L

  fname <- sub("\\.rds$", "", basename(b_path))
  m <- regmatches(fname,
                  regexec(paste0(CLASSIFIER, "_rep([0-9]+)_fold([0-9]+)"),
                          fname))[[1]]
  rep_id  <- as.integer(m[2])
  fold_id <- as.integer(m[3])

  per_fold_ct_rows[[length(per_fold_ct_rows) + 1L]] <- data.frame(
    classifier = CLASSIFIER, repeat_id = rep_id, fold_id = fold_id,
    celltype = names(tapply(fold_abs_per_feat, feat_celltype, mean)),
    mean_abs_shap = unname(tapply(fold_abs_per_feat, feat_celltype, mean)),
    stringsAsFactors = FALSE)
  per_fold_pw_rows[[length(per_fold_pw_rows) + 1L]] <- data.frame(
    classifier = CLASSIFIER, repeat_id = rep_id, fold_id = fold_id,
    pathway = names(tapply(fold_abs_per_feat, feat_pathway, mean)),
    mean_abs_shap = unname(tapply(fold_abs_per_feat, feat_pathway, mean)),
    stringsAsFactors = FALSE)
}
log_msg("OOF coverage: min=", min(oof_count), " max=", max(oof_count),
        " (expect 5 per donor)")
stopifnot(all(oof_count == 5L))

ct_per_fold <- do.call(rbind, per_fold_ct_rows)
pw_per_fold <- do.call(rbind, per_fold_pw_rows)
write.table(ct_per_fold,
            file.path(SHAP_OUTDIR, "shap_importance_per_celltype_per_fold.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pw_per_fold,
            file.path(SHAP_OUTDIR, "shap_importance_per_pathway_per_fold.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_celltype_per_fold.tsv (",
        nrow(ct_per_fold), " rows)")
log_msg("Wrote: shap_importance_per_pathway_per_fold.tsv (",
        nrow(pw_per_fold), " rows)")

shap_oof <- shap_sum / array(oof_count, dim = dim(shap_sum))
log_msg("SHAP array assembled: ",
        paste(dim(shap_oof), collapse = " × "))

saveRDS(list(shap = shap_oof, classes = classes,
             classifier = CLASSIFIER, denoise = DENOISE_METHOD,
             oof_count = oof_count),
        file.path(SHAP_OUTDIR, "shap_values_oof.rds"))
log_msg("Wrote: ", file.path(SHAP_OUTDIR, "shap_values_oof.rds"))

abs_shap <- abs(shap_oof)
per_feature <- apply(abs_shap, 2, mean)

per_feat_df <- data.frame(
  feature       = feats,
  pathway       = feat_pathway,
  celltype      = feat_celltype,
  mean_abs_shap = unname(per_feature),
  stringsAsFactors = FALSE
)
per_feat_df <- per_feat_df[order(-per_feat_df$mean_abs_shap), ]
write.table(per_feat_df,
            file.path(SHAP_OUTDIR, "shap_importance_per_feature.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_feature.tsv (", nrow(per_feat_df), " rows)")

per_ct <- aggregate(mean_abs_shap ~ celltype, per_feat_df, mean)
per_ct$n_features <- as.integer(table(per_feat_df$celltype)[per_ct$celltype])
colnames(per_ct)[2] <- "mean_mean_abs_shap"
per_ct <- per_ct[order(-per_ct$mean_mean_abs_shap), ]
write.table(per_ct,
            file.path(SHAP_OUTDIR, "shap_importance_per_celltype.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_celltype.tsv (", nrow(per_ct), " rows)")

per_pw <- aggregate(mean_abs_shap ~ pathway, per_feat_df, mean)
per_pw$n_celltypes <- as.integer(table(per_feat_df$pathway)[per_pw$pathway])
colnames(per_pw)[2] <- "mean_mean_abs_shap"
per_pw <- per_pw[order(-per_pw$mean_mean_abs_shap), ]
write.table(per_pw,
            file.path(SHAP_OUTDIR, "shap_importance_per_pathway.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_pathway.tsv (", nrow(per_pw), " rows)")

# ---- Per-class importance ----------------------------------------------------
# The same three aggregations without averaging across classes: for each
# feature, mean |SHAP| toward each class across all donors, then the mean of
# those values over the 50 pathways within each cell type and over the 25 cell
# types within each pathway. Long format, one row per (unit, class), rows in
# the class-averaged order above. Averaging the three class rows of any unit
# reproduces its class-averaged value. These feed Fig 4b-d (per-class stacked
# bars, per-class heatmaps) and Supplementary Table 7. No per-class fold SDs
# are written: they are large by construction for the 14-donor SjS class and
# are not reported anywhere.
per_feat_class <- do.call(rbind, lapply(seq_along(classes), function(k) data.frame(
  feature       = feats,
  pathway       = feat_pathway,
  celltype      = feat_celltype,
  class         = classes[k],
  mean_abs_shap = unname(colMeans(abs_shap[, , k])),
  stringsAsFactors = FALSE)))
per_feat_class <- per_feat_class[order(match(per_feat_class$feature, per_feat_df$feature),
                                       match(per_feat_class$class, classes)), ]
write.table(per_feat_class,
            file.path(SHAP_OUTDIR, "shap_importance_per_feature_per_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_feature_per_class.tsv (", nrow(per_feat_class), " rows)")

per_ct_class <- aggregate(mean_abs_shap ~ celltype + class, per_feat_class, mean)
colnames(per_ct_class)[3] <- "mean_mean_abs_shap"
per_ct_class <- per_ct_class[order(match(per_ct_class$celltype, per_ct$celltype),
                                   match(per_ct_class$class, classes)), ]
write.table(per_ct_class,
            file.path(SHAP_OUTDIR, "shap_importance_per_celltype_per_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_celltype_per_class.tsv (", nrow(per_ct_class), " rows)")

per_pw_class <- aggregate(mean_abs_shap ~ pathway + class, per_feat_class, mean)
colnames(per_pw_class)[3] <- "mean_mean_abs_shap"
per_pw_class <- per_pw_class[order(match(per_pw_class$pathway, per_pw$pathway),
                                   match(per_pw_class$class, classes)), ]
write.table(per_pw_class,
            file.path(SHAP_OUTDIR, "shap_importance_per_pathway_per_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote: shap_importance_per_pathway_per_class.tsv (", nrow(per_pw_class), " rows)")

log_msg("=== Done ===")
