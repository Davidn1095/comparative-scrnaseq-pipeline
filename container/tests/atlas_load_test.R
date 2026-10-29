# Load every R package the atlas uses; print one "PKG <name> <version>" line each.
# Exits non-zero if any package fails to load or a required function/argument is missing.
pkgs <- c("Seurat", "SeuratObject", "Matrix", "hdf5r", "reticulate", "pROC", "future",
          "glmnet", "ranger", "xgboost", "openxlsx2", "msigdbr",
          "ggplot2", "cowplot", "patchwork", "ggrepel", "ggrastr", "ggchicklet", "circlize",
          "magick", "png", "scales", "gtable", "dplyr", "tidyr", "readr", "tibble",
          "SingleR", "celldex", "scDblFinder", "decontX", "DESeq2", "fgsea", "limma",
          "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "SingleCellExperiment",
          "SummarizedExperiment", "S4Vectors", "BiocParallel", "ExperimentHub",
          "AnnotationHub", "ComplexHeatmap", "BiocManager")
fail <- character()
for (p in pkgs) {
  ok <- suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE))
  if (ok) cat("PKG", p, as.character(packageVersion(p)), "\n") else { cat("FAIL", p, "\n"); fail <- c(fail, p) }
}
cat("PKG R", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("PKG Bioconductor", as.character(BiocManager::version()), "\n")
# Functions and arguments the decisions rely on
chk <- list(
  "scDblFinder has 'dbr'"             = "dbr" %in% names(formals(scDblFinder::scDblFinder)),
  "Seurat::AddModuleScore exists"     = is.function(Seurat::AddModuleScore),
  "fgsea::fgsea exists"               = is.function(fgsea::fgsea),
  "SingleR::SingleR exists"           = is.function(SingleR::SingleR),
  "DESeq2::DESeq exists"              = is.function(DESeq2::DESeq))
# decontX: run on a tiny simulated SingleCellExperiment with an empty-droplet background
chk[["decontX runs with background and bgBatch"]] <- tryCatch({
  set.seed(1)
  g <- 60; nc <- 200; nb <- 400
  prof <- matrix(rgamma(g * 2, 1), g, 2)
  z <- rep(1:2, each = nc / 2)
  cells <- sapply(z, function(k) rpois(g, prof[, k] * 5))
  bg <- matrix(rpois(g * nb, 0.3), g, nb)
  rownames(cells) <- rownames(bg) <- paste0("g", seq_len(g))
  colnames(cells) <- paste0("cell", seq_len(nc)); colnames(bg) <- paste0("empty", seq_len(nb))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = cells))
  bsce <- SingleCellExperiment::SingleCellExperiment(list(counts = bg))
  res <- decontX::decontX(sce, z = z, batch = rep(c("b1", "b2"), nc / 2),
                          background = bsce, bgBatch = rep(c("b1", "b2"), nb / 2),
                          maxIter = 20, verbose = FALSE)
  "decontXcounts" %in% SummarizedExperiment::assayNames(res)
}, error = function(e) { cat("decontX error:", conditionMessage(e), "\n"); FALSE })
for (n in names(chk)) { cat(if (isTRUE(chk[[n]])) "CHECK OK  " else "CHECK FAIL", n, "\n"); if (!isTRUE(chk[[n]])) fail <- c(fail, n) }
# reticulate must see the image's Python 3.12
py <- reticulate::py_config(); pyv <- as.character(py$version)
cat("CHECK reticulate python", py$python, pyv, "\n")
if (!startsWith(pyv, "3.12")) fail <- c(fail, "reticulate python")
if (length(fail)) { cat("FAILED:", paste(fail, collapse = ", "), "\n"); quit(status = 1) }
cat("ALL R LOAD TESTS PASSED\n")
