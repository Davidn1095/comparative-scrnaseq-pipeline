# R packages of atlas-final-2026.10: the same top-level installs as Dockerfile.atlas (atlas-2026.10), with
# CRAN pinned to the Posit Package Manager snapshot of 2026-10-22 (atlas-2026.10 used .../latest on
# 2026-10-23) and Bioconductor 3.23 (whose 77 packages in atlas-2026.10c were still the current 3.23
# versions on 2026-10-28). Afterwards every package of the lock (locks/atlas-2026.10c.R.tsv) must be
# installed at its locked version; a CRAN package that is not is installed at that version from the CRAN
# archive, and anything else that differs stops the build.
options(Ncpus = 4L, warn = 1, repos = c(CRAN = "https://p3m.dev/cran/__linux__/noble/2026-10-22"))
stopifnot(as.character(BiocManager::version()) == "3.23")

cran <- c("Seurat", "SeuratObject", "Matrix", "hdf5r", "reticulate", "pROC", "future",
          "glmnet", "ranger", "xgboost", "openxlsx2", "msigdbr",
          "ggplot2", "cowplot", "patchwork", "ggrepel", "ggrastr", "circlize", "magick",
          "png", "scales", "gtable", "dplyr", "tidyr", "readr", "tibble")
bioc <- c("SingleR", "celldex", "scDblFinder", "decontX", "DESeq2", "fgsea", "limma",
          "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "SingleCellExperiment",
          "SummarizedExperiment", "S4Vectors", "BiocParallel", "ExperimentHub",
          "AnnotationHub", "ComplexHeatmap")
BiocManager::install(c(cran, bioc), ask = FALSE, update = FALSE)
BiocManager::install("decontX", type = "source", force = TRUE, update = FALSE, ask = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("hrbrmstr/ggchicklet@64c468dd0900153be1690dbfc5cfb35710da8183", upgrade = "never")

lock <- read.delim("/opt/atlas-locks/atlas-2026.10c.R.tsv", colClasses = "character")
have <- function() {
  ip <- installed.packages()
  setNames(ip[, "Version"], ip[, "Package"])
}
off <- function() {
  h <- have()
  lock[is.na(h[lock$Package]) | h[lock$Package] != lock$Version, ]
}
o <- off()
cran_fix <- o[o$Repository %in% c("CRAN", "RSPM"), ]
for (k in seq_len(nrow(cran_fix))) {
  message("installing locked CRAN version ", cran_fix$Package[k], " ", cran_fix$Version[k])
  remotes::install_version(cran_fix$Package[k], version = cran_fix$Version[k],
                           repos = "https://cloud.r-project.org", upgrade = "never", dependencies = FALSE)
}
o <- off()
if (nrow(o)) {
  h <- have()
  stop("packages not at their locked version: ",
       paste0(o$Package, " ", o$Version, " (installed ", h[o$Package], ")", collapse = "; "))
}
extra <- setdiff(names(have()), lock$Package)
if (length(extra)) stop("packages not in the lock: ", paste(extra, collapse = ", "))
library(decontX)
cat("R packages: all", nrow(lock), "at their locked versions\n")
