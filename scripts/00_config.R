################################################################################
# 00_config.R - Centralized Configuration for the scTSig Atlas Pipeline
#
# This file defines all configurable parameters, thresholds, and paths used
# across the pipeline. Each R script sources this file at the top.
#
# Parameters can be overridden via environment variables.
################################################################################

# --- Suppress Rplots.pdf creation from stray implicit prints ---
if (!interactive()) pdf(NULL)

# --- Helper: read environment variable with typed default ---
get_env <- function(key, default = "") {
  val <- Sys.getenv(key, default)
  if (identical(val, "")) default else val
}

get_env_numeric <- function(key, default) {
 as.numeric(get_env(key, as.character(default)))
}

get_env_integer <- function(key, default) {
  as.integer(get_env(key, as.character(default)))
}

# --- Paths ---
# Repository root: set ATLAS_ROOT, or run the scripts from the repository root.
ATLAS_ROOT <- get_env("ATLAS_ROOT", getwd())

# --- Library path ---
lib_path <- Sys.getenv("R_LIBS_USER", "")
if (lib_path != "" && dir.exists(lib_path)) {
  .libPaths(c(lib_path, .libPaths()))
}

# --- Basilisk / reticulate ---
# Basilisk env is baked into the container (BASILISK_USE_SYSTEM_DIR=1 in Dockerfile).
# No external BASILISK_EXTERNAL_DIR or RETICULATE_CONDA needed.

# --- Per-dataset environment ---
ACC_DIR    <- get_env("ACC_DIR", getwd())
ACCESSION  <- trimws(get_env("ACCESSION", basename(ACC_DIR)))
DISEASE    <- get_env("DISEASE", basename(dirname(ACC_DIR)))

# --- QC thresholds (script 01) ---
MIN_FEATURES     <- get_env_integer("MIN_FEATURES", 200)
MAX_MT_PCT       <- get_env_numeric("MAX_MT_PCT", 30)
MAX_RIBO_PCT     <- get_env_numeric("MAX_RIBO_PCT", 60)
IQR_MULT         <- get_env_numeric("IQR_MULT", 1.5)
# --- Analysis parameters ---
N_VARIABLE_FEATURES <- get_env_integer("N_VARIABLE_FEATURES", 2000)

# Inclusion thresholds (used by both DESeq2 DE and singleDeep prep):
#   MIN_DONORS         = 8   per condition
#   MIN_CELLS          = 100 total cells per cell type (singleDeep default)
#   MIN_CELLS_CLASS    = 10  cells per condition class (singleDeep default)
#   MIN_CELLS_PER_DONOR= 0   no per-donor cell floor (follows singleDeep
#                            default, avoids artificial confounding)
#
# MIN_DONORS=8 is the minimum required by singleDeep's nested cross-validation
# (3 outer x 2 inner stratified folds need >=6 donors per class; +2 margin
# protects against unbalanced random splits causing zero-class divisions in
# inverse-frequency class weighting). Applying the strictest threshold to
# both DE and singleDeep ensures consistent cell-type coverage across all
# analyses (27 cell types). DE benchmarks (Crowell 2020 muscat) recommend
# >=3 donors per condition; 8 is conservative.
PADJ_THRESH         <- get_env_numeric("PADJ_THRESH", 0.05)
LFC_THRESH          <- get_env_numeric("LFC_THRESH", 1.0)

# Genes excluded from DE / GSEA / classifier feature space.
# - MT-* : mitochondrial transcripts (technical / stress correlate)
# - RPL*/RPS*: ribosomal protein genes (housekeeping noise)
# - HB[ABDEGMQZ] : hemoglobin chains and pseudogenes (red blood cell contamination)
# - SLC4A1, ALAS2, AHSP, CA1: further erythroid genes, the same red-cell contamination
#   that the haemoglobin pattern misses (SLC4A1 was the largest SLE fold change)
# - Y-chromosome + XIST/TSIX: sex-chromosome genes (confounded by cohort sex
#   imbalance; CXG SLE has 9.3% male donors vs CXG HC 2.0% → drives spurious
#   DDX3Y/EIF1AY upregulation in SLE that is not biology).
EXCLUDED_GENE_REGEX <- "^MT-|^RPL|^RPS|^HB[ABDEGMQZ]|^(SLC4A1|ALAS2|AHSP|CA1)$|^(DDX3Y|EIF1AY|UTY|USP9Y|RPS4Y1|KDM5D|NLGN4Y|TMSB4Y|ZFY|TBL1Y|TXLNGY|PRKY|RPS4Y2|XIST|TSIX)$"
TOP_N_GENES         <- get_env_integer("TOP_N_GENES", 50)
MIN_DONORS          <- get_env_integer("MIN_DONORS", 8)
MIN_CELLS           <- get_env_integer("MIN_CELLS", 100)
MIN_CELLS_CLASS     <- get_env_integer("MIN_CELLS_CLASS", 10)
MIN_CELLS_PER_DONOR <- get_env_integer("MIN_CELLS_PER_DONOR", 0)

# --- Common cell-type set ---------------------------------------------------
# The 25 cell types analysed across all paper figures and tables. This is the
# intersection of: SLE-DE pass, SjS-DE pass, singleDeep CV-stable
# (MIN_DONORS=8). The binding constraint is SjS (only 14 SjS donors), which
# excludes 4 rare populations: Low-density basophils, Low-density neutrophils,
# Plasmablasts, Plasmacytoid dendritic cells.
#
# Names are the canonical atlas labels with slashes normalised to spaces
# (matches DE/GSEA/signatures/manuscript figure conventions). Use
# is_common_celltype() to test membership for raw atlas labels (which may
# carry slashes, e.g. "Th1/Th17 cells").
COMMON_CELL_TYPES <- c(
  "Central memory CD8 T cells",
  "Classical monocytes",
  "Effector memory CD8 T cells",
  "Exhausted B cells",
  "Follicular helper T cells",
  "Intermediate monocytes",
  "MAIT cells",
  "Myeloid dendritic cells",
  "Naive B cells",
  "Naive CD4 T cells",
  "Naive CD8 T cells",
  "Natural killer cells",
  "Non classical monocytes",
  "Non switched memory B cells",
  "Non Vd2 gd T cells",
  "Progenitor cells",
  "Switched memory B cells",
  "T regulatory cells",
  "Terminal effector CD4 T cells",
  "Terminal effector CD8 T cells",
  "Th1 cells",
  "Th1 Th17 cells",
  "Th17 cells",
  "Th2 cells",
  "Vd2 gd T cells"
)

is_common_celltype <- function(cts) {
  # Normalise both "/" and "-" to spaces so atlas raw labels
  # ("Th1/Th17 cells", "Non-Vd2 gd T cells") match the space-form
  # COMMON_CELL_TYPES list.
  gsub("[-/]", " ", cts) %in% COMMON_CELL_TYPES
}

# --- Doublet removal: design-appropriate, not procedurally uniform ---
# Every cohort must end in the same state, doublets removed once, which is not the
# same as every cohort receiving the same command. Accessions listed here arrive
# already doublet-filtered by the original authors and MUST NOT be passed through
# scDblFinder a second time.
#
# CXG_436154da (Perez 2022) pooled a median of 16 donors per 10x lane (88
# library_uuid values over 261 donors) and demultiplexed with demuxlet. The
# deposited CELLxGENE object retains only confidently assigned singlets: donor_id
# has 261 clean categories with no doublet, ambiguous or unassigned level. demuxlet
# has therefore already removed cross-donor doublets, which in a 16-donor pool are
# both the most frequent and the only reliably detectable kind.
#
# scDblFinder cannot add to that. It detects doublets by simulating them within a
# grouping variable, and this object carries no technical grouping column that
# survives ingest (library_uuid is not propagated; sample_id is constant), so it can
# only ever be grouped by the biological donor. That splits each pooled lane across
# ~16 groups and makes exactly the doublets demuxlet already removed invisible to
# the simulation, leaving it to remove mostly singlets. Two attempts bear this out:
# ungrouped it removed 39.8% of cells, grouped by donor_id 5.7%, the latter below
# every non-multiplexed cohort here (7.2-8.5%) despite the pooled design. The
# populations recovered between those two runs were strongly non-uniform across
# cell types and displaced the disease-level differential expression ranking.
#
# Do not re-add scDblFinder for these accessions. If a technical capture column is
# ever plumbed through ingest, that still would not justify it: the cross-donor
# doublets are already gone, and a second pass can only subtract singlets.
DOUBLETS_PREFILTERED <- c("CXG_436154da-bcf1-4130-9c8b-120ff9a888f2")

# --- Integration ---
DISEASES_EXCLUDE   <- get_env("DISEASES_EXCLUDE", "RA")
CONDITIONS_EXCLUDE <- get_env("CONDITIONS_EXCLUDE", "sicca")

# --- Figure parameters ---
FIG_WIDTH  <- get_env_numeric("FIG_WIDTH", 12)
FIG_HEIGHT <- get_env_numeric("FIG_HEIGHT", 10)
DPI        <- get_env_integer("DPI", 450)

# --- Color palettes ---
CONDITION_COLORS <- c(
  "Healthy" = "#009E73",
  "SLE"     = "#D55E00",
  "SjS"     = "#0072B2"
)

CONDITION_MAP <- c(
  "healthy" = "Healthy",
  "sle"     = "SLE",
  "sjs"     = "SjS"
)

CONDITION_ORDER <- c("Healthy", "SLE", "SjS")

# --- Lineage grouping (cell types -> 9 lineages) ---
# Includes both hyphenated and non-hyphenated cell-type variants so lookups match
# whatever the upstream annotation pipeline emitted. Plasmablasts, NK cells and
# DCs are separate lineages (used by 11_manuscript_figures.R).
LINEAGE_MAP <- c(
  "Naive CD4 T cells"             = "CD4 T",
  "Th1 cells"                     = "CD4 T",
  "Th2 cells"                     = "CD4 T",
  "Th17 cells"                    = "CD4 T",
  "Th1/Th17 cells"                = "CD4 T",
  "Th1 Th17 cells"                = "CD4 T",
  "Follicular helper T cells"     = "CD4 T",
  "Terminal effector CD4 T cells" = "CD4 T",
  "T regulatory cells"            = "CD4 T",
  "Naive CD8 T cells"             = "CD8 T",
  "Central memory CD8 T cells"    = "CD8 T",
  "Effector memory CD8 T cells"   = "CD8 T",
  "Terminal effector CD8 T cells" = "CD8 T",
  "MAIT cells"                    = "Unconventional T",
  "Vd2 gd T cells"                = "Unconventional T",
  "Non-Vd2 gd T cells"            = "Unconventional T",
  "Non Vd2 gd T cells"            = "Unconventional T",
  "Naive B cells"                 = "B cells",
  "Switched memory B cells"       = "B cells",
  "Non-switched memory B cells"   = "B cells",
  "Non switched memory B cells"   = "B cells",
  "Exhausted B cells"             = "B cells",
  "Plasmablasts"                  = "Plasmablasts",
  "Classical monocytes"           = "Monocytes",
  "Intermediate monocytes"        = "Monocytes",
  "Non-classical monocytes"       = "Monocytes",
  "Non classical monocytes"       = "Monocytes",
  "Myeloid dendritic cells"       = "DCs",
  "Plasmacytoid dendritic cells"  = "DCs",
  "Natural killer cells"          = "NK cells",
  "Low-density neutrophils"       = "Other",
  "Low density neutrophils"       = "Other",
  "Low-density basophils"         = "Other",
  "Low density basophils"         = "Other",
  "Progenitor cells"              = "Other"
)

LINEAGE_ORDER <- c("CD4 T", "CD8 T", "Unconventional T", "B cells",
                   "Plasmablasts", "NK cells", "Monocytes", "DCs", "Other")

LINEAGE_COLORS <- c(
  "CD4 T"            = "#E41A1C",
  "CD8 T"            = "#FF7F00",
  "Unconventional T" = "#FDBF6F",
  "B cells"          = "#377EB8",
  "Plasmablasts"     = "#984EA3",
  "NK cells"         = "#4DAF4A",
  "Monocytes"        = "#A65628",
  "DCs"              = "#F781BF",
  "Other"            = "#999999"
)


# --- Diseases analyzed (order matches CONDITION_ORDER: SLE before SjS) ---
DISEASES_ACTIVE <- c("SLE", "SjS")

# --- Neural network classification (script 11) ---
# N_FOLDS, EPOCHS, BATCH_SIZE, MIN_CELLS_PER_TYPE
# are set via env vars and read directly by the Python script

# --- Parallelization ---
# Auto-detect cores: SLURM_CPUS_PER_TASK > N_CORES env var > 1
N_CORES <- get_env_integer("N_CORES",
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
)
cat("Cores available:", N_CORES, "\n")

# Set up future plan for Seurat parallel operations (ScaleData, etc.)
if (N_CORES > 1 && requireNamespace("future", quietly = TRUE)) {
  future::plan("multicore", workers = N_CORES)
  # Large atlas integration needs more than 8GB for globals
  options(future.globals.maxSize = 100 * 1024^3)  # 100 GB
}

# Set up BiocParallel for Bioconductor packages (DESeq2, fgsea, etc.)
if (N_CORES > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::MulticoreParam(workers = N_CORES))
}

# --- Output directories ---
OBJ_DIR <- file.path(ATLAS_ROOT, "results", "objects")

get_env_logical <- function(key, default) {
  val <- tolower(get_env(key, as.character(default)))
  val %in% c("true", "1", "yes")
}

RUN_SUPPLEMENTARY <- get_env_logical("RUN_SUPPLEMENTARY", FALSE)

ATLAS_PATH <- file.path(OBJ_DIR, "atlas_integrated.rds")

if (RUN_SUPPLEMENTARY) {
  cat("*** SUPPLEMENTARY ANALYSES ENABLED ***\n")
}

# --- Reference data ---
SINGLER_REF <- file.path(ATLAS_ROOT, "data", "references", "singler_monaco.rds")
