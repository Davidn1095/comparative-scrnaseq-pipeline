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

# Genes excluded from DE / GSEA / classifier feature space, by exact rule
# (is_excluded_gene() below):
# - mitochondrial genes: the MT- prefix (technical / stress correlate)
# - ribosomal protein genes: the HGNC "L ribosomal proteins" (729) and "S ribosomal
#   proteins" (728) groups as HGNC lists them (hgnc_complete_set.txt, 2026-09-28).
#   Genes that merely share the RPL/RPS prefix, such as the RPS6KA kinases, are kept.
# - haemoglobin genes: the ten globin chains listed below (red-cell contamination).
#   HBEGF is not a haemoglobin gene and is kept.
# - Y-chromosome genes: every atlas feature whose gene is annotated to chromosome Y,
#   by Ensembl gene ID where the dataset deposits one and by gene name otherwise
#   (GSE253568), in Ensembl 116, GENCODE 32, 28 and 24 or, for the GRCh37 Perez
#   annotation, Ensembl 75. Pseudoautosomal genes are annotated to X and are kept.
# - XIST, TSIX. Sex-chromosome genes are confounded by the cohort sex imbalance: CXG
#   SLE has 9.3% male donors vs CXG HC 2.0%, which drives spurious DDX3Y/EIF1AY
#   upregulation in SLE that is not biology.
# - erythroid markers SLC4A1, ALAS2, AHSP, CA1: the red-cell contamination the
#   haemoglobin genes miss (SLC4A1 was the largest SLE fold change).
# - dissociation-induced genes: the van den Brink et al. 2017 list (Nat Methods 14:935,
#   doi:10.1038/nmeth.4437, Supplementary Table 5, 140 mouse genes), its 154 human orthologs
#   from the MGI homology report and the older atlas names of four of them, plus HSPA6, a
#   stress-inducible human HSP70 with no mouse ortholog. Provenance and mapping:
#   scripts/reference/dissociation_genes_vandenbrink2017.tsv. The healthy controls of
#   GSE157278 carry about twice the dissociation-gene reads of every other group (5.2% of
#   pseudobulk reads against 2.0-3.1%), which turned their stress response into apparent
#   SjD decreases.
EXCLUDED_RIBOSOMAL_GENES <- c(
  "FAU", "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14",
  "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23",
  "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3",
  "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL",
  "RPL37", "RPL37A", "RPL38", "RPL39", "RPL39L", "RPL3L", "RPL4", "RPL41", "RPL5",
  "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPLP0", "RPLP1", "RPLP2", "RPS10",
  "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18",
  "RPS19", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27",
  "RPS27A", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2",
  "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPSA", "UBA52"
)
EXCLUDED_HAEMOGLOBIN_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ",
                                "HBM", "HBQ1")
EXCLUDED_Y_GENES <- c(
  "AC006040.1", "AC006157.1", "AC006328.4", "AC006386.1", "AC007244.1", "AC007359.1",
  "AC007359.6", "AC007876.1", "AC008175.1", "AC009491.1", "AC009491.2", "AC009494.2",
  "AC009977.1", "AC010084.1", "AC010086.3", "AC010722.1", "AC010723.1", "AC010737.1",
  "AC010889.1", "AC010889.2", "AC010891.1", "AC010891.2", "AC011297.1", "AC011751.1",
  "AC012005.3", "AC012005.4", "AC012078.2", "AC022486.1", "AC024236.1", "AC064829.1",
  "AC244213.1", "AMELY", "BPY2", "BPY2B", "BPY2C", "CDY1", "CDY1B", "CDY2A", "CDY2B",
  "CSPG4P1Y", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DDX3Y", "EIF1AY", "ENSG00000223517",
  "ENSG00000228379", "ENSG00000229308", "ENSG00000235059", "ENSG00000236951",
  "ENSG00000251510", "ENSG00000254488", "ENSG00000260197", "FAM197Y1", "FAM197Y2",
  "FAM197Y3", "FAM197Y5", "FAM197Y5P", "FAM197Y6", "FAM197Y7", "FAM197Y8", "FAM224A",
  "FAM224B", "FAM41AY1", "FAM41AY2", "HSFY1", "HSFY2", "KDM5D", "LINC00266-4P",
  "LINC00278", "LINC00279", "LINC00280", "NLGN4Y", "NLGN4Y-AS1", "PCDH11Y", "PRKY",
  "PRORY", "PRY", "PRY2", "PRYP3", "RBMY1A1", "RBMY1B", "RBMY1D", "RBMY1E", "RBMY1F",
  "RBMY1J", "RP11-122L9.1", "RP11-414C23.1", "RP11-424G14.1", "RP11-65G9.1", "RPS4Y1",
  "RPS4Y2", "SEPTIN14P23", "SLC9B1P1", "SRY", "TBL1Y", "TGIF2LY", "TMSB4Y", "TSPY1",
  "TSPY10", "TSPY2", "TSPY3", "TSPY4", "TSPY8", "TSPY9", "TSPY9P", "TTTY1", "TTTY10",
  "TTTY11", "TTTY12", "TTTY13", "TTTY14", "TTTY15", "TTTY16", "TTTY17A", "TTTY17B",
  "TTTY17C", "TTTY18", "TTTY19", "TTTY1B", "TTTY2", "TTTY20", "TTTY21", "TTTY21B",
  "TTTY22", "TTTY23", "TTTY23B", "TTTY2B", "TTTY3", "TTTY3B", "TTTY4", "TTTY4B",
  "TTTY4C", "TTTY5", "TTTY6", "TTTY6B", "TTTY7", "TTTY7B", "TTTY8", "TTTY8B", "TTTY9B",
  "USP9Y", "UTY", "VCY", "VCY1B", "ZFY", "ZFY-AS1"
)
EXCLUDED_DISSOCIATION_GENES <- c(
  "ACTG1", "ANKRD1", "ARID5A", "ATF3", "ATF4", "BAG3", "BHLHE40", "BRD2", "BTG1",
  "BTG2", "CCN1", "CCNL1", "CEBPB", "CEBPD", "CEBPG", "CSRNP1", "CXCL1", "CXCL2",
  "CXCL3", "CYR61", "DCN", "DDX3X", "DDX5", "DES", "DNAJA1", "DNAJB1", "DNAJB4",
  "DUSP1", "DUSP8", "EGR1", "EGR2", "EIF1", "EIF5", "ERF", "ERFE", "ERRFI1", "FAM132B",
  "FOS", "FOSB", "FOSL2", "GADD45A", "GADD45G", "GCC1", "GEM", "H3-3B", "H3-5", "H3F3B",
  "H3F3C", "HIPK3", "HSP90AA1", "HSP90AA4P", "HSP90AB1", "HSP90AB3P", "HSPA1A",
  "HSPA1B", "HSPA5", "HSPA6", "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IDI1", "IER2",
  "IER3", "IER5", "IFRD1", "IL6", "IRF1", "IRF8", "ITPKC", "JUN", "JUNB", "JUND",
  "KCNE4", "KLF2", "KLF4", "KLF6", "KLF9", "LITAF", "LMNA", "MAFF", "MAFK", "MCL1",
  "MIDN", "MIR22HG", "MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H", "MT1L", "MT1M",
  "MT1X", "MT2A", "MYADM", "MYC", "MYD88", "NCKAP5L", "NCOA7", "NFKBIA", "NFKBIZ",
  "NOCT", "NOP58", "NPPC", "NR4A1", "ODC1", "OSGIN1", "OXNAD1", "PCF11", "PDE4B",
  "PER1", "PHLDA1", "PNP", "PNRC1", "PPP1CC", "PPP1R15A", "PXDC1", "RAP1B", "RAP1BL",
  "RASSF1", "RHOB", "RHOH", "RIPK1", "SAT1", "SBNO2", "SDC4", "SERPINE1", "SKIL",
  "SLC10A6", "SLC38A2", "SLC41A1", "SOCS3", "SQSTM1", "SRF", "SRSF5", "SRSF7", "STAT3",
  "TAGLN2", "TIPARP", "TNFAIP3", "TNFAIP6", "TPM3", "TPPP3", "TRA2A", "TRA2B", "TRIB1",
  "TUBB4B", "TUBB6", "UBC", "USP2", "WAC", "ZC3H12A", "ZFAND5", "ZFP36", "ZFP36L1",
  "ZFP36L2", "ZYX"
)
EXCLUDED_GENES <- unique(c(EXCLUDED_RIBOSOMAL_GENES, EXCLUDED_HAEMOGLOBIN_GENES,
                          EXCLUDED_Y_GENES, "XIST", "TSIX",
                          "SLC4A1", "ALAS2", "AHSP", "CA1",
                          EXCLUDED_DISSOCIATION_GENES))
is_excluded_gene <- function(genes) {
  genes <- as.character(genes)
  startsWith(genes, "MT-") | genes %in% EXCLUDED_GENES | genes %in% excluded_genes_hgnc()
}

# --- Gene symbol harmonisation (HGNC complete set, 2026-09-28) ---
# The four datasets name some genes differently (RIGI/DDX58, ATP5F1E/ATP5E, SELENOF/SEP15),
# so their gene names are mapped onto one reference symbol set before any cross-dataset step
# (harmonise_genes() in 00_utils.R, applied at the start of 03 and 04a). The reference is the
# HGNC complete set downloaded on 2026-09-28 (scripts/reference/hgnc_complete_set_2026-09-28.README).
# Symbols map by approved symbol, then by a previous symbol and then by an alias, each only when
# it points to exactly one approved symbol; anything else is left as it is. The same mapping is
# applied to the Hallmark gene sets, the exclusion lists above and the Toro-Dominguez signature.
HGNC_FILE <- file.path(ATLAS_ROOT, "scripts", "reference", "hgnc_complete_set_2026-09-28.txt.gz")
.hgnc <- new.env(parent = emptyenv())
hgnc_table <- function() {
  if (is.null(.hgnc$approved)) {
    if (!file.exists(HGNC_FILE)) stop("HGNC reference not found: ", HGNC_FILE, call. = FALSE)
    h <- utils::read.delim(gzfile(HGNC_FILE), colClasses = "character", quote = "",
                           na.strings = "", comment.char = "")
    h <- h[h$status == "Approved", c("hgnc_id", "symbol", "prev_symbol", "alias_symbol",
                                      "ensembl_gene_id")]
    # multi-valued fields are quoted in the HGNC file ("H1FV|H1F0")
    for (col in c("prev_symbol", "alias_symbol")) h[[col]] <- gsub("\"", "", h[[col]], fixed = TRUE)
    # previous symbols and aliases that name exactly one approved entry and are not themselves
    # an approved symbol; ambiguous ones are left unmapped
    unique_keys <- function(col) {
      v <- strsplit(ifelse(is.na(h[[col]]), "", h[[col]]), "|", fixed = TRUE)
      d <- unique(data.frame(key = unlist(v), symbol = rep(h$symbol, lengths(v)),
                             stringsAsFactors = FALSE))
      d <- d[nzchar(d$key) & !d$key %in% h$symbol, ]
      n <- table(d$key)
      d <- d[d$key %in% names(n)[n == 1], ]
      stats::setNames(d$symbol, d$key)
    }
    ens <- h[!is.na(h$ensembl_gene_id), ]
    dup_ens <- ens$ensembl_gene_id[duplicated(ens$ensembl_gene_id)]
    ens <- ens[!ens$ensembl_gene_id %in% dup_ens, ]
    .hgnc$approved <- h$symbol
    .hgnc$ensembl  <- stats::setNames(ens$symbol, ens$ensembl_gene_id)
    .hgnc$previous <- unique_keys("prev_symbol")
    .hgnc$alias    <- unique_keys("alias_symbol")
  }
  .hgnc
}
# One row per input symbol: the approved HGNC symbol it maps to (NA if none) and how.
hgnc_map_symbols <- function(sym) {
  H <- hgnc_table()
  sym <- as.character(sym)
  target <- rep(NA_character_, length(sym)); how <- rep("unmapped", length(sym))
  i <- sym %in% H$approved;                       target[i] <- sym[i];               how[i] <- "approved"
  j <- is.na(target) & sym %in% names(H$previous); target[j] <- H$previous[sym[j]]; how[j] <- "previous"
  k <- is.na(target) & sym %in% names(H$alias);    target[k] <- H$alias[sym[k]];    how[k] <- "alias"
  data.frame(original = sym, target = unname(target), how = how, stringsAsFactors = FALSE)
}
# The HGNC symbol of each input, or the input itself where it does not map.
harmonise_symbols <- function(sym) {
  m <- hgnc_map_symbols(sym)
  ifelse(is.na(m$target), m$original, m$target)
}
# The exclusion lists above in the HGNC namespace, so that a gene stays excluded under the
# symbol it is harmonised to.
excluded_genes_hgnc <- function() {
  if (is.null(.hgnc$excluded)) .hgnc$excluded <- unique(harmonise_symbols(EXCLUDED_GENES))
  .hgnc$excluded
}
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

# Display names for cell types in figures (script 10) and supplementary tables
# (script 11), matching the manuscript text. Keys are the internal names with
# "-", "_" and "/" read as spaces; names not listed are shown as they are, and
# any other string (a gene, a pathway) is returned untouched.
CELLTYPE_DISPLAY <- c(
  "Non classical monocytes"     = "Non-classical monocytes",
  "Non switched memory B cells" = "Non-switched memory B cells",
  "Non Vd2 gd T cells"          = "Non-V\u03b42 \u03b3\u03b4 T cells",
  "Vd2 gd T cells"              = "V\u03b42 \u03b3\u03b4 T cells",
  "Th1 Th17 cells"              = "Th1/Th17 cells"
)
display_celltype <- function(x) {
  x <- as.character(x)
  key <- gsub("[-_/]", " ", x)
  hit <- !is.na(key) & key %in% names(CELLTYPE_DISPLAY)
  x[hit] <- unname(CELLTYPE_DISPLAY[key[hit]])
  x
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

# --- One sample per donor ---
# CXG_436154da (Perez 2022) holds 274 samples from 261 donors: 11 SLE patients were
# sampled two or three times, at a flare, after its treatment and, for three of them,
# at a stable visit. Every analysis here aggregates by donor, so those visits would be
# pooled into one profile. Each patient keeps the first visit only
# (keep_first_visit() in 00_utils.R), applied after QC and annotation, in
# 03_pseudobulk.R and 04a_integrate.R: 13 samples and 32,637 cells are removed.
FIRST_VISIT_ONLY <- c("CXG_436154da-bcf1-4130-9c8b-120ff9a888f2")

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
