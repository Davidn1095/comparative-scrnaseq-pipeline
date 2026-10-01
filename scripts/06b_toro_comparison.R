#!/usr/bin/env Rscript
################################################################################
# 06b_toro_comparison.R
#
# Recovery of the shared SLE/RA/SjS bulk-PBMC signature of Toro-Dominguez et al. 2014
# (Arthritis Res Ther 16:489, doi:10.1186/s13075-014-0489-x) in the atlas DEGs.
#
# Inputs
#   - scripts/reference/toro2014/ (Additional file 1 of the paper; provenance in its README):
#       toro2014_signature.csv   the signature genes with direction; asterisk genes dropped,
#                                leaving 371 (188 up, 183 down)
#       toro2014_go_terms.tsv    the 115 GO biological-process terms of the paper; the 68
#                                with P < 0.01 (TORO_GO_P) are used, up and down
#       toro2014_gain_genes.csv  the 132 genes found only by the meta-analysis
#     Every symbol is mapped onto the HGNC reference with harmonise_symbols() (00_config.R).
#   - Atlas SLE and SjS DEGs: results/plots/differential_expression/de_{sle,sjs}/deseq2_<CT>.tsv (06)
#   - Atlas GSEA: results/plots/differential_expression/gsea_significant.tsv (07)
#
# A signature gene is recovered in a disease when it is a DEG (padj < 0.05, |log2FC| >= 1) in
# the signature's direction in at least one of the 25 cell types; in both diseases when that
# holds for each disease, not necessarily in the same cell type.
#
# Outputs (results/manuscript1_figures/toro_comparison/, or TORO_OUT_DIR)
#   - toro_go_pathway_recovery.tsv   per GO term: genes, testable genes, recovered in SLE,
#                                    SjS and both (Supplementary Table 8, script 11)
#   - toro_shared_hits.tsv           per signature gene: DEG status by disease and direction
#   - toro_shared_localisation.tsv   cell types of the genes recovered in both diseases
#   - overlap_summary.tsv            gene-level counts, gain genes included
#   - toro_hallmark_concordance.tsv, atlas_shared_GOBP_hypergeometric.tsv,
#     toro_term_overlap_with_atlas_hypergeometric.tsv   supporting tables, not used in main.tex
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

ROOT    <- Sys.getenv("ATLAS_ROOT", getwd())   # repository root
# harmonise_symbols(): the 2014 symbols onto the HGNC reference the atlas genes use
source(file.path(ROOT, "scripts", "00_config.R"))
DE_DIR  <- file.path(ROOT, "results/plots/differential_expression")
OUT_DIR <- Sys.getenv("TORO_OUT_DIR", file.path(ROOT, "results/manuscript1_figures/toro_comparison"))
TORO_REF <- file.path(ROOT, "scripts", "reference", "toro2014")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ  <- 0.05
LFC   <- 1

# -----------------------------------------------------------------------------
# Toro-Dominguez 2014 "gain" genes (Additional file 1, sheet "Gain and Loss genes", n=132)
# Subset of the shared signature detectable only via meta-analysis (weak signal)
# -----------------------------------------------------------------------------
gain_csv <- file.path(TORO_REF, "toro2014_gain_genes.csv")
if (file.exists(gain_csv)) {
  gain_lines <- readLines(gain_csv)
  start <- grep("^BRCA1;", gain_lines)[1]
  gain_block <- gain_lines[start:length(gain_lines)]
  gain_block <- gain_block[!grepl("^Loss genes", gain_block) &
                           !grepl("^Loss genes that", gain_block)]
  # Loss block would have started by now, but file only contains gain section
  toro_gain <- unlist(strsplit(gain_block, ";"))
  toro_gain <- trimws(toro_gain)
  toro_gain <- toro_gain[nzchar(toro_gain)]
  toro_gain <- unique(harmonise_symbols(toro_gain))
  cat(sprintf("Toro gain genes (Additional file 1): %d\n", length(toro_gain)))
} else {
  toro_gain <- character(0)
}

# -----------------------------------------------------------------------------
# Toro-Dominguez 2014 shared signature (Additional file 1, sheet "GENE EXPRESSION SIGNATURE")
# -----------------------------------------------------------------------------
sig_csv <- file.path(TORO_REF, "toro2014_signature.csv")
if (file.exists(sig_csv)) {
  raw <- readLines(sig_csv)
  # Skip header rows; find "EntrezID;Name;..." line as start
  hdr_idx <- grep("^EntrezID;", raw)[1]
  data_lines <- raw[(hdr_idx + 1):length(raw)]
  # Each data line has semicolon-delimited: up_entrez;up_name;up_es;up_pval;;;down_entrez;down_name;down_es;down_pval
  parts <- strsplit(data_lines, ";", fixed = TRUE)
  extract <- function(p, side) {
    cols <- if (side == "up") 1:4 else 7:10
    vals <- lapply(p, function(r) {
      if (length(r) >= max(cols)) r[cols] else rep(NA, 4)
    })
    do.call(rbind, vals)
  }
  up_mat   <- extract(parts, "up")
  down_mat <- extract(parts, "down")
  parse_block <- function(m, direction) {
    df <- data.frame(entrez = m[,1], gene = m[,2],
                     es = as.numeric(sub(",", ".", m[,3])),
                     pval = as.numeric(sub(",", ".", m[,4])),
                     stringsAsFactors = FALSE)
    df <- df[nzchar(trimws(df$gene)) & !is.na(df$gene), ]
    # asterisk = not DE in all 3 diseases — strip marker, flag, then remove (paper does exclude them)
    df$incomplete <- grepl("\\*$", df$gene)
    df$gene <- harmonise_symbols(sub("\\*$", "", df$gene))
    df$direction <- direction
    df
  }
  toro_up_df   <- parse_block(up_mat,   "up")
  toro_down_df <- parse_block(down_mat, "down")
  # Paper's "shared signature" excludes asterisk genes (not DE in all diseases)
  toro_up   <- unique(toro_up_df$gene[!toro_up_df$incomplete])
  toro_down <- unique(toro_down_df$gene[!toro_down_df$incomplete])
  toro_up_all   <- unique(toro_up_df$gene)
  toro_down_all <- unique(toro_down_df$gene)
  cat(sprintf("Toro shared signature (Additional file 1): %d up, %d down after asterisk filter\n",
              length(toro_up), length(toro_down)))
  cat(sprintf("  (including asterisk-flagged: %d up, %d down)\n",
              length(toro_up_all), length(toro_down_all)))
} else {
  stop("Toro-Dominguez signature not found: ", sig_csv)
}

toro_all  <- sort(unique(c(toro_up, toro_down)))
cat(sprintf("Toro shared signature (final, asterisk-excluded): %d up + %d down = %d\n",
            length(toro_up), length(toro_down), length(toro_all)))

# -----------------------------------------------------------------------------
# Atlas DEGs per cell type, per disease
# -----------------------------------------------------------------------------
load_de <- function(disease) {
  dir <- file.path(DE_DIR, paste0("de_", disease))
  files <- list.files(dir, pattern = "^deseq2_.*\\.tsv$", full.names = TRUE)
  files <- files[!grepl("_annotated", basename(files))]
  lapply(files, function(f) {
    ct <- sub("^deseq2_", "", sub("\\.tsv$", "", basename(f)))
    df <- suppressMessages(read_tsv(f, show_col_types = FALSE))
    df$cell_type <- ct
    df$disease   <- disease
    df
  }) |> bind_rows()
}

sle <- load_de("sle") |> filter(!is.na(padj))
sjs <- load_de("sjs") |> filter(!is.na(padj))

cat(sprintf("SLE: %d cell types, %d total tests | SjS: %d cell types, %d total tests\n",
            length(unique(sle$cell_type)), nrow(sle),
            length(unique(sjs$cell_type)), nrow(sjs)))

sig <- function(df, direction = c("up","down","any")) {
  direction <- match.arg(direction)
  x <- df |> filter(padj < PADJ, abs(log2FoldChange) >= LFC)
  if (direction == "up")   x <- x |> filter(log2FoldChange >=  LFC)
  if (direction == "down") x <- x |> filter(log2FoldChange <= -LFC)
  x
}

# Genes DE in BOTH diseases (shared) in at least one cell type, same direction
shared_up   <- intersect(unique(sig(sle,"up")$gene),   unique(sig(sjs,"up")$gene))
shared_down <- intersect(unique(sig(sle,"down")$gene), unique(sig(sjs,"down")$gene))
shared_any  <- unique(c(shared_up, shared_down))

cat(sprintf("\nAtlas SLE n SjS shared DEGs (padj<%.2f, |log2FC|>=%g, same direction, any cell type):\n",
            PADJ, LFC))
cat(sprintf("  up:   %d\n  down: %d\n  total unique: %d\n",
            length(shared_up), length(shared_down), length(shared_any)))

# -----------------------------------------------------------------------------
# Overlap with Toro-Dominguez
# -----------------------------------------------------------------------------
ov_up   <- intersect(toro_up,   shared_up)
ov_down <- intersect(toro_down, shared_down)
ov_any  <- intersect(toro_all,  shared_any)

# Also: Toro gene DE in SLE alone, SjS alone, both, neither
toro_hits <- tibble(
  gene      = toro_all,
  toro_dir  = ifelse(toro_all %in% toro_up, "up",
                     ifelse(toro_all %in% toro_down, "down", NA)),
  in_sle_up   = toro_all %in% unique(sig(sle,"up")$gene),
  in_sle_down = toro_all %in% unique(sig(sle,"down")$gene),
  in_sjs_up   = toro_all %in% unique(sig(sjs,"up")$gene),
  in_sjs_down = toro_all %in% unique(sig(sjs,"down")$gene)
) |>
  mutate(
    sle_concordant = (toro_dir == "up"   & in_sle_up)   | (toro_dir == "down" & in_sle_down),
    sjs_concordant = (toro_dir == "up"   & in_sjs_up)   | (toro_dir == "down" & in_sjs_down),
    both_concordant = sle_concordant & sjs_concordant
  )

write_tsv(toro_hits, file.path(OUT_DIR, "toro_shared_hits.tsv"))

# Cell-type localisation for concordant genes: which cell types drive them
ct_localisation <- function(gene, disease_df, dir) {
  sub <- disease_df |> filter(gene == !!gene, padj < PADJ,
                              (dir == "up"   & log2FoldChange >=  LFC) |
                              (dir == "down" & log2FoldChange <= -LFC))
  if (nrow(sub) == 0) return(NA_character_)
  paste(sort(unique(sub$cell_type)), collapse = ";")
}

loc_rows <- toro_hits |>
  filter(both_concordant) |>
  rowwise() |>
  mutate(
    sle_cell_types = ct_localisation(gene, sle, toro_dir),
    sjs_cell_types = ct_localisation(gene, sjs, toro_dir)
  ) |>
  ungroup()

write_tsv(loc_rows, file.path(OUT_DIR, "toro_shared_localisation.tsv"))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
# Additional: gain-gene recovery
# (Toro gain genes that are DE in SLE, SjS, or both in our atlas, any direction)
if (length(toro_gain) > 0) {
  any_de_sle <- unique(c(sig(sle,"up")$gene, sig(sle,"down")$gene))
  any_de_sjs <- unique(c(sig(sjs,"up")$gene, sig(sjs,"down")$gene))
  gain_in_sle   <- intersect(toro_gain, any_de_sle)
  gain_in_sjs   <- intersect(toro_gain, any_de_sjs)
  gain_in_both  <- intersect(gain_in_sle, gain_in_sjs)
  gain_shared_concordant <- intersect(toro_gain, shared_any)
} else {
  gain_in_sle <- gain_in_sjs <- gain_in_both <- gain_shared_concordant <- character(0)
}

summary_df <- tibble(
  metric = c(
    "Toro up-regulated (transcribed)",
    "Toro down-regulated (transcribed)",
    "Atlas SLE-SjS shared up (>=1 cell type)",
    "Atlas SLE-SjS shared down (>=1 cell type)",
    "Overlap: Toro up ∩ Atlas shared up",
    "Overlap: Toro down ∩ Atlas shared down",
    "Overlap (any direction, sign-concordant)",
    "Toro genes concordant in SLE only",
    "Toro genes concordant in SjS only",
    "Toro genes concordant in BOTH diseases",
    "Toro GAIN genes (Additional file 1)",
    "  ...DE in SLE (any cell type)",
    "  ...DE in SjS (any cell type)",
    "  ...DE in both diseases",
    "  ...shared and sign-concordant"
  ),
  value = c(
    length(toro_up),
    length(toro_down),
    length(shared_up),
    length(shared_down),
    length(ov_up),
    length(ov_down),
    length(ov_up) + length(ov_down),
    sum(toro_hits$sle_concordant & !toro_hits$sjs_concordant, na.rm = TRUE),
    sum(toro_hits$sjs_concordant & !toro_hits$sle_concordant, na.rm = TRUE),
    sum(toro_hits$both_concordant, na.rm = TRUE),
    length(toro_gain),
    length(gain_in_sle),
    length(gain_in_sjs),
    length(gain_in_both),
    length(gain_shared_concordant)
  )
)

write_tsv(summary_df, file.path(OUT_DIR, "overlap_summary.tsv"))
print(summary_df, n = Inf)

# -----------------------------------------------------------------------------
# Hallmark pathway concordance
# -----------------------------------------------------------------------------
gsea <- suppressMessages(read_tsv(file.path(DE_DIR, "gsea_significant.tsv"),
                                  show_col_types = FALSE))
# clean disease labels + strip duplicated "_annotated" suffixes
gsea$disease   <- toupper(sub("^DE_", "", gsea$disease))
gsea$disease   <- ifelse(gsea$disease == "SJS", "SjS", gsea$disease)
gsea$cell_type <- trimws(gsub("( annotated)+$", "", gsea$cell_type))
gsea <- gsea |> distinct(pathway, disease, cell_type, .keep_all = TRUE) |>
  filter(disease %in% c("SLE", "SjS"))

# Toro-reported enriched themes mapped to MSigDB Hallmark tokens
toro_pathways <- list(
  "Type I IFN"         = c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  "Cytokine signaling" = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                           "HALLMARK_INFLAMMATORY_RESPONSE"),
  "Response to virus"  = c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  "Mitotic cell cycle" = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE"),
  "Apoptosis"          = c("HALLMARK_APOPTOSIS"),
  "Translation (down)" = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                           "HALLMARK_MTORC1_SIGNALING")
)

conc <- lapply(names(toro_pathways), function(theme) {
  tokens <- toro_pathways[[theme]]
  hits <- gsea |> filter(pathway %in% tokens)
  tibble(
    theme     = theme,
    token     = tokens,
    n_sle_ct  = sapply(tokens, function(t) length(unique(hits$cell_type[hits$pathway == t & hits$disease == "SLE"]))),
    n_sjs_ct  = sapply(tokens, function(t) length(unique(hits$cell_type[hits$pathway == t & hits$disease == "SjS"]))),
    median_nes_sle = sapply(tokens, function(t) {
      v <- hits$NES[hits$pathway == t & hits$disease == "SLE"]
      if (length(v)) median(v) else NA
    }),
    median_nes_sjs = sapply(tokens, function(t) {
      v <- hits$NES[hits$pathway == t & hits$disease == "SjS"]
      if (length(v)) median(v) else NA
    })
  )
}) |> bind_rows()

write_tsv(conc, file.path(OUT_DIR, "toro_hallmark_concordance.tsv"))
cat("\nHallmark concordance:\n"); print(conc, n = Inf)

# -----------------------------------------------------------------------------
# Toro GO categories (Additional file 1, sheet "GO Analysis") — per-term gene recovery
# Every GO biological-process term the original study reports as enriched, P < 0.01
# (TORO_GO_P), for the over- and the under-expressed genes, with no minimum size.
# A member gene is testable when it has an adjusted P in at least one cell type of both
# diseases; it is recovered in a disease when it is a DEG in the term's direction there.
# -----------------------------------------------------------------------------
TORO_GO_P <- as.numeric(Sys.getenv("TORO_GO_P", "0.01"))
go_all <- read.delim(file.path(TORO_REF, "toro2014_go_terms.tsv"), stringsAsFactors = FALSE,
                     colClasses = c(p_value = "numeric"))
go_df <- go_all[go_all$p_value < TORO_GO_P, ]
cat(sprintf("\nToro GO terms with P < %g: %d of %d (%d up, %d down)\n", TORO_GO_P, nrow(go_df),
            nrow(go_all), sum(go_df$direction == "up"), sum(go_df$direction == "down")))

sle_up   <- unique(sig(sle,"up")$gene);   sle_down <- unique(sig(sle,"down")$gene)
sjs_up   <- unique(sig(sjs,"up")$gene);   sjs_down <- unique(sig(sjs,"down")$gene)
testable <- intersect(unique(sle$gene), unique(sjs$gene))

pw_recov <- do.call(rbind, lapply(seq_len(nrow(go_df)), function(i) {
  r <- go_df[i, ]
  g <- unique(harmonise_symbols(strsplit(r$genes, " ", fixed = TRUE)[[1]]))
  up <- r$direction == "up"
  in_sle <- g[g %in% if (up) sle_up else sle_down]
  in_sjs <- g[g %in% if (up) sjs_up else sjs_down]
  both <- intersect(in_sle, in_sjs)
  data.frame(direction = r$direction, go_id = r$go_id, term = r$term, p_value = r$p_value,
             n_genes = length(g), n_testable = sum(g %in% testable),
             n_recovered_sle = length(in_sle), n_recovered_sjs = length(in_sjs),
             n_recovered_both = length(both), recovered_both = paste(both, collapse = ";"),
             stringsAsFactors = FALSE)
}))
write_tsv(pw_recov, file.path(OUT_DIR, "toro_go_pathway_recovery.tsv"))
cat("\nToro GO term recovery (by fraction recovered in both diseases):\n")
print(as_tibble(pw_recov) |>
        mutate(frac_both = round(n_recovered_both / n_genes, 2)) |>
        arrange(desc(frac_both)) |>
        select(direction, go_id, term, n_genes, n_testable, n_recovered_both, frac_both) |>
        head(20), n = 20)

# -----------------------------------------------------------------------------
# Apples-to-apples: hypergeometric GO BP enrichment on our shared DEGs
# (matches Toro-Dominguez's GeneCodis hypergeometric GO BP analysis)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

run_enrich <- function(genes, universe, label) {
  eg <- bitr(genes,   fromType = "SYMBOL", toType = "ENTREZID",
             OrgDb = org.Hs.eg.db, drop = TRUE)
  bg <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
             OrgDb = org.Hs.eg.db, drop = TRUE)
  res <- enrichGO(gene = eg$ENTREZID, universe = bg$ENTREZID,
                  OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", qvalueCutoff = 0.05,
                  readable = TRUE)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble())
  as.data.frame(res) |> as_tibble() |> mutate(gene_set = label)
}

# Universe: all genes tested in DESeq2 across cell types
universe <- unique(c(sle$gene, sjs$gene))

cat("\nRunning hypergeometric GO BP enrichment on shared DEGs...\n")
enr_up    <- run_enrich(shared_up,   universe, "shared_up")
enr_down  <- run_enrich(shared_down, universe, "shared_down")
enr_all   <- bind_rows(enr_up, enr_down)
write_tsv(enr_all, file.path(OUT_DIR, "atlas_shared_GOBP_hypergeometric.tsv"))

# The selected Toro GO terms (P < TORO_GO_P), by direction
toro_terms_up   <- go_df$go_id[go_df$direction == "up"]
toro_terms_down <- go_df$go_id[go_df$direction == "down"]

shared_term_overlap <- tibble(
  toro_term = c(toro_terms_up, toro_terms_down),
  toro_direction = c(rep("up", length(toro_terms_up)),
                     rep("down", length(toro_terms_down))),
  in_our_up   = c(toro_terms_up, toro_terms_down) %in% enr_up$ID,
  in_our_down = c(toro_terms_up, toro_terms_down) %in% enr_down$ID
)
write_tsv(shared_term_overlap, file.path(OUT_DIR, "toro_term_overlap_with_atlas_hypergeometric.tsv"))

cat("\nApples-to-apples: Toro GO BP terms recovered in our atlas hypergeometric?\n")
print(shared_term_overlap, n = Inf)

cat(sprintf("\nOur shared-UP DEGs: %d significant GO BP terms (q<0.05)\n",
            nrow(enr_up)))
cat(sprintf("Our shared-DOWN DEGs: %d significant GO BP terms (q<0.05)\n",
            nrow(enr_down)))

cat("\nOutputs written to: ", OUT_DIR, "\n")
