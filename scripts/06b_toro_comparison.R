#!/usr/bin/env Rscript
################################################################################
# 06b_toro_comparison.R
#
# Compare the shared SLE+RA+SjS bulk-PBMC gene signature from Toro-Dominguez
# et al. 2014 (Arthritis Res Ther, PMC4295333) against our single-cell atlas
# SLE and SjS DEGs and Hallmark GSEA.
#
# Inputs
#   - Toro-Dominguez gene list: hard-coded from supplementary screenshot
#     (approximate OCR of Additional file 2; top hits are reliable, long tail
#     contains OCR noise — treat overlap counts as lower bounds).
#   - Atlas SLE DEGs: results/plots/differential_expression/de_sle/deseq2_<CT>.tsv
#   - Atlas SjS DEGs: results/plots/differential_expression/de_sjs/deseq2_<CT>.tsv
#   - Atlas GSEA: results/plots/differential_expression/gsea_significant.tsv
#
# Outputs  (results/manuscript1_figures/toro_comparison/)
#   - overlap_summary.tsv
#   - toro_shared_hits.tsv
#   - toro_hallmark_concordance.tsv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

ROOT    <- Sys.getenv("ATLAS_ROOT", getwd())   # repository root
DE_DIR  <- file.path(ROOT, "results/plots/differential_expression")
OUT_DIR <- file.path(ROOT, "results/manuscript1_figures/toro_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ  <- 0.05
LFC   <- 1

# -----------------------------------------------------------------------------
# Toro-Dominguez 2014 "gain" genes (Additional file 3, n=132)
# Subset of the shared signature detectable only via meta-analysis (weak signal)
# -----------------------------------------------------------------------------
gain_csv <- file.path(OUT_DIR, "toro2014_addfile3_gain_loss.csv")
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
  toro_gain <- unique(toro_gain)
  cat(sprintf("Toro gain genes (Additional file 3): %d\n", length(toro_gain)))
} else {
  toro_gain <- character(0)
}

# -----------------------------------------------------------------------------
# Toro-Dominguez 2014 shared signature from Additional File 1, Spreadsheet 1
# -----------------------------------------------------------------------------
sig_csv <- file.path(OUT_DIR, "toro2014_sheet1_shared_signature.csv")
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
    df$gene <- sub("\\*$", "", df$gene)
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
  cat(sprintf("Toro shared signature (Add. file 1, sheet 1): %d up, %d down after asterisk filter\n",
              length(toro_up), length(toro_down)))
  cat(sprintf("  (including asterisk-flagged: %d up, %d down)\n",
              length(toro_up_all), length(toro_down_all)))
} else {
  stop("Spreadsheet 1 CSV not found: ", sig_csv)
}

# (legacy OCR list retained as fallback; no longer used)
toro_up_ocr <- c(
  "HERC5","RTP4","DDX58","SRBD1","PSMB9","SAP30","VAMP5","CCNA2","AZI2","MT1E",
  "C12orf4","GTPBP2","AIM2","BRCA1","BAZ1A","C1GALT1","AURKB","OIP5","TRIM22",
  "TMX1","TRAFD1","HESX1","NMI","BLVRA","PSMA6","IFITM1","CCNB2","LGALS9",
  "ATP6V1C1","CISH","TOR1B","C1GALT1C1","ETNK1","PKD2","TNFSF10","HOXB7",
  "C21orf91","SCO2","MYD88","PRDM1","MRPL15","KIAA0101","GNG5","CHMP5","CASP1",
  "PRDX4","SP110","MELK","LDHA","RACGAP1","SRGAP2","OASL","SAR1B","LRRN1",
  "CBFB","GBP1","KPNB1","PSMA4","TNFAIP6","TRAF3","RSAD2","IFI16","PML","RECQL",
  "NDUFAF1","STOM","TYMP","PSMA2","CDK1","GLRX","PLSCR1","UCHL3","WSB2",
  "ACOT13","IFITM3","BARD1","UGP2","ANXA4","TMEM165","RALA","DLAT","VRK2",
  "BATF","NDUFV2","COMMD3","SPTLC2","RNF170","RARS","RNF31","CAV1","CLDN5",
  "TDRD7","MTDH","CCP110","IFI35","GINS2","PTPLA","IFIT1","UFD1","MUC1",
  "C20orf21","DAPP1","POLE2","EFR3A","H1F0","FAM13A","CDC45","ATP5H","MINPP1",
  "ARFIP1","RPA3","C3orf14","CD69","IGF2BP3","MRPL40","REXO2","CDC20","OBFC2A",
  "COX7B","FAS","ALG6","C11orf73","DLGAP5","PSMC2","RGL1","HJURP","IFI44",
  "ICT1","CMC2","RHOT1","HMMR","C4orf27","TK1","TTK","ISG15","RPS27L","RFK",
  "GLRX2","SIGLEC1","UGCG","TIPRL","CD58","NDUFAB1","BUB1","IFI44L","TRIP6",
  "ARL5A","FKBP3","CTBS","ACP2","ERP44","EDEM3","CD164","IL15RA","USP18",
  "ERGIC2","COPS5","PSMG1","CENPA","UBE2C","CENPM","SQRDL","SUZ12","NDUFA4",
  "STK38","C20orf4","DRAP1","SP140","MCTS1","SSR1","HIGD1A","STRADB","OTOF",
  "PSMD14","RBCK1","RBX1","CDKN3","CAPZA2","NUSAP1","RPL26L1","TRMT2A","ABCA1",
  "SRPK1","IRF7","EIF2AK2","STAT1","MX1","MX2","OAS1","OAS2","OAS3","IFI27"
)

toro_down_ocr <- c(
  "EIF3F","TMEM8B","SLC41A3","EIF4B","ADCK3","MARCKSL1","FCGBP","RPL29",
  "ARFGAP3","PRPF8","PPP1R9B","HNRNPA0","RPS4X","C11orf2","FAM160B","PGAP3",
  "SCARB1","PCBP2","ZNF250","CCNI","CAMK2N1","MAPKAPK3","SNTA1","RPL3L",
  "MEF2D","PPP1R13B","ZNF274","TCTN1","RAB40C","HABP4","TRIM25","CCDC101",
  "EIF3L","MEPCE","FBLN2","C12orf10","COL6A2","TIMM22","MAST2","SSBP2","CRIP2",
  "WNT2B","FOXO1","MED25","KDM1A","EIF3D","MAP3K4","RRP8","SLC7A6","TNFRSF10B",
  "AUTS2","ATP5G2","PIK3IP1","GRM2","SH3PXD2A","ASTE1","POLR1E","RPL19",
  "CHMP7","DNAJB1","TGFBR2","NDEL1","ARFGAP2","TNRC6B","SPATA20","CYTH3","HDC",
  "NMB","RASGRP2","MECP2","USE1","UBE2N","ARHGAP1","ARHGAP4","MAN1B1","SMCR7L",
  "ARHGAP2","MRC2","LFRN3","PHYH","OSGEP","PEX16","PBX2","GLTSCR2","SNRPN",
  "NEO1","SFMBT1","NR3C2","SMG9","LAPTM5","CAMKK2","TRAPPC4","S100A5","SH3BP4",
  "MDH2","MED24","FXYD5","LRRC23","FBXW4","CCR9","SCML1","FBLIM1","TTK",
  "EEF2","RPL32","RPS15","RPS16","RPL9","RPL10","RPS27","EIF3G","EIF3H",
  "EIF3I","EIF3J"
)

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
    "Toro GAIN genes (Add. file 3)",
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
# Toro GO categories (Additional file 1, Spreadsheet 4) — per-pathway gene recovery
# -----------------------------------------------------------------------------
go_csv <- file.path(OUT_DIR, "toro2014_sheet4_go_enrichment.csv")
if (file.exists(go_csv)) {
  go_lines <- readLines(go_csv)
  in_up <- FALSE; in_down <- FALSE
  go_rows <- list()
  for (line in go_lines) {
    if (grepl("^Pathways related with Up", line)) { in_up <- TRUE;  in_down <- FALSE; next }
    if (grepl("^Pathways related with Down", line)) { in_up <- FALSE; in_down <- TRUE;  next }
    if (grepl("^Genes;", line) || !nzchar(trimws(line)) || grepl("^;+$", line)) next
    flds <- strsplit(line, ";", fixed = TRUE)[[1]]
    if (length(flds) < 4) next
    genes_str <- paste(flds[4:length(flds)], collapse = " ")
    genes <- trimws(unlist(strsplit(genes_str, "[ ,]+")))
    genes <- genes[nzchar(genes) & genes != "\""]
    go_rows[[length(go_rows) + 1]] <- tibble(
      direction    = if (in_up) "up" else "down",
      term         = gsub("\"", "", flds[3]),
      n_genes_toro = length(genes),
      toro_genes   = paste(genes, collapse = ";")
    )
  }
  go_df <- bind_rows(go_rows)

  sle_up   <- unique(sig(sle,"up")$gene);   sle_down <- unique(sig(sle,"down")$gene)
  sjs_up   <- unique(sig(sjs,"up")$gene);   sjs_down <- unique(sig(sjs,"down")$gene)

  pw_recov <- go_df |>
    rowwise() |>
    mutate(
      genes_vec = list(strsplit(toro_genes, ";", fixed = TRUE)[[1]]),
      n_in_sle_concordant = sum(genes_vec %in% if (direction == "up") sle_up else sle_down),
      n_in_sjs_concordant = sum(genes_vec %in% if (direction == "up") sjs_up else sjs_down),
      n_shared_concordant = sum(genes_vec %in% if (direction == "up") shared_up else shared_down),
      recovered_sle = paste(intersect(genes_vec, if (direction == "up") sle_up else sle_down), collapse = ";"),
      recovered_shared = paste(intersect(genes_vec, if (direction == "up") shared_up else shared_down), collapse = ";")
    ) |>
    ungroup() |>
    select(direction, term, n_genes_toro, n_in_sle_concordant, n_in_sjs_concordant,
           n_shared_concordant, recovered_shared)

  write_tsv(pw_recov, file.path(OUT_DIR, "toro_go_pathway_recovery.tsv"))
  cat("\nToro GO pathway recovery (top by fraction in shared):\n")
  print(pw_recov |>
          mutate(frac_shared = round(n_shared_concordant / n_genes_toro, 2)) |>
          arrange(desc(frac_shared)) |>
          select(direction, term, n_genes_toro, n_shared_concordant, frac_shared) |>
          head(15), n = 15)
}

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

# Toro's enriched terms (top ones from sheet 4)
toro_terms_up <- c(
  "GO:0000278", "GO:0019221", "GO:0009615", "GO:0060337", "GO:0044419",
  "GO:0000075", "GO:0006915", "GO:0006260", "GO:0045087", "GO:0051437",
  "GO:0006281", "GO:0002474", "GO:0060333", "GO:0007165", "GO:0007259",
  "GO:0051607"
)
toro_terms_down <- c(
  "GO:0010467", "GO:0044267", "GO:0006412", "GO:0006414", "GO:0006415",
  "GO:0006413"
)

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
