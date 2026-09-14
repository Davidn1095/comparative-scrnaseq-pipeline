#!/usr/bin/env Rscript
################################################################################
# 10_manuscript_figures.R - Compose figures for Manuscript 1
#
# Manuscript 1: Comparative single-cell transcriptomic landscape of SLE and SjS
# Target journal: Nature Communications
#
# Manuscript figures (180 mm / 7.087 in double-column width):
#   Figure 1: Atlas Overview              (build_fig1)
#   Figure 2: Disease-Specific Programs   (build_fig2)
#   Figure 3: Pathway Enrichment          (build_fig3)
#   Figure 4: Donor-level Classification  (build_fig4)
#
# Output: results/manuscript1_figures/
################################################################################

# --- Load configuration ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  return(file.path(Sys.getenv("ATLAS_ROOT", normalizePath("~/projects/autoimmune-atlas")), "scripts"))
}
SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
source(file.path(SCRIPT_DIR, "00_utils.R"))

# --- Load packages ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

# Optional packages
has_magick <- requireNamespace("magick", quietly = TRUE)
has_cowplot <- requireNamespace("cowplot", quietly = TRUE)
has_png <- requireNamespace("png", quietly = TRUE)
has_ComplexHeatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)

if (has_magick) library(magick)
if (has_cowplot) library(cowplot)
if (has_ComplexHeatmap) {
  library(ComplexHeatmap)
  library(circlize)
}

has_Seurat <- requireNamespace("Seurat", quietly = TRUE)
if (has_Seurat) suppressPackageStartupMessages(library(Seurat))

if (!requireNamespace("ggrastr", quietly = TRUE)) {
  warning("ggrastr not installed — UMAP rasterisation will be skipped")
}

# Lineage groupings (LINEAGE_MAP, LINEAGE_ORDER, LINEAGE_COLORS) come from 00_config.R.

# Disease colors used by ComplexHeatmap annotations in figure 3 (pathways)
disease_colors <- c("SLE" = "#D55E00", "SjS" = "#0072B2")

# --- Configuration ---
ATLAS_OVERVIEW <- file.path(ATLAS_ROOT, "results", "plots", "atlas_overview")
DIFFERENTIAL_EXPRESSION <- file.path(ATLAS_ROOT, "results", "plots", "differential_expression")
OUTDIR <- file.path(ATLAS_ROOT, "results", "manuscript1_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Nature Communications figure specifications
# Double-column: 180 mm = 7.087 in; single-column: 88 mm = 3.465 in
# Text: 5-7 pt body, 8 pt bold lowercase panel labels
# Font: Helvetica or Arial (sans-serif)
# DPI: 300 min, 450 recommended
FIG_WIDTH        <- 180 / 25.4   # 7.087 in (Nature double-column max)
FIG_WIDTH_SINGLE <- 88 / 25.4    # 3.465 in (Nature single-column)
BASE_SIZE  <- 7                   # Nature max text size
LABEL_SIZE <- 8                   # Nature panel label size
FONT_FAMILY <- "Helvetica"        # Nature required sans-serif

# =============================================================================
# Helper functions
# =============================================================================

load_img <- function(path) {
  if (!file.exists(path)) {
    log_msg("  WARNING: ", path, " not found")
    return(NULL)
  }
  if (file.info(path)$size == 0) {
    log_msg("  WARNING: ", path, " is empty (0 bytes)")
    return(NULL)
  }
  if (has_magick) {
    img <- tryCatch(
      magick::image_read(path),
      error = function(e) {
        log_msg("  WARNING: Failed to read ", path, ": ", e$message)
        NULL
      }
    )
    return(img)
  } else if (has_png && grepl("\\.png$", path, ignore.case = TRUE)) {
    img <- tryCatch(
      png::readPNG(path),
      error = function(e) {
        log_msg("  WARNING: Failed to read ", path, ": ", e$message)
        NULL
      }
    )
    return(img)
  }
  log_msg("  WARNING: Cannot load image (install magick or png package)")
  NULL
}

img_to_grob <- function(img) {
  if (is.null(img)) {
    return(grid::textGrob("Image not found", gp = grid::gpar(col = "grey50")))
  }
  if (has_magick && inherits(img, "magick-image")) {
    return(grid::rasterGrob(as.raster(img), interpolate = TRUE))
  } else if (is.array(img)) {
    return(grid::rasterGrob(img, interpolate = TRUE))
  }
  grid::textGrob("Cannot display image", gp = grid::gpar(col = "grey50"))
}

# Save figure as PDF (vector, Nature Communications preferred)
save_nature_fig <- function(plot, outpath, width = FIG_WIDTH, height = 8) {
  pdf_path <- sub("\\.png$", ".pdf", outpath)
  ggsave(pdf_path, plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  log_msg("  Saved: ", pdf_path)
}

theme_nature <- function(base_size = BASE_SIZE, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, color = "grey40"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ---- Figure 4 helpers and palette ----
COND_COLOR   <- c(healthy = "#009E73", sle = "#D55E00", sjs = "#0072B2")
COND_LABEL   <- c(healthy = "HC",      sle = "SLE",     sjs = "SjS")
CLASS_FILL   <- setNames(unname(COND_COLOR), unname(COND_LABEL))  # fig 4 stacked bars

# ---- Bar convention ----
# BAR_RADIUS is the corner radius of EVERY bar in the manuscript, stacked or
# single (Fig 2a, the Jaccard bars under 2b and 3c). Single bars are plain
# chicklets: fill and outline are one shape and cannot invert.
# Stacked bars (Figs 1e, 2b top, 3c top, 4b, 4c) use stack_outline(): plain
# geom_col segments with no stroke under ONE no-fill rounded outline per bar,
# so the rounding lives on the silhouette and a thin segment can never invert
# (a per-segment chicklet inverts once the segment is shorter than its corner
# radius). The radius is 0.75 pt: the top segment's square corners poke
# 0.29 x r beyond the arc, which a 0.2 stroke hides up to r ~ 1 pt but not at
# 1.5 pt. A rounded outline inverts the same way when the WHOLE bar is shorter
# than ~2 x r, so bars below STACK_MIN_FRAC of the tallest bar get a square
# outline instead, indistinguishable at that height.
BAR_RADIUS   <- grid::unit(0.75, "pt")
STACK_MIN_FRAC <- 0.03
stack_outline <- function(tot, x, y, width, size = 0.2) {
  cut   <- STACK_MIN_FRAC * max(tot[[y]])
  tall  <- tot[tot[[y]] >= cut, , drop = FALSE]
  short <- tot[tot[[y]] <  cut, , drop = FALSE]
  list(
    ggchicklet::geom_chicklet(data = tall, aes(x = .data[[x]], y = .data[[y]]),
                              inherit.aes = FALSE, width = width, fill = NA,
                              color = "black", size = size, radius = BAR_RADIUS),
    if (nrow(short) > 0)
      geom_col(data = short, aes(x = .data[[x]], y = .data[[y]]),
               inherit.aes = FALSE, width = width, fill = NA,
               color = "black", linewidth = size)
  )
}

clean_ct <- function(x) gsub("_", " ", x)

# Shared pathway-label formatter used by build_fig3 and build_fig4. Strips
# HALLMARK_, swaps _ for space, title-cases and re-uppercases known acronyms
# so the same Hallmark ID renders identically across figures.
format_pathway_label <- function(x) {
  s <- sub("^HALLMARK_", "", x); s <- tolower(gsub("_", " ", s))
  s <- tools::toTitleCase(s)
  acro <- c("Ifn"="IFN", "Il2"="IL2", "Il6"="IL6", "Uv"="UV", "Dna"="DNA",
            "Myc"="MYC", "Mtorc1"="MTORC1", "Mtor"="MTOR", "Pi3k"="PI3K",
            "Akt"="AKT", "Tnfa"="TNFA", "Nfkb"="NFKB", "Tgf"="TGF",
            "E2f"="E2F", "G2m"="G2M", "P53"="P53", "Wnt"="WNT",
            "Stat3"="STAT3", "Stat5"="STAT5", "Jak"="JAK", "Kras"="KRAS",
            "Ros"="ROS", "V1"="V1", "V2"="V2")
  for (k in names(acro)) s <- gsub(paste0("\\b", k, "\\b"), acro[k], s)
  s
}
# Back-compat alias for the three call sites already using clean_path.
clean_path <- format_pathway_label

# =============================================================================
# Figure 1: Pipeline + Atlas Overview (6 panels, 4 rows)
#   a: Pipeline overview (external PDF from TikZ/Overleaf)
#   b: UMAP by dataset (accession)
#   c: UMAP by condition (Healthy / SjS / SLE)
#   d: UMAP by lineage (8 groups) + centroid labels
#   e: Cell lineage composition by condition
#   f: Per-donor cell-type composition (propeller bubbles)
# =============================================================================

# Helper: bubble plot of per-cell-type differential abundance results.
# Reads the three pairwise propeller TSVs written by 05_composition.R and
# renders cell types on the x-axis ordered by Fig 2's DEG count.
build_composition_bubble <- function() {
  comp_dir <- file.path(ATLAS_ROOT, "results", "composition")
  load_one <- function(fn, label) {
    df <- read.delim(file.path(comp_dir, fn), stringsAsFactors = FALSE)
    data.frame(
      cell_type   = df$cell_type,
      comparison  = label,
      log2_fc     = df$log2_fc,
      fdr         = df$fdr,
      significant = df$significant == "TRUE" | df$significant == TRUE,
      stringsAsFactors = FALSE
    )
  }
  d <- rbind(
    load_one("propeller_HC_vs_SLE.tsv",  "HC vs SLE"),
    load_one("propeller_HC_vs_SjS.tsv",  "HC vs SjS"),
    load_one("propeller_SLE_vs_SjS.tsv", "SLE vs SjS")
  )
  d$cell_type <- gsub("[-/]", " ", d$cell_type)

  # Cell-type order: same as Fig 2 panel A (DEG count descending).
  sig_file <- file.path(DIFFERENTIAL_EXPRESSION, "signatures",
                        "disease_signatures_all.tsv")
  ct_order_fig2 <- character(0)
  if (file.exists(sig_file)) {
    ct_order_fig2 <- read.delim(sig_file, stringsAsFactors = FALSE) %>%
      mutate(disease = case_when(
        disease %in% c("DE_SLE", "SLE") ~ "SLE",
        disease %in% c("DE_SJS", "SjS") ~ "SjS",
        TRUE ~ NA_character_
      )) %>%
      filter(disease %in% c("SjS", "SLE"), !is.na(padj),
             padj < 0.05, abs(log2FoldChange) >= 1.0) %>%
      count(cell_type, name = "n_total") %>%
      arrange(desc(n_total)) %>%
      pull(cell_type)
  }
  ct_levels <- c(ct_order_fig2,
                 setdiff(unique(d$cell_type), ct_order_fig2))
  d$cell_type <- factor(d$cell_type, levels = ct_levels)
  d$comparison <- factor(d$comparison,
                         levels = rev(c("HC vs SLE", "HC vs SjS", "SLE vs SjS")))

  SIZE_CAP <- 15
  d$neg_log10_fdr <- pmin(-log10(d$fdr), SIZE_CAP)

  dominant_group <- function(cmp, lfc) {
    mapply(function(c, l) {
      if (c == "HC vs SLE")  { if (l > 0) "SLE" else "HC" }
      else if (c == "HC vs SjS") { if (l > 0) "SjS" else "HC" }
      else                       { if (l > 0) "SjS" else "SLE" }
    }, cmp, lfc, USE.NAMES = FALSE)
  }
  d$dominant <- dominant_group(as.character(d$comparison), d$log2_fc)
  d$dominant <- ifelse(d$significant, d$dominant, NA_character_)
  d$dominant <- factor(d$dominant, levels = c("HC", "SLE", "SjS"))

  ggplot(d, aes(x = cell_type, y = comparison)) +
    geom_point(data = d[!d$significant, ],
               aes(size = neg_log10_fdr),
               fill = "white", colour = "grey55", shape = 21, stroke = 0.25) +
    geom_point(data = d[d$significant, ],
               aes(size = neg_log10_fdr, fill = dominant),
               shape = 21, colour = "black", stroke = 0.3) +
    scale_x_discrete(limits = ct_levels) +
    scale_size_continuous(
      name = expression(-log[10] * "(FDR)"),
      range = c(0.4, 3.5),
      breaks = c(2, 5, 10, 15),
      labels = c("2", "5", "10", expression("" >= "15")),
      limits = c(0, SIZE_CAP)
    ) +
    scale_fill_manual(
      name = "Higher in",
      values = c("HC" = "#009E73", "SLE" = "#D55E00", "SjS" = "#0072B2"),
      labels = c("HC", "SLE", "SjS"),
      na.value = "white", drop = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    theme_nature(base_size = 7) +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7, face = "bold"),
      panel.grid.major.x = element_line(color = "grey95", linewidth = 0.2),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(1, "pt"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5),
      legend.key.size = unit(0.2, "cm"),
      plot.margin = margin(2, 2, 2, 10)
    )
}

build_fig1 <- function() {
  log_msg("--- M1 Figure 1: Atlas Overview ---")

  if (!has_Seurat) {
    log_msg("  Seurat not available — skipping Figure 1 (requires atlas)")
    return()
  }

  # Load atlas and keep only active (non-excluded) datasets
  atlas <- load_atlas()

  # Active datasets: 3 SLE + 2 SjS
  active_datasets <- basename(list.dirs(
    file.path(ATLAS_ROOT, "data", "SLE"), recursive = FALSE
  ))
  active_datasets <- c(active_datasets, basename(list.dirs(
    file.path(ATLAS_ROOT, "data", "SjS"), recursive = FALSE
  )))
  log_msg("  Active datasets: ", paste(active_datasets, collapse = ", "))

  atlas <- subset(atlas, dataset_id %in% active_datasets)

  if (CONDITIONS_EXCLUDE != "") {
    conds_excl <- trimws(strsplit(CONDITIONS_EXCLUDE, ",")[[1]])
    atlas <- subset(atlas, !condition %in% conds_excl)
  }
  atlas <- map_plot_condition(atlas)
  log_msg("  Atlas: ", ncol(atlas), " cells")

  # UMAP coordinates: scVI-latent UMAP (trained on full 869k atlas x 2000 HVGs,
  # n_latent=30, see Methods). Cells absent from scVI output (non-common cell
  # types) are filtered out by the join.
  scvi_path <- file.path(ATLAS_ROOT, "results", "integration",
                         "scvi_full_gpu", "scvi_umap.tsv")
  scvi_df <- read.delim(scvi_path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(scvi_df)[1] <- "barcode"
  atlas_bcs <- colnames(atlas)
  match_idx <- match(atlas_bcs, scvi_df$barcode)
  umap_df <- data.frame(
    UMAP_1     = scvi_df$UMAP_1[match_idx],
    UMAP_2     = scvi_df$UMAP_2[match_idx],
    condition  = atlas$plot_condition,
    cell_type  = as.character(atlas$cell_type),
    dataset_id = atlas$dataset_id,
    row.names  = atlas_bcs,
    stringsAsFactors = FALSE
  )
  umap_df <- umap_df[!is.na(umap_df$UMAP_1), ]
  log_msg("  scVI UMAP joined: ", nrow(umap_df), " cells (",
          length(atlas_bcs) - nrow(umap_df), " atlas cells absent from scVI)")

  # Restrict to the 25 common cell types used for all downstream analyses (DE,
  # GSEA, composition, classifier). Drops plasmacytoid DCs, plasmablasts,
  # low-density basophils and low-density neutrophils so that every panel of
  # Fig 1 b-e visualises the same analysis cohort.
  umap_df <- umap_df[is_common_celltype(umap_df$cell_type), ]
  log_msg("  Atlas restricted to 25 common cell types: ", nrow(umap_df), " cells")

  # Assign lineage. After the 25-cell-type filter, the "Other" lineage in
  # LINEAGE_MAP contains only Progenitor cells, so relabel for clarity. The
  # Plasmablasts lineage is empty after filtering and is dropped from the
  # 8-lineage order used by panels d and e.
  umap_df$lineage <- LINEAGE_MAP[umap_df$cell_type]
  umap_df$lineage[umap_df$lineage == "Other"] <- "Progenitor"
  fig1_lineage_order <- c("CD4 T", "CD8 T", "Unconventional T", "B cells",
                          "NK cells", "Monocytes", "DCs", "Progenitor")
  fig1_lineage_colors <- setNames(
    c(unname(LINEAGE_COLORS[c("CD4 T", "CD8 T", "Unconventional T", "B cells",
                              "NK cells", "Monocytes", "DCs")]),
      unname(LINEAGE_COLORS["Other"])),
    fig1_lineage_order
  )
  umap_df$lineage <- factor(umap_df$lineage, levels = fig1_lineage_order)

  # Shuffle to avoid overplotting bias
  set.seed(42)
  umap_df <- umap_df[sample(nrow(umap_df)), ]

  # Dataset colour palette. The CXG identifier is kept in full (no shortening
  # anywhere in the manuscript); it is broken over two lines for the legend.
  dataset_labels <- setNames(
    ifelse(grepl("^CXG_", umap_df$dataset_id),
           sub("^(CXG_[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-)", "\\1\n", umap_df$dataset_id),
           umap_df$dataset_id),
    umap_df$dataset_id
  )
  umap_df$dataset_label <- dataset_labels[umap_df$dataset_id]
  datasets <- sort(unique(umap_df$dataset_label))
  n_ds <- length(datasets)
  dataset_colors <- setNames(scales::hue_pal()(n_ds), datasets)

  # ---- Panel B: UMAP by dataset ----
  panel_b <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = dataset_label)) +
    ggrastr::rasterise(geom_point(size = 0.2, alpha = 0.4, stroke = 0), dpi = 300) +
    scale_color_manual(values = dataset_colors, name = "Dataset") +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_nature() +
    theme(legend.position = c(0.82, 0.15),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 6),
          axis.text = element_text(size = 7),
          plot.margin = margin(12, 12, 12, 12)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                ncol = 1))

  # ---- Panel C: UMAP by condition ----
  panel_c <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
    ggrastr::rasterise(geom_point(size = 0.2, alpha = 0.4, stroke = 0), dpi = 300) +
    scale_color_manual(values = CONDITION_COLORS, name = NULL) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_nature() +
    theme(legend.position = c(0.82, 0.15),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 7),
          axis.text = element_text(size = 7),
          plot.margin = margin(12, 12, 12, 12)) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  # ---- Panel D: UMAP by lineage (8 groups) + centroid labels ----
  centroids <- umap_df %>%
    group_by(lineage) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")

  panel_d <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = lineage)) +
    ggrastr::rasterise(geom_point(size = 0.2, alpha = 0.4, stroke = 0), dpi = 300) +
    ggrepel::geom_text_repel(
      data = centroids, aes(x = x, y = y, label = lineage),
      size = 2.5, fontface = "bold", color = "grey20",
      bg.color = "white", bg.r = 0.15,
      box.padding = 0.4, point.padding = 0.2,
      min.segment.length = 0.3, seed = 42, max.overlaps = 20,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = fig1_lineage_colors, name = NULL) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_nature() +
    theme(legend.position = "none",
          axis.text = element_text(size = 7),
          plot.margin = margin(12, 12, 12, 12))

  # ---- Panel E: Composition by condition (8 lineage groups) ----
  comp_df <- as.data.frame(table(umap_df$condition, umap_df$lineage))
  names(comp_df) <- c("condition", "lineage", "n")
  comp_df <- comp_df %>%
    filter(n > 0) %>%
    group_by(condition) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  comp_df$lineage <- factor(comp_df$lineage, levels = rev(fig1_lineage_order))

  comp_tot <- aggregate(pct ~ condition, comp_df, sum)
  panel_e <- ggplot(comp_df, aes(x = condition, y = pct)) +
    geom_col(aes(fill = lineage), width = 0.7, color = NA,
             position = position_stack(reverse = TRUE)) +
    stack_outline(comp_tot, "condition", "pct", width = 0.7, size = 0.3) +
    scale_fill_manual(values = fig1_lineage_colors, name = NULL,
                      breaks = rev(fig1_lineage_order)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_nature() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          axis.text = element_text(size = 7),
          plot.margin = margin(12, 12, 12, 12)) +
    guides(fill = guide_legend(nrow = 2, reverse = TRUE,
                               override.aes = list(color = NA)))

  # ---- Panel F: per-donor cell-type composition (propeller bubbles) ----
  panel_f <- build_composition_bubble()

  # ---- Compose: 3 rows (b+c / d+e / f) — panel a is pipeline PDF, merged after
  row1 <- cowplot::plot_grid(
    panel_b, panel_c,
    ncol = 2,
    labels = c("b", "c"), label_size = LABEL_SIZE, label_fontface = "bold"
  )
  row2 <- cowplot::plot_grid(
    panel_d, panel_e,
    ncol = 2, rel_widths = c(1, 1),
    labels = c("d", "e"), label_size = LABEL_SIZE, label_fontface = "bold"
  )
  row3 <- cowplot::plot_grid(
    panel_f,
    ncol = 1,
    labels = c("f"), label_size = LABEL_SIZE, label_fontface = "bold"
  )
  fig <- cowplot::plot_grid(
    row1, row2, row3,
    ncol = 1, rel_heights = c(1, 1, 1)
  )

  panels_pdf <- file.path(OUTDIR, "M1_Figure1_panels_bcef.pdf")
  save_nature_fig(fig, panels_pdf, width = FIG_WIDTH, height = 8.0)

  # Final Figure 1 is assembled by fig1_combined.tex (pipeline TikZ + this PDF)
  combined_pdf <- file.path(OUTDIR, "M1_Figure1_atlas_overview.pdf")
  if (!file.exists(combined_pdf)) {
    file.copy(panels_pdf, combined_pdf, overwrite = TRUE)
    log_msg("  WARNING: pipeline PDF not found, panels only → ", combined_pdf)
  }

  # Free memory
  rm(atlas, umap_df, comp_df)
  gc(verbose = FALSE)
}

# =============================================================================
# Figure 2: Disease-Specific Programs (3 panels)
#   a: DEG counts barplot (per cell type, split by disease and direction)
#   b: Shared vs disease-specific DEGs stacked bar
#   c: Jaccard similarity (SLE vs SjS DEG overlap per cell type)
# =============================================================================

build_fig2 <- function() {
  log_msg("--- M1 Figure 2: Disease-Specific Programs ---")

  # ---- Load DEG signatures (shared by panels A and C) ----
  panel_a <- NULL
  sig_file <- file.path(DIFFERENTIAL_EXPRESSION, "signatures", "disease_signatures_all.tsv")

  sigs <- NULL
  if (file.exists(sig_file)) {
    sigs <- read.delim(sig_file, stringsAsFactors = FALSE) %>%
      mutate(disease = case_when(
        disease == "DE_SJS" ~ "SjS",
        disease == "DE_SLE" ~ "SLE",
        TRUE ~ disease
      )) %>%
      filter(disease %in% c("SjS", "SLE"), !is.na(padj),
             padj < 0.05, abs(log2FoldChange) >= 1.0)
  }

  # Canonical cell-type ordering used across ALL three fig 2 panels: total
  # DEG count per cell type (descending). Computed once here so panels a/b/c
  # share identical column order; readers can scan a single cell type
  # vertically across the figure.
  ct_order <- NULL
  if (!is.null(sigs) && nrow(sigs) > 0) {
    ct_order <- sigs %>%
      count(cell_type, name = "n_total") %>%
      arrange(desc(n_total)) %>%
      pull(cell_type)
  }

  # ---- Panel A: DEG counts barplot (vertical bars, cell types on X-axis) ----
  if (!is.null(sigs)) {
    tryCatch({

      # Count DEGs by cell type, disease, and direction
      # Make down-regulated negative for diverging bar chart
      deg_counts <- sigs %>%
        mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
        group_by(cell_type, disease, direction) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(
          n_plot = ifelse(direction == "Down", -n, n),
          fill_group = paste(disease, direction, sep = "_")
        )

      # Order cell types using the canonical ct_order computed at function
      # scope (shared with panels b and c).
      deg_counts$cell_type <- factor(deg_counts$cell_type, levels = ct_order)
      deg_counts$disease <- factor(deg_counts$disease, levels = c("SLE", "SjS"))
      deg_counts$fill_group <- factor(deg_counts$fill_group,
                                       levels = c("SLE_Down", "SLE_Up", "SjS_Down", "SjS_Up"))

      # Colors: SLE = warm (orange family), SjS = cool (blue family)
      # Up = saturated, Down = lighter variant of same hue
      fill_colors <- c(
        "SLE_Up"   = "#D55E00",
        "SLE_Down" = "#E69F00",
        "SjS_Up"   = "#0072B2",
        "SjS_Down" = "#56B4E9"
      )

      # Calculate x-axis limits for symmetry
      x_max <- max(abs(deg_counts$n_plot)) * 1.1

      # Vertical diverging barplot: cell types on x-axis
      panel_a <- ggplot(deg_counts, aes(x = cell_type, y = n_plot, fill = fill_group,
                                         group = disease)) +
        ggchicklet::geom_chicklet(position = position_dodge(width = 1.0, preserve = "single"),
                                  color = "black", size = 0.15, width = 0.95,
                                  radius = BAR_RADIUS) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        scale_fill_manual(
          values = fill_colors,
          labels = c("SLE_Up" = "SLE Up", "SLE_Down" = "SLE Down",
                     "SjS_Up" = "SjS Up", "SjS_Down" = "SjS Down"),
          name = NULL,
          guide = guide_legend(nrow = 1, override.aes = list(color = NA))
        ) +
        scale_y_continuous(
          limits = c(-x_max, x_max),
          labels = function(y) abs(y),
          expand = expansion(mult = c(0.02, 0.02))
        ) +
        annotate("text", x = 0.5, y = x_max * 0.85, label = "Up",
                 hjust = 0, size = 2.5, fontface = "bold") +
        annotate("text", x = 0.5, y = -x_max * 0.85, label = "Down",
                 hjust = 0, size = 2.5, fontface = "bold") +
        labs(x = NULL, y = "Number of DEGs",
             title = "DEG Counts") +
        theme_nature(base_size = 8) +
        theme(
          axis.text.x  = element_text(size = 5, angle = 45, hjust = 1),
          plot.margin  = margin(5, 5, 5, 30),
          legend.position = "top",
          legend.justification = "right",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major.x = element_blank()
        )
    }, error = function(e) {
      log_msg("  WARNING: DEG counts panel failed: ", e$message)
    })
  }

  if (is.null(panel_a)) {
    panel_a <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "DEG data not found")
  }

  # ---- Panel B: Jaccard similarity (cross-disease DEG overlap) per cell type ----
  # Higher similarity = more shared DEGs between SLE and SjS
  panel_b <- NULL
  sim_dir <- file.path(DIFFERENTIAL_EXPRESSION, "similarity")
  if (dir.exists(sim_dir)) {
    tryCatch({
      # Canonical files only: similarity_<ct>.tsv. Stale legacy copies from
      # earlier runs (similarity_<ct>_annotated.tsv and double-annotated)
      # would produce duplicate bars per cell type.
      sim_files <- list.files(sim_dir,
                              pattern = "^similarity_[^.]+\\.tsv$",
                              full.names = TRUE)
      sim_files <- sim_files[!grepl("_annotated", basename(sim_files))]
      sim_data <- list()

      for (f in sim_files) {
        ct <- gsub("similarity_", "", basename(f))
        ct <- gsub("\\.tsv$", "", ct)
        ct <- gsub("_", " ", ct)

        df <- read.delim(f, row.names = 1, stringsAsFactors = FALSE)
        if ("SjS" %in% colnames(df) && "SLE" %in% rownames(df)) {
          sim_data[[ct]] <- df["SLE", "SjS"]
        } else if (ncol(df) >= 2 && nrow(df) >= 1) {
          sim_data[[ct]] <- df[1, 2]
        }
      }

      if (length(sim_data) > 0) {
        sim_df <- data.frame(
          cell_type = names(sim_data),
          jaccard = unlist(sim_data),
          stringsAsFactors = FALSE
        )
        # Apply canonical ordering (panel-a/c match). Fall back to Jaccard-
        # descending if ct_order is unavailable for any reason.
        if (!is.null(ct_order)) {
          sim_df$cell_type <- factor(sim_df$cell_type, levels = ct_order)
          sim_df <- sim_df[order(sim_df$cell_type), ]
        } else {
          sim_df <- sim_df %>%
            arrange(desc(jaccard)) %>%
            mutate(cell_type = factor(cell_type, levels = cell_type))
        }

        panel_b <- ggplot(sim_df, aes(x = cell_type, y = jaccard)) +
          ggchicklet::geom_chicklet(fill = "#009E73", color = "black", size = 0.3,
                                    radius = BAR_RADIUS) +
          scale_y_continuous(limits = c(0, max(sim_df$jaccard) * 1.05),
                             expand = c(0, 0)) +
          labs(y = "Jaccard similarity", x = NULL) +
          theme_nature(base_size = 8) +
          theme(
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
            plot.margin = margin(5, 5, 5, 30),
            panel.grid.major.x = element_blank()
          )
      }
    }, error = function(e) {
      log_msg("  WARNING: Jaccard panel failed: ", e$message)
    })
  }

  if (is.null(panel_b)) {
    panel_b <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Jaccard data not found")
  }

  # ---- Panel C: Shared vs disease-specific DEGs stacked bar plot ----
  # Source: full DESeq2 outputs (same as Jaccard panel) — direction-resolved.
  # Categories: SLE-specific, Shared (concordant log2FC sign), Discordant
  # (DE in both with opposite sign), SjS-specific.
  panel_c <- NULL

  tryCatch({
    load_de_full <- function(disease) {
      dir <- file.path(DIFFERENTIAL_EXPRESSION, paste0("de_", disease))
      files <- list.files(dir, pattern = "^deseq2_.*\\.tsv$", full.names = TRUE)
      files <- files[!grepl("_annotated", basename(files))]
      do.call(rbind, lapply(files, function(f) {
        ct <- sub("^deseq2_", "", sub("\\.tsv$", "", basename(f)))
        df <- tryCatch(read.table(f, header = TRUE, sep = "\t"), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) return(NULL)
        df$cell_type <- ct
        df$disease   <- toupper(disease)
        df
      }))
    }

    full_de <- rbind(load_de_full("sle"), load_de_full("sjs")) %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1.0,
             !grepl(EXCLUDED_GENE_REGEX, gene, ignore.case = TRUE)) %>%
      mutate(cell_type = gsub("_", " ", cell_type),
             direction = ifelse(log2FoldChange > 0, "up", "down"))

    # Supplementary table: full DEG counts per cell type x disease x direction
    # (no top-N cap; same filters as panel c — padj<0.05, |log2FC|>=1.0,
    # MT/RPL/RPS/HB genes excluded).
    deg_table <- full_de %>%
      group_by(cell_type, disease, direction) %>%
      summarise(n = n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = c(disease, direction),
                         values_from = n,
                         values_fill = 0L,
                         names_sep = "_") %>%
      mutate(SLE_Total = SLE_up + SLE_down,
             SJS_Total = SJS_up + SJS_down,
             Total     = SLE_Total + SJS_Total) %>%
      select(cell_type,
             SLE_Up = SLE_up, SLE_Down = SLE_down, SLE_Total,
             SjS_Up = SJS_up, SjS_Down = SJS_down, SjS_Total = SJS_Total,
             Total) %>%
      arrange(desc(Total))
    deg_table_path <- file.path(OUTDIR, "SupplementaryTable_DEG_counts.tsv")
    write.table(deg_table, deg_table_path,
                sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("  Supplementary DEG counts written to: ", deg_table_path)

    cell_types_c <- sort(unique(full_de$cell_type))

    deg_counts_c <- do.call(rbind, lapply(cell_types_c, function(ct) {
      s <- full_de[full_de$cell_type == ct, ]
      sle_g <- s[s$disease == "SLE", c("gene", "direction")]
      sjs_g <- s[s$disease == "SJS", c("gene", "direction")]
      if (nrow(sle_g) > 0) names(sle_g)[2] <- "dir_sle"
      if (nrow(sjs_g) > 0) names(sjs_g)[2] <- "dir_sjs"

      if (nrow(sle_g) == 0 || nrow(sjs_g) == 0) {
        n_concordant <- 0
        n_discordant <- 0
      } else {
        inter <- merge(sle_g, sjs_g, by = "gene")
        n_concordant <- sum(inter$dir_sle == inter$dir_sjs)
        n_discordant <- sum(inter$dir_sle != inter$dir_sjs)
      }
      n_sle_only <- nrow(sle_g) - n_concordant - n_discordant
      n_sjs_only <- nrow(sjs_g) - n_concordant - n_discordant

      data.frame(
        cell_type        = ct,
        `SLE-specific`   = n_sle_only,
        Shared           = n_concordant,
        Discordant       = n_discordant,
        `SjS-specific`   = n_sjs_only,
        check.names      = FALSE
      )
    }))

    # Order using canonical ct_order (shared with panels a, b). Fall back to
    # within-panel total if ct_order is missing.
    deg_counts_c$total <- deg_counts_c$`SLE-specific` + deg_counts_c$Shared +
                          deg_counts_c$`SjS-specific`
    if (!is.null(ct_order)) {
      deg_counts_c$cell_type <- factor(deg_counts_c$cell_type, levels = ct_order)
      deg_counts_c <- deg_counts_c[order(deg_counts_c$cell_type), ]
    } else {
      deg_counts_c <- deg_counts_c[order(-deg_counts_c$total), ]
      deg_counts_c$cell_type <- factor(deg_counts_c$cell_type, levels = deg_counts_c$cell_type)
    }

    # Discordant counts (DE in both diseases with opposite log2FC sign) are
    # included as a fourth category for visual + categorical consistency with
    # Fig 3c, even though counts are universally small (n=4 total across all
    # cell types). Per-cell-type breakdown is also written to a side TSV.
    discordant_tsv <- file.path(OUTDIR, "M1_Figure2_discordant_counts.tsv")
    write.table(
      deg_counts_c[, c("cell_type", "Discordant")],
      discordant_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("  Discordant DEG counts written to: ", discordant_tsv)

    deg_long_c <- deg_counts_c %>%
      select(-total) %>%
      pivot_longer(cols = c("SLE-specific", "Shared", "Discordant", "SjS-specific"),
                   names_to = "category", values_to = "count")
    deg_long_c$category <- factor(deg_long_c$category,
                                  levels = c("SLE-specific", "Shared",
                                             "Discordant", "SjS-specific"))

    cat_colors_c <- c("SLE-specific" = "#D55E00",
                      "Shared"       = "#CC79A7",   # Okabe-Ito reddish purple (palette family); green is HC elsewhere
                      "Discordant"   = "#BBBBBB",   # lighter grey: moves away from Shared #CC79A7 under deuteranopia
                      "SjS-specific" = "#0072B2")

    deg_tot_c <- aggregate(count ~ cell_type, deg_long_c, sum)
    panel_c <- ggplot(deg_long_c, aes(x = cell_type, y = count)) +
      geom_col(aes(fill = category), width = 0.9, color = NA,
               position = position_stack(reverse = TRUE)) +
      stack_outline(deg_tot_c, "cell_type", "count", width = 0.9, size = 0.2) +
      scale_fill_manual(values = cat_colors_c, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(y = "Number of DEGs", x = NULL,
           title = "Shared vs Specific DEGs") +
      theme_nature(base_size = 8) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(5, 5, 0, 30),
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.x = element_blank()
      )
  }, error = function(e) {
    log_msg("  WARNING: Shared/specific DEG plot failed: ", e$message)
  })

  if (is.null(panel_c)) {
    panel_c <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Signature data not found")
  }

  # ---- Compose: panel a (standalone with its own X axis) +
  #               b/c stacked sub-figure with shared X axis ----
  bc_subfig <- cowplot::plot_grid(
    panel_c,
    panel_b,
    ncol = 1,
    rel_heights = c(0.7, 0.6),
    align = "v",
    axis  = "lr"
  )
  fig <- cowplot::plot_grid(
    panel_a,
    bc_subfig,
    ncol = 1,
    rel_heights = c(1, 1),
    labels = c("a", "b"),
    label_size = LABEL_SIZE,
    label_fontface = "bold"
  )

  outpath <- file.path(OUTDIR, "M1_Figure2_disease_programs.pdf")
  save_nature_fig(fig, outpath, width = FIG_WIDTH, height = 6.69)
}


# =============================================================================
# Figure 4: Pathway-level Shared vs Distinct Programs
#   - Divergent plots (common + disease-specific pathways)
#   - Heatmap of pathways x cell types
# =============================================================================


# Helper: lineage-grouped column ordering for cell types
# Lineage priority: T cell, NK cell, B cell, Monocyte, Dendritic, Other
ct_lineage_order <- function(cts) {
  lineage <- dplyr::case_when(
    grepl("monocyte", cts, ignore.case = TRUE) ~ "4_Monocyte",
    grepl("dendritic", cts, ignore.case = TRUE) ~ "5_Dendritic",
    grepl("memory B|Naive B|Plasmablast|Exhausted B|switched B|B cells?$",
          cts, ignore.case = TRUE) ~ "3_B cell",
    grepl("NK|killer", cts, ignore.case = TRUE) ~ "2_NK",
    grepl("T cell|Th1|Th2|Th17|Treg|MAIT|Tfh|regulatory|Follicular|Vd2|gd T|CD4|CD8|effector",
          cts, ignore.case = TRUE) ~ "1_T cell",
    TRUE ~ "6_Other"
  )
  cts[order(lineage, cts)]
}

# Helper: combined SLE/SjS heatmap with hierarchical pathway clustering
# gsea_norm: normalised long-form GSEA data
# pathways: character vector of pathways to include (rows)
# excluded_cts: cell types to drop (rare types failing thresholds)
build_combined_pathway_heatmap <- function(gsea_norm, pathways, excluded_cts,
                                            ct_order = NULL) {
  cts <- setdiff(unique(gsea_norm$cell_type), excluded_cts)
  if (is.null(ct_order)) {
    ct_order <- ct_lineage_order(cts)
  } else {
    # Honour caller's order, append any cell types not covered
    ct_order <- c(intersect(ct_order, cts), setdiff(cts, ct_order))
  }

  # Build (pathways x cell_types) matrix per disease, NA where padj >= 0.05
  build_mat <- function(dis) {
    sub <- gsea_norm %>%
      filter(disease == dis, pathway_clean %in% pathways, cell_type %in% ct_order) %>%
      mutate(NES = ifelse(!is.na(padj) & padj < 0.05, NES, NA_real_))
    mat <- matrix(NA_real_, nrow = length(pathways), ncol = length(ct_order),
                  dimnames = list(pathways, ct_order))
    for (i in seq_len(nrow(sub))) {
      mat[sub$pathway_clean[i], sub$cell_type[i]] <- sub$NES[i]
    }
    mat
  }
  mat_sle <- build_mat("SLE")
  mat_sjs <- build_mat("SjS")

  # Row order: preserve the order given in `pathways` (caller orders by
  # total absolute NES, matching panel a's bubble-plot x-axis).
  mat_sle <- mat_sle[pathways, , drop = FALSE]
  mat_sjs <- mat_sjs[pathways, , drop = FALSE]

  rownames(mat_sle) <- paste0("SLE__", rownames(mat_sle))
  rownames(mat_sjs) <- paste0("SjS__", rownames(mat_sjs))
  mat <- rbind(mat_sle, mat_sjs)
  disease_split <- factor(c(rep("SLE", nrow(mat_sle)), rep("SjS", nrow(mat_sjs))),
                          levels = c("SLE", "SjS"))
  row_labels <- gsub("^(SLE|SjS)__", "", rownames(mat))

  finite_vals <- mat[is.finite(mat)]
  max_abs <- if (length(finite_vals) > 0) min(max(abs(finite_vals)), 4.0) else 4.0
  col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs),
                                   c("#2166AC", "white", "#B2182B"))

  ha_right <- rowAnnotation(
    Disease = anno_block(
      gp = gpar(fill = disease_colors[c("SLE", "SjS")]),
      labels = c("SLE", "SjS"),
      labels_gp = gpar(col = "white", fontsize = 8, fontface = "bold")
    ),
    show_annotation_name = FALSE, show_legend = FALSE
  )

  ht <- Heatmap(mat, name = "NES", col = col_fun,
                na_col = "white",
                rect_gp = gpar(col = "white", lwd = 0.3),
                border = TRUE,
                right_annotation = ha_right,
                row_split = disease_split,
                column_title = "Per-cell-type pathway enrichment",
                column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                row_title = NULL,
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, show_column_names = TRUE,
                row_names_side = "left", column_names_side = "bottom",
                column_names_rot = 45,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 6),
                row_labels = row_labels,
                row_gap = unit(3, "mm"),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                            labels_gp = gpar(fontsize = 6),
                                            legend_height = unit(2.5, "cm")))

  grid.grabExpr(draw(ht, heatmap_legend_side = "right",
                     padding = unit(c(4, 6, 2, 2), "mm")))
}

build_fig3 <- function() {
  log_msg("--- M1 Figure 3: Pathway-level Shared vs Distinct ---")

  panel_a <- NULL
  heatmap_grob <- NULL

  gsea_file <- file.path(DIFFERENTIAL_EXPRESSION, "gsea_all_results.tsv")
  if (file.exists(gsea_file)) {
    tryCatch({
      # Suppress Rplots.pdf creation in non-interactive sessions
      pdf(NULL)
      on.exit(dev.off(), add = TRUE)

      gsea_raw <- read.delim(gsea_file, stringsAsFactors = FALSE)

      # Normalise: harmonise disease labels, strip variable "annotated" suffixes
      gsea_norm <- gsea_raw %>%
        mutate(disease = case_when(
          disease %in% c("DE_SLE", "SLE") ~ "SLE",
          disease %in% c("DE_SJS", "SjS") ~ "SjS",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(disease), grepl("^HALLMARK_", pathway)) %>%
        mutate(cell_type = trimws(gsub("(\\sannotated)+$", "", cell_type)),
               pathway_clean = vapply(pathway, format_pathway_label,
                                       character(1)))

      # Exclude rare cell types failing MIN_DONORS=5 / MIN_CELLS=50 thresholds.
      # After the DE rerun, the GSEA TSV will not contain these and this filter
      # becomes a no-op; the explicit list below is a fallback for preview.
      # This read was pointed at de_skipped_celltypes.tsv, which no script has ever
      # written - 06 writes de_skipped_summary.tsv. The branch therefore always fell
      # through to the hardcoded list below, so the exclusions applied here came from
      # a literal rather than from the floors that DE actually enforced. It happened to
      # be harmless only because 06 restricts DE to COMMON_CELL_TYPES before the floors
      # run, so those four types never reach GSEA in the first place. It would stop
      # being harmless the moment COMMON_CELL_TYPES changes: the literal would keep
      # stripping cell types from Fig 3 that DE had included.
      skip_path <- file.path(DIFFERENTIAL_EXPRESSION, "de_skipped_summary.tsv")
      if (!file.exists(skip_path)) {
        stop("de_skipped_summary.tsv not found at ", skip_path, ".\n",
             "  It records which cell types the DE floors excluded, and Fig 3 must ",
             "apply the same exclusions as the DE it summarises.\n",
             "  06_differential_analysis.R always writes it, empty if nothing was ",
             "skipped, so its absence means 06 has not run against the current atlas.\n",
             "  Refusing to substitute a hardcoded list, because that silently ",
             "decouples the figure from the analysis.")
      }
      excluded_cts <- unique(read.delim(skip_path, stringsAsFactors = FALSE)$cell_type)
      gsea_norm <- gsea_norm %>% filter(!cell_type %in% excluded_cts)

      # Deduplicate (pathway, disease, cell_type): keep most-significant variant
      gsea_norm <- gsea_norm %>%
        group_by(pathway_clean, disease, cell_type) %>%
        arrange(padj, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup()

      # Mean NES across ALL cell types (no padj filter — honest across-cell-type mean)
      pw_means <- gsea_norm %>%
        group_by(pathway_clean, disease) %>%
        summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop")

      # Number of significant cell types per pathway per disease (padj < 0.05)
      sig_counts <- gsea_norm %>%
        filter(!is.na(padj), padj < 0.05) %>%
        group_by(pathway_clean, disease) %>%
        summarise(n_sig = n(), .groups = "drop")

      pw <- pw_means %>%
        pivot_wider(names_from = disease, values_from = mean_NES,
                    names_prefix = "mean_NES_") %>%
        left_join(sig_counts %>% pivot_wider(names_from = disease,
                                              values_from = n_sig,
                                              names_prefix = "n_sig_"),
                  by = "pathway_clean") %>%
        mutate(n_sig_SLE = ifelse(is.na(n_sig_SLE), 0, n_sig_SLE),
               n_sig_SjS = ifelse(is.na(n_sig_SjS), 0, n_sig_SjS),
               mean_NES_SLE = ifelse(is.na(mean_NES_SLE), 0, mean_NES_SLE),
               mean_NES_SjS = ifelse(is.na(mean_NES_SjS), 0, mean_NES_SjS),
               total_abs_NES = abs(mean_NES_SLE) + abs(mean_NES_SjS))

      # Inclusion: significant in at least one cell type in at least one disease
      pw_keep <- pw %>%
        filter(n_sig_SLE > 0 | n_sig_SjS > 0) %>%
        arrange(desc(total_abs_NES))

      pathway_levels <- pw_keep$pathway_clean

      plot_long <- pw_keep %>%
        select(pathway_clean, total_abs_NES,
               mean_SLE = mean_NES_SLE, n_SLE = n_sig_SLE,
               mean_SjS = mean_NES_SjS, n_SjS = n_sig_SjS) %>%
        tidyr::pivot_longer(
          cols = -c(pathway_clean, total_abs_NES),
          names_to = c(".value", "disease"),
          names_sep = "_"
        ) %>%
        rename(mean_NES = mean, n_sig = n) %>%
        mutate(disease = factor(disease, levels = c("SLE", "SjS")),
               is_sig = n_sig > 0,
               fill_color = ifelse(is_sig,
                                   ifelse(disease == "SLE", "#D55E00", "#0072B2"),
                                   "white"),
               pathway_clean = factor(pathway_clean, levels = pathway_levels))

      panel_a <- ggplot(plot_long, aes(x = pathway_clean, y = mean_NES,
                                        group = disease)) +
        geom_hline(yintercept = 0, colour = "gray60",
                   linetype = "dashed", linewidth = 0.3) +
        geom_point(aes(size = n_sig, colour = disease, fill = fill_color),
                   shape = 21, stroke = 0.6, alpha = 0.6) +
        scale_colour_manual(values = c("SLE" = "#D55E00", "SjS" = "#0072B2"),
                            name = NULL) +
        scale_fill_identity() +
        scale_size_continuous(range = c(0.5, 3),
                              limits = c(0, NA),
                              breaks = c(0, 5, 15, 25),
                              name = "Sig. cell types") +
        labs(x = NULL, y = "Mean NES",
             title = "Hallmark pathway enrichment") +
        theme_nature(base_size = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
              plot.title = element_text(size = 9, face = "bold"),
              plot.margin = margin(t = 4, r = 4, b = 4, l = 12, unit = "mm"),
              legend.position = "right",
              legend.box = "vertical",
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 7),
              legend.key.size = unit(0.3, "cm"))

      # --- Panel B (combined heatmap): Hallmark pathways significant in at
      # least one cell type / disease, ordered by total absolute NES
      # (|mean_NES_SLE| + |mean_NES_SjS|) descending — identical to panel a's
      # x-axis. SLE on top, SjS on bottom; same row and column order in both.
      # Non-significant cells (padj >= 0.05) shown white. Pathways universally
      # non-significant in both diseases (4 of 50: bile acid metabolism, notch,
      # hedgehog, spermatogenesis) are dropped to match panel a's inclusion.
      all_hallmarks_ordered <- pathway_levels

      # Cell-type order: shared with Fig 2 panel A (DEG count descending,
      # SLE + SjS combined at padj < 0.05 & |log2FC| >= 1).
      ct_order_fig2 <- NULL
      sig_file_fig3 <- file.path(DIFFERENTIAL_EXPRESSION, "signatures",
                                 "disease_signatures_all.tsv")
      if (file.exists(sig_file_fig3)) {
        ct_order_fig2 <- read.delim(sig_file_fig3, stringsAsFactors = FALSE) %>%
          mutate(disease = case_when(
            disease %in% c("DE_SLE", "SLE") ~ "SLE",
            disease %in% c("DE_SJS", "SjS") ~ "SjS",
            TRUE ~ NA_character_
          )) %>%
          filter(disease %in% c("SjS", "SLE"), !is.na(padj),
                 padj < 0.05, abs(log2FoldChange) >= 1.0) %>%
          count(cell_type, name = "n_total") %>%
          arrange(desc(n_total)) %>%
          pull(cell_type)
      }

      heatmap_grob <- build_combined_pathway_heatmap(
        gsea_norm, all_hallmarks_ordered, excluded_cts,
        ct_order = ct_order_fig2)

      # --- Panel C: per-cell-type pathway sharing categories + Jaccard ---
      # Mirrors Fig 2 panel b structure (stacked shared/specific top + Jaccard
      # bottom), but applied to Hallmark pathways instead of DEGs.
      ct_levels_c <- c(ct_order_fig2,
                       setdiff(unique(gsea_norm$cell_type), ct_order_fig2))
      cats <- gsea_norm %>%
        filter(!is.na(padj)) %>%
        mutate(sig = padj < 0.05,
               sign = ifelse(NES > 0, "+", "-")) %>%
        group_by(cell_type, pathway_clean) %>%
        summarise(
          sle_sig  = any(disease == "SLE" & sig),
          sjs_sig  = any(disease == "SjS" & sig),
          sle_sign = ifelse(any(disease == "SLE"), sign[disease == "SLE"][1], NA_character_),
          sjs_sign = ifelse(any(disease == "SjS"), sign[disease == "SjS"][1], NA_character_),
          .groups = "drop")

      ct_counts <- cats %>%
        group_by(cell_type) %>%
        summarise(
          `SLE-specific` = sum(sle_sig & !sjs_sig),
          `SjS-specific` = sum(sjs_sig & !sle_sig),
          Shared         = sum(sle_sig & sjs_sig & sle_sign == sjs_sign),
          Discordant     = sum(sle_sig & sjs_sig & sle_sign != sjs_sign),
          n_union        = sum(sle_sig | sjs_sig),
          .groups = "drop") %>%
        mutate(jaccard = ifelse(n_union > 0, Shared / n_union, 0),
               cell_type = factor(cell_type, levels = ct_levels_c)) %>%
        arrange(cell_type)

      cat_long <- ct_counts %>%
        select(cell_type, `SLE-specific`, Shared, Discordant, `SjS-specific`) %>%
        tidyr::pivot_longer(cols = -cell_type,
                            names_to = "category", values_to = "count") %>%
        mutate(category = factor(category,
                                  levels = c("SLE-specific", "Shared",
                                             "Discordant", "SjS-specific")))

      cat_colors_pw <- c("SLE-specific" = "#D55E00",
                         "Shared"       = "#CC79A7",   # Okabe-Ito reddish purple (palette family); green is HC elsewhere
                         "Discordant"   = "#BBBBBB",   # lighter grey: moves away from Shared #CC79A7 under deuteranopia
                         "SjS-specific" = "#0072B2")

      ct_tot <- aggregate(count ~ cell_type, cat_long, sum)
      panel_c_top <- ggplot(cat_long, aes(x = cell_type, y = count)) +
        geom_col(aes(fill = category), width = 0.9, color = NA,
                 position = position_stack(reverse = TRUE)) +
        stack_outline(ct_tot, "cell_type", "count", width = 0.9, size = 0.2) +
        scale_fill_manual(values = cat_colors_pw, name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(y = "Number of pathways", x = NULL,
             title = "Shared vs Specific Pathways") +
        theme_nature(base_size = 8) +
        theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin  = margin(5, 5, 0, 30),
          legend.position = "top",
          legend.justification = "right",
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.25, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major.x = element_blank()
        )

      panel_c_bot <- ggplot(ct_counts, aes(x = cell_type, y = jaccard)) +
        ggchicklet::geom_chicklet(fill = "#009E73", color = "black", size = 0.3,
                                   radius = BAR_RADIUS) +
        scale_y_continuous(limits = c(0, max(ct_counts$jaccard) * 1.05),
                           expand = c(0, 0)) +
        labs(y = "Jaccard similarity", x = NULL) +
        theme_nature(base_size = 8) +
        theme(
          axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
          plot.margin = margin(5, 5, 5, 30),
          panel.grid.major.x = element_blank()
        )

      panel_c <- cowplot::plot_grid(panel_c_top, panel_c_bot, ncol = 1,
                                     rel_heights = c(0.7, 0.6),
                                     align = "v", axis = "lr")
    }, error = function(e) {
      log_msg("  WARNING: Pathway analysis failed: ", e$message)
    })
  }

  # Fallbacks
  if (is.null(panel_a)) panel_a <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No pathway data")
  if (is.null(heatmap_grob)) heatmap_grob <- grid::textGrob("Pathway heatmap\nnot available", gp = grid::gpar(col = "grey50"))
  if (!exists("panel_c") || is.null(panel_c)) panel_c <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Panel c not available")

  # Compose: 3 rows (a = bubble plot, b = combined SLE/SjS heatmap, c = pathway
  # sharing + Jaccard mirroring Fig 2 panel b structure)
  row1 <- cowplot::plot_grid(panel_a, labels = "a",
                              label_size = LABEL_SIZE, label_fontface = "bold")
  row2 <- cowplot::plot_grid(heatmap_grob, labels = "b",
                              label_size = LABEL_SIZE, label_fontface = "bold")
  row3 <- cowplot::plot_grid(panel_c, labels = "c",
                              label_size = LABEL_SIZE, label_fontface = "bold")
  fig <- cowplot::plot_grid(row1, row2, row3, ncol = 1,
                             rel_heights = c(0.16, 0.64, 0.20))

  outpath <- file.path(OUTDIR, "M1_Figure3_pathways.png")
  # 14.0 in: panel b holds 100 label rows (50 pathways x 2 diseases) at 5 pt;
  # 0.64 x 14.0 in gives ~2.0 mm per row inside the PDF (label height 1.76 mm).
  save_nature_fig(fig, outpath, width = FIG_WIDTH, height = 14.0)
}

# =============================================================================
# Figure 4: Donor-level Classification (4 panels, per-class SHAP)
#   a: Out-of-fold confusion matrix (elastic net, per-donor consensus)
#   b: Cell type ranking: one stacked bar per cell type, split by class
#      (HC / SLE / SjS); segment height = per-class mean |SHAP|. Ordered by the
#      class-averaged mean, which is one third of the stacked height.
#   c: Pathway ranking, stacked as in b
#   d: Cell type × pathway |SHAP| toward each class: three heatmaps, one per
#      class, sharing the row/column order of b and c and one sqrt colour scale
#
# Bars in b and c follow the shared stacked-bar convention (stack_outline():
# plain geom_col segments under one no-fill rounded outline per bar, 0.75 pt).
# No fold-SD error bars: per-class fold SDs are large by construction for the
# 14-donor SjS class and are not reported anywhere.
#
# Inputs: results/donor_classifiers_<method>/ (DENOISE_METHOD env var, mandatory)
# — cv_predictions.tsv; shap_<clf>/shap_importance_per_{celltype,pathway}.tsv
# (class-averaged, used for ordering only) and
# shap_importance_per_{feature,celltype,pathway}_per_class.tsv (per-class
# values). Produced upstream by scripts 08 (training) and 09 (SHAP).
# =============================================================================

build_fig4 <- function() {
  log_msg("--- M1 Figure 4: Donor-level Classification (per-class SHAP) ---")

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
  SRC_DIR <- switch(DENOISE_METHOD,
    scvi              = file.path(ATLAS_ROOT, "results", "donor_classifiers_scvi"),
    limma             = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma"),
    limma_modulescore = file.path(ATLAS_ROOT, "results", "donor_classifiers_limma_modulescore"),
    file.path(ATLAS_ROOT, "results", "donor_classifiers"))
  HEADLINE_CLF <- switch(DENOISE_METHOD,
    limma_modulescore = "glmnet",
    "xgb")
  log_msg("  SRC_DIR: ", SRC_DIR)
  log_msg("  Headline classifier: ", HEADLINE_CLF)

  # ---- Load data ----
  cv_preds <- read.delim(file.path(SRC_DIR, "cv_predictions.tsv"),
                         stringsAsFactors = FALSE)
  cv_preds <- cv_preds[cv_preds$classifier == HEADLINE_CLF, ]
  SHAP_DIR <- file.path(SRC_DIR, paste0("shap_", HEADLINE_CLF))
  # Class-averaged rankings: used ONLY to fix the cell-type / pathway ordering.
  per_ct   <- read.delim(file.path(SHAP_DIR, "shap_importance_per_celltype.tsv"),
                         stringsAsFactors = FALSE)
  per_pw   <- read.delim(file.path(SHAP_DIR, "shap_importance_per_pathway.tsv"),
                         stringsAsFactors = FALSE)
  # Per-class values (script 09): one row per unit × class.
  per_ct_class   <- read.delim(file.path(SHAP_DIR, "shap_importance_per_celltype_per_class.tsv"),
                               stringsAsFactors = FALSE)
  per_pw_class   <- read.delim(file.path(SHAP_DIR, "shap_importance_per_pathway_per_class.tsv"),
                               stringsAsFactors = FALSE)
  per_feat_class <- read.delim(file.path(SHAP_DIR, "shap_importance_per_feature_per_class.tsv"),
                               stringsAsFactors = FALSE)
  stopifnot(setequal(unique(per_feat_class$class), names(COND_LABEL)))

  # Sanity check: the class mean of the per-class values must equal the
  # class-averaged values that define the ordering.
  chk <- aggregate(mean_mean_abs_shap ~ celltype, per_ct_class, mean)
  chk <- merge(chk, per_ct, by = "celltype")
  stopifnot(max(abs(chk$mean_mean_abs_shap.x - chk$mean_mean_abs_shap.y)) < 1e-10)
  chk <- aggregate(mean_mean_abs_shap ~ pathway, per_pw_class, mean)
  chk <- merge(chk, per_pw, by = "pathway")
  stopifnot(max(abs(chk$mean_mean_abs_shap.x - chk$mean_mean_abs_shap.y)) < 1e-10)

  # ---- Panel A: confusion matrix (per-donor consensus) ----
  class_order <- c("healthy", "sle", "sjs")
  donor_preds <- cv_preds %>%
    dplyr::group_by(donor_id, true) %>%
    dplyr::summarise(prob_healthy = mean(prob_healthy),
                     prob_sle     = mean(prob_sle),
                     prob_sjs     = mean(prob_sjs),
                     .groups = "drop")
  prob_mat <- as.matrix(donor_preds[, c("prob_healthy", "prob_sle", "prob_sjs")])
  donor_preds$pred <- class_order[max.col(prob_mat, ties.method = "first")]
  stopifnot(nrow(donor_preds) == 291L)

  cm <- table(true = factor(donor_preds$true, levels = class_order),
              pred = factor(donor_preds$pred, levels = class_order))
  cm_df <- as.data.frame(cm)
  cm_df$row_total <- as.numeric(rowSums(cm)[as.character(cm_df$true)])
  cm_df$row_pct   <- cm_df$Freq / cm_df$row_total
  cm_df$true_lab  <- factor(COND_LABEL[as.character(cm_df$true)],
                            levels = unname(COND_LABEL[class_order]))
  cm_df$pred_lab  <- factor(COND_LABEL[as.character(cm_df$pred)],
                            levels = unname(COND_LABEL[class_order]))

  panel_a <- ggplot(cm_df, aes(x = pred_lab, y = true_lab, fill = row_pct)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Freq), size = 3.0, fontface = "bold", vjust = -0.5) +
    geom_text(aes(label = sprintf("%.1f%%", row_pct * 100)),
              size = 2.3, vjust = 1.5, color = "grey25") +
    scale_fill_gradient(low = "white", high = "#B71C1C",
                        limits = c(0, 1), guide = "none") +
    coord_fixed() +
    labs(x = "Predicted", y = "True") +
    theme_nature() +
    theme(axis.text.x  = element_text(face = "bold"),
          axis.text.y  = element_text(face = "bold"),
          panel.grid   = element_blank(),
          axis.ticks   = element_blank(),
          plot.margin  = margin(8, 12, 8, 12))

  # ---- Panel B: cell type ranking (stacked by class; class-averaged ordering) ----
  # Stack total = sum over the three classes = 3 × the class-averaged mean |SHAP|,
  # so ordering by the class-averaged value orders by stack height.
  per_ct <- per_ct[order(-per_ct$mean_mean_abs_shap), ]
  ct_clean_order <- clean_ct(per_ct$celltype)
  per_ct_class$celltype_clean <- factor(clean_ct(per_ct_class$celltype),
                                        levels = ct_clean_order)
  per_ct_class$class_lab <- factor(COND_LABEL[per_ct_class$class],
                                   levels = unname(COND_LABEL))
  per_ct_tot <- aggregate(mean_mean_abs_shap ~ celltype_clean, per_ct_class, sum)
  panel_b <- ggplot(per_ct_class, aes(x = celltype_clean, y = mean_mean_abs_shap)) +
    geom_col(aes(fill = class_lab), width = 0.75, color = NA,
             position = position_stack(reverse = TRUE)) +
    stack_outline(per_ct_tot, "celltype_clean", "mean_mean_abs_shap", width = 0.75) +
    scale_fill_manual(values = CLASS_FILL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Mean |SHAP|") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
          axis.text.x = element_text(size = BASE_SIZE - 2,
                                     angle = 45, hjust = 1, vjust = 1),
          legend.position = "inside",
          legend.position.inside = c(0.92, 0.85),
          legend.key.size = unit(3, "mm"),
          legend.background = element_blank())

  # ---- Panel C: pathway ranking (stacked by class; class-averaged ordering) ----
  per_pw <- per_pw[order(-per_pw$mean_mean_abs_shap), ]
  path_clean_order <- vapply(per_pw$pathway, clean_path, character(1))
  per_pw_class$pathway_clean <- factor(vapply(per_pw_class$pathway, clean_path, character(1)),
                                       levels = path_clean_order)
  per_pw_class$class_lab <- factor(COND_LABEL[per_pw_class$class],
                                   levels = unname(COND_LABEL))
  per_pw_tot <- aggregate(mean_mean_abs_shap ~ pathway_clean, per_pw_class, sum)
  panel_c <- ggplot(per_pw_class, aes(x = pathway_clean, y = mean_mean_abs_shap)) +
    geom_col(aes(fill = class_lab), width = 0.75, color = NA,
             position = position_stack(reverse = TRUE)) +
    stack_outline(per_pw_tot, "pathway_clean", "mean_mean_abs_shap", width = 0.75) +
    scale_fill_manual(values = CLASS_FILL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Mean |SHAP|") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
          axis.text.x = element_text(size = BASE_SIZE - 2.5,
                                     angle = 45, hjust = 1, vjust = 1),
          plot.margin = margin(8, 16, 8, 16),
          legend.position = "inside",
          legend.position.inside = c(0.95, 0.85),
          legend.key.size = unit(3, "mm"),
          legend.background = element_blank())

  # ---- Panel D: heatmap, one per class, shared order and colour scale ----
  heat_df <- per_feat_class
  heat_df$celltype_clean <- factor(clean_ct(heat_df$celltype), levels = ct_clean_order)
  heat_df$pathway_clean  <- factor(vapply(heat_df$pathway, clean_path, character(1)),
                                   levels = rev(path_clean_order))
  heat_df$class_lab      <- factor(COND_LABEL[heat_df$class], levels = unname(COND_LABEL))

  panel_d <- ggplot(heat_df, aes(x = celltype_clean, y = pathway_clean,
                                 fill = mean_abs_shap)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "rocket", direction = -1,
                         name = "Mean\n|SHAP|", trans = "sqrt",
                         guide = guide_colorbar(barwidth = unit(1.5, "mm"),
                                                barheight = unit(30, "mm"),
                                                ticks.colour = "white",
                                                frame.colour = NA)) +
    facet_wrap(~ class_lab, ncol = 3) +
    labs(x = NULL, y = NULL) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = BASE_SIZE - 3),
          axis.text.y = element_text(size = BASE_SIZE - 3),
          panel.grid  = element_blank(),
          axis.ticks  = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          panel.spacing.x = unit(1, "mm"),
          strip.text = element_text(size = BASE_SIZE, face = "bold"),
          plot.margin = margin(2, 0, 2, 0),
          legend.position = "right",
          legend.title = element_text(size = BASE_SIZE - 2.5),
          legend.text  = element_text(size = BASE_SIZE - 3),
          legend.margin = margin(2, 2, 2, 2),
          legend.box.spacing = unit(2, "mm"))

  # ---- Compose ----
  row1 <- cowplot::plot_grid(panel_a, NULL, panel_b, ncol = 3,
                             rel_widths = c(1, 0.08, 2),
                             labels = c("a", "", "b"),
                             label_size = LABEL_SIZE,
                             label_fontface = "bold", label_y = 1.0)
  panel_c_labeled <- cowplot::plot_grid(panel_c, labels = "c",
                                        label_size = LABEL_SIZE,
                                        label_fontface = "bold", label_y = 1.0)
  panel_d_labeled <- cowplot::plot_grid(panel_d, labels = "d",
                                        label_size = LABEL_SIZE,
                                        label_fontface = "bold", label_y = 1.0)
  fig4 <- cowplot::plot_grid(row1, panel_c_labeled, panel_d_labeled,
                             ncol = 1, rel_heights = c(60, 60, 130))

  HEIGHT_MM <- 60 + 60 + 130 + 8
  HEIGHT_IN <- HEIGHT_MM / 25.4
  outpath <- file.path(OUTDIR, "M1_Figure4_classification.pdf")
  save_nature_fig(fig4, outpath, width = FIG_WIDTH, height = HEIGHT_IN)
  # PNG for review
  ggsave(sub("\\.pdf$", ".png", outpath), fig4,
         width = FIG_WIDTH, height = HEIGHT_IN, dpi = 300, bg = "white")
}

# =============================================================================
# Copy PDFs to manuscript/figures/
# =============================================================================

FIG_MAP <- list(
  c("M1_Figure2_disease_programs.pdf",      "fig2_disease_programs.pdf"),
  c("M1_Figure3_pathways.pdf",              "fig3_pathways.pdf"),
  c("M1_Figure4_classification.pdf",        "fig4_classification.pdf")
)

copy_to_manuscript <- function(fig_num = NULL) {
  MANUSCRIPT_FIGDIR <- file.path(ATLAS_ROOT, "manuscript", "figures")
  dir.create(MANUSCRIPT_FIGDIR, recursive = TRUE, showWarnings = FALSE)
  # When fig_num is set (i.e. called from --figure N), restrict the copy to
  # the FIG_MAP entry whose source filename references that figure number,
  # so a single-figure rebuild can never propagate a stale sibling PDF.
  pattern <- if (!is.null(fig_num)) paste0("^M1_Figure", fig_num, "_") else NULL
  for (fm in FIG_MAP) {
    if (!is.null(pattern) && !grepl(pattern, fm[1])) next
    src <- file.path(OUTDIR, fm[1])
    dst <- file.path(MANUSCRIPT_FIGDIR, fm[2])
    if (file.exists(src)) {
      file.copy(src, dst, overwrite = TRUE)
      log_msg("  Copied: ", fm[1], " -> manuscript/figures/", fm[2])
    }
  }
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  log_msg("=" %>% rep(60) %>% paste(collapse = ""))
  log_msg("  Composing Manuscript 1 Figures (Nature Comms: 180mm width, 450 DPI, Helvetica)")
  log_msg("=" %>% rep(60) %>% paste(collapse = ""))
  log_msg("  Output: ", OUTDIR)

  build_fig1()
  build_fig2()
  build_fig3()
  build_fig4()

  copy_to_manuscript()

  log_msg("")
  log_msg("=" %>% rep(60) %>% paste(collapse = ""))
  log_msg("  Manuscript 1 figures saved (180mm width, PNG + PDF)")
  log_msg("  Main: figs 1, 2, 3, 4")
  log_msg("  ", OUTDIR)
  log_msg("=" %>% rep(60) %>% paste(collapse = ""))
}

# Run main with optional figure selection
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2 && args[1] == "--figure") {
    fig_num <- as.integer(args[2])
    log_msg("=" %>% rep(60) %>% paste(collapse = ""))
    log_msg("  Running single figure: ", fig_num)
    log_msg("=" %>% rep(60) %>% paste(collapse = ""))
    log_msg("  Output: ", OUTDIR)
    switch(as.character(fig_num),
      `1` = build_fig1(),
      `2` = build_fig2(),
      `3` = build_fig3(),
      `4` = build_fig4(),
      stop("Unknown figure number: ", fig_num, ". Valid: 1, 2, 3, 4")
    )
    copy_to_manuscript(fig_num = fig_num)
  } else {
    main()
  }
}
