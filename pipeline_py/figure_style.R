# Figure style shared by 08_figures.R (manuscript figures) and 10_supplementary.R (supplementary figures and tables):
# constants, theme, panel titles, bar and strip conventions, colours, dataset labels and the cell-type order. Sourced after the caller sets PIPELINE_DIR, DRAFT, DRAFT_TEXT, OUT and
# log_msg.

# =============================================================================
# Style (Springer Nature artwork guide: text 5-7 pt sans serif, 85 or 175 mm wide)
# =============================================================================
FIG_WIDTH      <- 175 / 25.4   # two-column width 175 mm
FIG_MAX_HEIGHT <- 240 / 25.4   # tallest figure that fits a page
MIN_PT         <- 7
BLOCK_MM       <- 7            # depth of the coloured group blocks (heatmap blocks, facet strips)
FRAME_LWD      <- 1            # frame of the heatmaps and their blocks (grid lwd; 1 = 0.75 pt)
BASE_SIZE      <- 7
LABEL_SIZE     <- 7            # text at most 7 pt
TITLE_PT       <- 7            # text at most 7 pt
FONT_FAMILY    <- "Helvetica"
TEXT_COL       <- "#000000"    # every text element: set in theme_nature() and the heatmap gpars

theme_nature <- function(base_size = BASE_SIZE, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # theme_minimal() draws axis text grey30 and strip text grey10; all text is black
      text = element_text(family = base_family, colour = TEXT_COL),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, colour = TEXT_COL),
      plot.subtitle = element_text(size = base_size, colour = TEXT_COL),
      axis.title = element_text(size = base_size, colour = TEXT_COL),
      axis.text = element_text(size = base_size, colour = TEXT_COL),
      legend.title = element_text(size = base_size, colour = TEXT_COL),
      legend.text = element_text(size = base_size, colour = TEXT_COL),
      strip.text = element_text(colour = TEXT_COL),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Panel title, the same in every figure: TITLE_PT bold, centred on the panel's plotting area, 4 pt above it (above any
# legend placed at the top). Added last, so no later theme can restyle it.
panel_title <- function(text) {
  list(labs(title = text),
       theme(plot.title = element_text(size = TITLE_PT, face = "bold", hjust = 0.5, margin = margin(b = 4)),
             plot.title.position = "panel"))
}

# Bar convention: BAR_RADIUS is the corner radius of every bar. Single bars are plain chicklets. Stacked bars use
# stack_outline(): plain geom_col segments with no stroke under one no-fill rounded outline per bar, so a thin segment
# can never invert; bars below STACK_MIN_FRAC of the tallest bar get a square outline.
BAR_RADIUS     <- grid::unit(0.75, "pt")
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

# Colour each facet strip by its label, drawn like the group blocks of the pathway heatmap: filled, FRAME_LWD black
# frame, white bold label, BLOCK_MM deep, 1 mm from the panel. ggplot fills every strip alike, so the fills and the
# strip width are set on the gtable.
class_strips <- function(p, fills, STRIP_MM = BLOCK_MM, GAP_MM = 1) {
  g <- ggplotGrob(p)
  labels_of <- function(x) c(if (is.character(x$label)) x$label,
                             unlist(lapply(x$grobs, labels_of)),
                             unlist(lapply(x$children, labels_of)))
  refill <- function(x, fill) {
    if (inherits(x, "rect")) x$gp$fill <- fill
    if (!is.null(x$grobs)) x$grobs <- lapply(x$grobs, refill, fill = fill)
    if (length(x$children)) {
      for (nm in names(x$children)) x$children[[nm]] <- refill(x$children[[nm]], fill)
    }
    x
  }
  idx <- grep("^strip-r", g$layout$name)
  for (i in idx) {
    g$grobs[[i]] <- refill(g$grobs[[i]], fills[[labels_of(g$grobs[[i]])[1]]])
    g$grobs[[i]]$widths <- unit(STRIP_MM, "mm")
    g$grobs[[i]]$layout$clip <- "off"
  }
  g$layout$clip[idx] <- "off"
  col <- unique(g$layout$l[idx])
  stopifnot(length(col) == 1)
  g$widths[col] <- unit(STRIP_MM, "mm")
  gtable::gtable_add_cols(g, unit(GAP_MM, "mm"), pos = col - 1)
}

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

draft_mark <- function(p) {
  if (!DRAFT) return(p)
  ggdraw(p) + draw_label(DRAFT_TEXT, x = 0.995, y = 0.998, hjust = 1, vjust = 1, colour = "red",
                         fontface = "bold", size = MIN_PT, fontfamily = FONT_FAMILY)
}
# PDF device of every figure (08_figures.R, 10_supplementary.R): Cairo::CairoPDF, not grDevices::cairo_pdf, which (pango
# 1.52 in the image) rounds glyph positions, closing spaces and splitting words at 5-7 pt; same embedded fonts
pdf_device <- function(filename, width, height, ...) Cairo::CairoPDF(filename, width, height, family = FONT_FAMILY, bg = "white")
save_fig <- function(p, name, height_mm) {
  stopifnot(height_mm <= FIG_MAX_HEIGHT * 25.4)
  f <- file.path(OUT, name)
  ggsave(f, draft_mark(p), width = FIG_WIDTH, height = height_mm / 25.4, device = pdf_device, bg = "white")
  log_msg("saved", f)
}

# =============================================================================
# Groups, datasets, cell types
# =============================================================================
CONDITION_COLORS <- c("HC" = "#009E73", "SLE" = "#D55E00", "SjD" = "#0072B2")
DISEASE_COLORS <- c("SLE" = "#D55E00", "SjD" = "#0072B2")
cond_label <- function(x) ifelse(x == "control", "HC", x)
datasets <- fread(file.path(PIPELINE_DIR, "datasets.tsv"), skip = "dataset\t")
# Datasets are named by their accession (GEO or GSA-Human), in the figures and in Table 1
DATASET_LABEL <- c("Perez" = "GSE174188", "GSE135779_adult" = "GSE135779", "GSE253568" = "GSE253568",
                   "GSE316482" = "GSE316482", "HRA003613" = "HRA003613")
stopifnot(setequal(names(DATASET_LABEL), datasets$dataset))
DATASET_LABEL <- DATASET_LABEL[datasets$dataset]
DATASET_COLORS <- setNames(scales::hue_pal()(length(DATASET_LABEL)), DATASET_LABEL)
# The annotation label "MAIT and gamma-delta T" (the key of every step output) is shown as "MAIT and γδ T", written
# with escapes so the scripts do not depend on the locale; read_ct() is fread() with the cell-type columns relabelled
ct_display <- function(x) replace(x, x == "MAIT and gamma-delta T", "MAIT and \u03b3\u03b4 T")
read_ct <- function(...) {
  x <- fread(...)
  for (c in intersect(c("cell_type", "Cell Type"), names(x))) set(x, j = c, value = ct_display(x[[c]]))
  x
}
# Cell types in lineage order (legend of Fig 2c); colours are shades of their lineage colour
LINEAGE <- c("Naive CD4 T" = "CD4 T", "Memory CD4 T" = "CD4 T", "Regulatory T" = "CD4 T",
             "Naive CD8 T" = "CD8 T", "CD8 effector memory T" = "CD8 T", "CD8 effector T (TEMRA)" = "CD8 T",
             "MAIT and \u03b3\u03b4 T" = "Unconventional T", "NK" = "NK cells",
             "Cycling T&NK" = "NK cells",  # 2026-11-03: the 16th type of the intersect run, next to NK
             "Naive B" = "B cells", "Memory B" = "B cells",
             "Classical monocytes" = "Monocytes", "Non-classical monocytes" = "Monocytes",
             "cDC2" = "DCs", "pDC" = "DCs", "HSPC" = "Progenitor")
CT_ORDER <- names(LINEAGE)
LINEAGE_ORDER <- c("CD4 T", "CD8 T", "Unconventional T", "B cells", "NK cells", "Monocytes", "DCs", "Progenitor")
LINEAGE_COLORS <- c("CD4 T" = "#E41A1C", "CD8 T" = "#FF7F00", "Unconventional T" = "#FDBF6F", "B cells" = "#377EB8",
                    "NK cells" = "#4DAF4A", "Monocytes" = "#A65628", "DCs" = "#F781BF", "Progenitor" = "#999999")
CT_COLORS <- c("Naive CD4 T" = "#FB9A99", "Memory CD4 T" = "#E41A1C", "Regulatory T" = "#99000D",
               "Naive CD8 T" = "#FDD0A2", "CD8 effector memory T" = "#FF7F00", "CD8 effector T (TEMRA)" = "#A63603",
               "MAIT and \u03b3\u03b4 T" = "#E6C229", "NK" = "#4DAF4A", "Cycling T&NK" = "#A1D99B",
               "Naive B" = "#9ECAE1", "Memory B" = "#08519C",
               "Classical monocytes" = "#A65628", "Non-classical monocytes" = "#D9A066",
               "cDC2" = "#F781BF", "pDC" = "#984EA3", "HSPC" = "#999999")

