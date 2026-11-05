#!/usr/bin/env Rscript
# Every figure and table of the paper, main and supplementary, from the data 07_figure_data.py prepares; it reads
# nothing else. Not a pipeline step and not a notebook port. Style: theme_nature(), panel_title(), the bar conventions
# (BAR_RADIUS, stack_outline()), the coloured facet strips (class_strips()) and the heatmap layout, sized for the Springer
# Nature artwork guide (175 mm wide, text 5-7 pt). Fig 1 is TikZ (fig1_pipeline.tex).
# Usage: run by make_figures.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration; FIG_DRAFT=1 adds a
#        DRAFT mark)
# Input: $FIG_DATA, default $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/07_figures/data/ (07_figure_data.py)
# Output, overwritten on each run:
#   manuscript/figures/ ($FIG_OUT)        fig2_atlas.pdf, fig3_de.pdf, fig4_pathways.pdf
#   manuscript/tables/ ($TABLE_OUT)       table1.xlsx
#   manuscript/supplementary/ ($SUPP_OUT) one Excel workbook per table, flat values with a title row:
#     supp_table_qc.xlsx           cells per dataset after each quality-control filter
#     supp_table_composition.xlsx  scCODA joint model, the three contrasts, with the references credible with
#     supp_table_de_sle.xlsx       PyDESeq2, SLE vs HC, every gene tested
#     supp_table_de_sjd.xlsx       PyDESeq2, SjD vs HC, every gene tested
#     supp_table_sle_vs_sjd.xlsx   SLE vs SjD indirect comparison per gene and cell type, with the Fig 3b classes
#     supp_table_gsea.xlsx         fgsea, Hallmark 2025.1, the three contrasts
#     supp_fig_atlas.pdf                UMAPs before integration and of the controls, annotation clusters, marker dot plot
#     supp_fig_composition_cohorts.pdf  per-cohort composition heatmap
#     supp_fig_de_cohorts.pdf           per-cohort DE and the fixed-effect meta-analysis vs the joint model
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot); library(grid)
  library(ComplexHeatmap); library(circlize); library(data.table)
})

ATLAS_RUN <- Sys.getenv("ATLAS_RUN", "02_integration")
ATLAS_ROOT <- Sys.getenv("ATLAS_ROOT", ".")  # repository root (run from it, or set ATLAS_ROOT)
RUN  <- file.path(Sys.getenv("PIPELINE_OUT", file.path(ATLAS_ROOT, "out")), "py", "atlas", ATLAS_RUN)
DATA <- Sys.getenv("FIG_DATA", file.path(RUN, "07_figures", "data"))
FIG_OUT <- Sys.getenv("FIG_OUT", file.path(ATLAS_ROOT, "manuscript", "figures"))
TABLE_OUT <- Sys.getenv("TABLE_OUT", file.path(ATLAS_ROOT, "manuscript", "tables"))
SUPP_OUT <- Sys.getenv("SUPP_OUT", file.path(ATLAS_ROOT, "manuscript", "supplementary"))
for (d in c(FIG_OUT, TABLE_OUT, SUPP_OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
DRAFT <- Sys.getenv("FIG_DRAFT", "0") == "1"  # FIG_DRAFT=1 adds a red DRAFT | run | date mark
DRAFT_TEXT <- paste0("DRAFT | ", ATLAS_RUN, " | ", format(Sys.Date()))
log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
pdf(NULL)  # null default device: plot measurements otherwise write Rplots.pdf into the working directory (the repo)

# =============================================================================
# Style (Springer Nature artwork guide)
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
# PDF device of every figure: Cairo::CairoPDF, not grDevices::cairo_pdf, which (pango
# 1.52 in the image) rounds glyph positions, closing spaces and splitting words at 5-7 pt; same embedded fonts
pdf_device <- function(filename, width, height, ...) Cairo::CairoPDF(filename, width, height, family = FONT_FAMILY, bg = "white")

save_fig <- function(p, dir, name, height_mm) {
  stopifnot(height_mm <= FIG_MAX_HEIGHT * 25.4)
  f <- file.path(dir, name)
  ggsave(f, draft_mark(p), width = FIG_WIDTH, height = height_mm / 25.4, device = pdf_device, bg = "white")
  log_msg("saved", f)
}
# Tables: one flat workbook each, title in row 1 (and a DRAFT mark in row 2 with FIG_DRAFT=1), header, data, notes one
# blank row below; flat values, no formulas
write_table <- function(df, path, sheet, title, notes = character(), int_cols = character(), dec_cols = character(),
                        widths = "auto") {
  r0 <- if (DRAFT) 3 else 2
  wb <- openxlsx2::wb_workbook()
  wb <- openxlsx2::wb_add_worksheet(wb, sheet)
  wb <- openxlsx2::wb_add_data(wb, x = title, start_row = 1, start_col = 1)
  wb <- openxlsx2::wb_add_font(wb, dims = "A1", bold = TRUE)
  if (DRAFT) {
    wb <- openxlsx2::wb_add_data(wb, x = DRAFT_TEXT, start_row = 2, start_col = 1)
    wb <- openxlsx2::wb_add_font(wb, dims = "A2", bold = TRUE, color = openxlsx2::wb_color(hex = "FFFF0000"))
  }
  wb <- openxlsx2::wb_add_data(wb, x = as.data.frame(df), start_row = r0, start_col = 1, na.strings = "")
  wb <- openxlsx2::wb_add_font(wb, dims = openxlsx2::wb_dims(rows = r0, cols = seq_len(ncol(df))), bold = TRUE)
  rows <- r0 + seq_len(nrow(df))
  for (c in match(int_cols, names(df)))
    wb <- openxlsx2::wb_add_numfmt(wb, dims = openxlsx2::wb_dims(rows = rows, cols = c), numfmt = "#,##0")
  for (c in match(dec_cols, names(df)))
    wb <- openxlsx2::wb_add_numfmt(wb, dims = openxlsx2::wb_dims(rows = rows, cols = c), numfmt = "0.000")
  for (i in seq_along(notes))
    wb <- openxlsx2::wb_add_data(wb, x = notes[i], start_row = max(rows) + 1 + i, start_col = 1)
  wb <- openxlsx2::wb_set_col_widths(wb, cols = seq_len(ncol(df)), widths = widths)
  openxlsx2::wb_save(wb, path, overwrite = TRUE)
  log_msg("saved", path, "|", nrow(df), "rows")
}
# Supplementary tables, numbered in order of first citation in main.tex. manuscript/supplementary/numbering.tsv holds the
# same numbers for manuscript/tools/build_supplementary.py, which checks them against main.tex: change both together.
SUPP_TABLE_NO <- c(supp_table_qc.xlsx = 1, supp_table_composition.xlsx = 2, supp_table_de_sle.xlsx = 3,
                   supp_table_de_sjd.xlsx = 4, supp_table_sle_vs_sjd.xlsx = 5, supp_table_gsea.xlsx = 6)
write_supp_table <- function(df, file, title, ...)
  write_table(df, file.path(SUPP_OUT, file), "Table", sprintf("Supplementary Table %d. %s", SUPP_TABLE_NO[[file]], title), ...)

# =============================================================================
# Groups, datasets, cell types
# =============================================================================
CONDITION_COLORS <- c("HC" = "#009E73", "SLE" = "#D55E00", "SjD" = "#0072B2")
DISEASE_COLORS <- c("SLE" = "#D55E00", "SjD" = "#0072B2")
cond_label <- function(x) ifelse(x == "control", "HC", x)
# The annotation label "MAIT and gamma-delta T" (the key of every step output) is shown as "MAIT and γδ T", written
# with escapes so the scripts do not depend on the locale; read_ct() is fread() with the cell-type columns relabelled
ct_display <- function(x) replace(x, x == "MAIT and gamma-delta T", "MAIT and γδ T")
read_ct <- function(...) {
  x <- fread(...)
  for (c in intersect(c("cell_type", "Cell Type"), names(x))) set(x, j = c, value = ct_display(x[[c]]))
  x
}
read_gz <- function(f) read_ct(cmd = paste("zcat", shQuote(file.path(DATA, f))))  # no R.utils in the image
# Datasets: those of Table 1 (07_figure_data.py writes them in the order of pipeline_py/datasets.tsv), named by their
# accession (GEO or GSA-Human) in the figures and tables
t1 <- read_ct(file.path(DATA, "table1.tsv"))
datasets <- t1[, .(dataset, disease)]
DATASET_LABEL <- c("Perez" = "GSE174188", "GSE135779_adult" = "GSE135779", "GSE253568" = "GSE253568",
                   "GSE316482" = "GSE316482", "HRA003613" = "HRA003613")
stopifnot(setequal(names(DATASET_LABEL), datasets$dataset))
DATASET_LABEL <- DATASET_LABEL[datasets$dataset]
DATASET_COLORS <- setNames(scales::hue_pal()(length(DATASET_LABEL)), DATASET_LABEL)
SLE_DS <- datasets$dataset[datasets$disease == "SLE"]; SJD_DS <- datasets$dataset[datasets$disease == "SjD"]
# Cell types in lineage order (legend of Fig 2c), the one cell-type order of every figure and table; colours are shades of
# their lineage colour
LINEAGE <- c("Naive CD4 T" = "CD4 T", "Memory CD4 T" = "CD4 T", "Regulatory T" = "CD4 T",
             "Naive CD8 T" = "CD8 T", "CD8 effector memory T" = "CD8 T", "CD8 effector T (TEMRA)" = "CD8 T",
             "MAIT and γδ T" = "Unconventional T", "NK" = "NK cells",
             "Cycling T&NK" = "NK cells",  # 2026-11-03: the 16th type of the intersect run, next to NK
             "Naive B" = "B cells", "Memory B" = "B cells",
             "Classical monocytes" = "Monocytes", "Non-classical monocytes" = "Monocytes",
             "cDC2" = "DCs", "pDC" = "DCs", "HSPC" = "Progenitor")
CT_ORDER <- names(LINEAGE)
ct_rank <- function(x) match(x, CT_ORDER)
LINEAGE_ORDER <- c("CD4 T", "CD8 T", "Unconventional T", "B cells", "NK cells", "Monocytes", "DCs", "Progenitor")
LINEAGE_COLORS <- c("CD4 T" = "#E41A1C", "CD8 T" = "#FF7F00", "Unconventional T" = "#FDBF6F", "B cells" = "#377EB8",
                    "NK cells" = "#4DAF4A", "Monocytes" = "#A65628", "DCs" = "#F781BF", "Progenitor" = "#999999")
CT_COLORS <- c("Naive CD4 T" = "#FB9A99", "Memory CD4 T" = "#E41A1C", "Regulatory T" = "#99000D",
               "Naive CD8 T" = "#FDD0A2", "CD8 effector memory T" = "#FF7F00", "CD8 effector T (TEMRA)" = "#A63603",
               "MAIT and γδ T" = "#E6C229", "NK" = "#4DAF4A", "Cycling T&NK" = "#A1D99B",
               "Naive B" = "#9ECAE1", "Memory B" = "#08519C",
               "Classical monocytes" = "#A65628", "Non-classical monocytes" = "#D9A066",
               "cDC2" = "#F781BF", "pDC" = "#984EA3", "HSPC" = "#999999")

# UMAP panels (Fig 2a-c, supplementary atlas figure a-c), the same size: same cell size, no
# coord_fixed. Legends at 5 pt (LEG_PT).
LEG_PT <- 5
UMAP_MM <- 46  # height of a row of UMAP panels
no_legend <- theme(legend.position = "none")
umap_panel <- function(d, colour, values, title, axes = "UMAP", legend_name = NULL) {
  ggplot(d, aes(.data[[paste0(axes, "_1")]], .data[[paste0(axes, "_2")]], colour = .data[[colour]])) +
    ggrastr::rasterise(geom_point(size = 0.2, alpha = 0.4, stroke = 0), dpi = 600) +
    scale_colour_manual(values = values, name = legend_name) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_nature() +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.key.size = unit(0.2, "cm"),
          legend.key.spacing.y = unit(0.5, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = LEG_PT, margin = margin(l = 1.5)),
          legend.title = element_text(size = LEG_PT, hjust = 0.5),
          axis.text = element_blank(),   # UMAP coordinates carry no meaning: axis lines and titles, no tick labels
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          panel.grid.major = element_blank(),   # theme_nature() sets the major grid explicitly
          plot.margin = margin(12, 12, 12, 12)) +
    panel_title(title)
}

# =============================================================================
# Results read once, used by the figures and the tables
# =============================================================================
# DEG rule of the paper: adjusted p < 0.05 and |log2FC| >= 1
is_deg <- function(padj, lfc) !is.na(padj) & padj < 0.05 & abs(lfc) >= 1
# differential expression vs HC (step 04), every gene tested
de <- rbindlist(lapply(c("SLE", "SjD"), function(d) read_ct(file.path(DATA, paste0("de_", d, "_vs_control.tsv")))[, disease := d]))
deg <- de[is_deg(adj_p_value, log_fc)]
TESTED <- intersect(CT_ORDER, unique(de$cell_type))  # cell types with a DE fit (HSPC has too few controls)
# composition (step 05): the one scCODA model, three contrasts. SLE vs SjD is the fit with base level SLE
# (SjD_vs_SLE) with every estimate negated and the interval bounds swapped, so that positive means higher in SLE. With the
# references each effect is credible with (step 05b): all_references, the rule of the text, Fig 2e and the per-cohort
# figure.
COMP <- data.table(contrast = c("SLE vs HC", "SjD vs HC", "SLE vs SjD"), fit = c("SLE_vs_control", "SjD_vs_control", "SjD_vs_SLE"),
                   base = c("control", "control", "SLE"), level = c("SLE", "SjD", "SjD"), sign = c(1, 1, -1))
reference <- readLines(file.path(DATA, "composition_reference.txt"), n = 1)
comp <- rbindlist(lapply(seq_len(nrow(COMP)), function(i) {
  x <- COMP[i]; s <- x$sign
  ef <- read_ct(file.path(DATA, sprintf("composition_%s_effects.tsv", x$fit)))
  cr <- read_ct(file.path(DATA, sprintf("composition_%s_credible.tsv", x$fit)))[
    Covariate == sprintf("C(condition, Treatment('%s'))[T.%s]", x$base, x$level)]
  stopifnot(nrow(cr) == nrow(ef))
  lo <- ef[["ETI 89% lower"]]; hi <- ef[["ETI 89% upper"]]
  data.table(contrast = x$contrast, cell_type = ef$cell_type, final_parameter = s * ef[["Final Parameter"]],
             eti_lower = if (s > 0) lo else -hi, eti_upper = if (s > 0) hi else -lo, sd = ef$SD,
             inclusion = ef[["Inclusion probability"]], log2fc = s * ef[["log2-fold change"]],
             credible = as.logical(cr$credible[match(ef$cell_type, cr[["Cell Type"]])]))
}))
sens <- read_ct(file.path(DATA, "reference_sensitivity.tsv"))[, .(contrast, cell_type, credible_with, all_references = all_references == "yes")]
comp <- merge(comp, sens, by = c("contrast", "cell_type"), all.x = TRUE, sort = FALSE)
# SLE vs SjD indirect comparison with its equivalence test and classes (07_figure_data.py)
ind <- read_gz("indirect.tsv.gz")
ind[, class := fcoalesce(class, "")]
# gene set enrichment (step 06; SLE vs SjD: 07_figure_data.py)
gsea <- rbindlist(list(read_ct(file.path(DATA, "gsea_SLE_vs_control.tsv"))[, `:=`(contrast = "SLE vs HC", block = "SLE")],
                       read_ct(file.path(DATA, "gsea_SjD_vs_control.tsv"))[, `:=`(contrast = "SjD vs HC", block = "SjD")],
                       read_ct(file.path(DATA, "indirect_gsea.tsv"))[, `:=`(contrast = "SLE vs SjD", block = "SLE vs SjD")]),
                  use.names = TRUE, fill = TRUE)

# =============================================================================
# Table 1: datasets
# =============================================================================
ACCESSION <- c("Perez" = "GEO GSE174188; CELLxGENE 436154da", "GSE135779_adult" = "GEO GSE135779 (adult arm)",
               "GSE253568" = "GEO GSE253568", "GSE316482" = "GEO GSE316482*", "HRA003613" = "GSA-Human HRA003613")
# Reference: the publication's citation number (Vancouver), in order of first citation in manuscript/main.tex (perez2022,
# neharbelaid2020, mcdermott2025, xiang2023; checked 2026-11-05): change them when the citations of main.tex change.
# GSE316482 has no publication.
REFERENCE <- c("Perez" = "6", "GSE135779_adult" = "8", "GSE253568" = "11", "GSE316482" = "GEO deposit only", "HRA003613" = "13")
tab <- data.frame(`Repository/accession` = ACCESSION[t1$dataset], Disease = t1$disease, Cases = t1$cases,
                  Controls = t1$controls, Cells = t1$cells_analysed,
                  `Median genes per cell` = t1$median_genes_per_cell,
                  Reference = unname(REFERENCE[t1$dataset]), check.names = FALSE)
tab <- rbind(tab, data.frame(`Repository/accession` = "Total", Disease = NA, Cases = sum(t1$cases), Controls = sum(t1$controls),
                             Cells = sum(t1$cells_analysed), `Median genes per cell` = NA, Reference = NA,
                             check.names = FALSE))
write_table(tab, file.path(TABLE_OUT, "table1.xlsx"), "Table 1", "Table 1. scRNA-seq datasets of the atlas",
            c("Cases and controls are donors; cells are those of the atlas, after quality control and the removal of low-quality and platelet clusters.",
              "* Deposited counts include intronic reads; the other datasets count exonic reads only."),
            int_cols = c("Cases", "Controls", "Cells", "Median genes per cell"), widths = c(36, 9, 8, 9, 14, 21, 17))

# =============================================================================
# Supplementary tables
# =============================================================================
# --- QC per dataset: the recorded counts of step 01, the atlas and the annotation exclusion; the per-sample
# thresholds are not recorded, so the table has none. Doublet rate: of the cells that entered doublet detection;
# annotation exclusion: of the atlas cells. ---
qc <- read_ct(file.path(DATA, "qc_counts.tsv"))
qc <- qc[match(names(DATASET_LABEL), dataset)]
stopifnot(!anyNA(qc$dataset), qc$after_mt - qc$doublets - qc$empty_after_ambient == qc$kept,
          qc$in_atlas - qc$excluded_annotation == qc$final)
qt <- data.frame(Dataset = DATASET_LABEL[qc$dataset], `Samples in the atlas` = qc$samples, `Cells before QC` = qc$before_qc,
                 `After the outlier filter` = qc$after_outlier, `After the mitochondrial filter` = qc$after_mt,
                 `Doublets removed` = ifelse(qc$doublets_run, qc$doublets, NA),
                 `Doublet rate (%)` = ifelse(qc$doublets_run, round(100 * qc$doublets / qc$after_mt, 1), NA),
                 `Kept after QC` = qc$kept, `In the atlas` = qc$in_atlas, `Excluded at annotation` = qc$excluded_annotation,
                 `Excluded at annotation (%)` = round(100 * qc$excluded_annotation / qc$in_atlas, 1),
                 `Final cells` = qc$final, check.names = FALSE)
tot <- qt[1, ]; tot[] <- NA; tot$Dataset <- "Total"
for (k in setdiff(names(qt), c("Dataset", "Doublet rate (%)", "Excluded at annotation (%)"))) tot[[k]] <- sum(qt[[k]], na.rm = TRUE)
tot$`Excluded at annotation (%)` <- round(100 * tot$`Excluded at annotation` / tot$`In the atlas`, 1)
qt <- rbind(qt, tot)
write_supp_table(qt, "supp_table_qc.xlsx", "Quality control: cells per dataset after each filter",
                 c(paste0("Outlier filter: cells beyond 5 median absolute deviations (MADs) of the sample median in log total counts, log detected genes or the share of counts in the 20 most expressed genes. ",
                          "Mitochondrial filter: cells above 3 MADs of the sample median in the logit mitochondrial fraction and with more than 8% mitochondrial counts. ",
                          "Doublets: scDblFinder per library; doublet rate of the cells that passed both filters. Excluded at annotation: clusters labelled low quality or platelets."),
                   "HRA003613: doublet detection not run; the authors removed doublets before deposit.",
                   sprintf("GSE174188: %d samples passed quality control; the %d later-visit samples of donors with more than one visit were removed at integration, leaving the first visit of each donor.",
                           qc[dataset == "Perez", samples_qc], qc[dataset == "Perez", samples_qc - samples]),
                   "HRA003613: 1 cell with no counts in the genes of the atlas was removed at integration.",
                   "GSE316482: 1 cell empty after ambient RNA correction was removed before doublet detection."),
                 int_cols = setdiff(names(qt), c("Dataset", "Doublet rate (%)", "Excluded at annotation (%)")))
log_msg("QC table: final cells", sum(qc$final), "| samples", sum(qc$samples))

# --- composition ---
ct_tab <- comp[order(match(contrast, COMP$contrast), ct_rank(cell_type))]
write_supp_table(data.frame(Contrast = ct_tab$contrast, `Cell type` = ct_tab$cell_type,
                            `Reference cell type` = ifelse(ct_tab$cell_type == reference, "yes", ""),
                            `Final parameter` = ct_tab$final_parameter, `ETI 89% lower` = ct_tab$eti_lower,
                            `ETI 89% upper` = ct_tab$eti_upper, SD = ct_tab$sd, `Inclusion probability` = ct_tab$inclusion,
                            `log2 fold change` = ct_tab$log2fc, `Credible (FDR 0.2)` = ifelse(ct_tab$credible, "yes", "no"),
                            `Credible with references` = ifelse(ct_tab$cell_type == reference, "", ct_tab$credible_with),
                            check.names = FALSE),
                 "supp_table_composition.xlsx", "Cell-type composition",
                 c("SLE vs HC and SjD vs HC: fit with healthy controls (HC) as the base level; SLE vs SjD: fit with SLE as the base level, estimates negated so that positive means higher in SLE.",
                   "Final parameter: scCODA effect on the log scale, 0 when the effect is not included; ETI: equal-tailed interval; credible: inclusion at a false discovery rate of 0.2.",
                   "log2 fold change: scCODA's log2 fold change of the expected proportion, which includes the shift relative to the reference cell type.",
                   "Credible with references: references (cDC2, regulatory T, non-classical monocytes, NK, CD8 effector memory T, pDC; the cell type itself excluded) with which the effect is credible; effects credible with every reference are reported."),
                 dec_cols = c("Final parameter", "ETI 89% lower", "ETI 89% upper", "SD", "Inclusion probability", "log2 fold change"))

# --- differential expression: every gene tested, per disease ---
for (d in c("SLE", "SjD")) {
  x <- de[disease == d]
  t <- data.frame(`Cell type` = x$cell_type, Gene = x$variable, `Base mean` = x$baseMean, log2FC = x$log_fc,
                  `log2FC SE` = x$lfcSE, `Wald statistic` = x$stat, P = x$p_value, `Adjusted P` = x$adj_p_value,
                  DEG = ifelse(is_deg(x$adj_p_value, x$log_fc), "yes", "no"),
                  Cases = x$n_case, Controls = x$n_control, Design = x$design, check.names = FALSE) %>%
    arrange(ct_rank(`Cell type`), P)
  write_supp_table(t, sprintf("supp_table_de_%s.xlsx", tolower(d)), sprintf("Differential expression, %s vs healthy controls", d),
                   c("Every gene tested in each cell type, after the edgeR filterByExpr filter and the filter keeping genes expressed in at least 10% of the cells in at least two pseudobulk samples; Adjusted P: Benjamini-Hochberg within cell type after independent filtering (empty when the gene was filtered out).",
                     "DEG: adjusted P < 0.05 and |log2FC| >= 1. Cases and Controls: donors with a pseudobulk sample in the cell type. Design: the model of the cell type, terms with one level dropped."),
                   int_cols = c("Cases", "Controls"))
}

# --- SLE vs SjD: the indirect comparison per gene and cell type (Fig 3b; Methods), every gene tested in both diseases ---
it <- data.frame(`Cell type` = ind$cell_type, Gene = ind$variable, `log2FC SLE` = ind$log_fc_sle, `log2FC SE SLE` = ind$lfcSE_sle,
                 `log2FC SjD` = ind$log_fc_sjd, `log2FC SE SjD` = ind$lfcSE_sjd, `Difference (SLE - SjD)` = ind$diff,
                 `Difference SE` = ind$se, z = ind$z, P = ind$pvalue, `Adjusted P` = ind$padj, `Equivalence P` = ind$p_tost,
                 `Adjusted equivalence P` = ind$padj_tost, `DEG in SLE` = ifelse(as.logical(ind$deg_sle), "yes", "no"),
                 `DEG in SjD` = ifelse(as.logical(ind$deg_sjd), "yes", "no"), Class = ind$class, check.names = FALSE) %>%
  arrange(ct_rank(`Cell type`), P)
write_supp_table(it, "supp_table_sle_vs_sjd.xlsx", "SLE vs SjD comparison of gene expression",
                 c("Every gene tested in both diseases within a cell type. Difference: log2FC SLE - log2FC SjD, each against the healthy controls of its own datasets; Difference SE: square root of the sum of the two squared SE; z: Difference / Difference SE; P: two-sided, standard normal; Adjusted P: Benjamini-Hochberg within cell type.",
                   "Equivalence P: two one-sided tests with an equivalence margin of 1 log2FC (|Difference| < 1), the larger of the two one-sided P; Adjusted equivalence P: Benjamini-Hochberg within cell type.",
                   "DEG in SLE, DEG in SjD: adjusted P < 0.05 and |log2FC| >= 1 against healthy controls. Class, for genes differentially expressed in either disease: Higher in SLE or Higher in SjD, adjusted P < 0.05 and |Difference| >= 1; Equivalent, adjusted equivalence P < 0.05; Inconclusive, neither."),
                 dec_cols = c("log2FC SLE", "log2FC SE SLE", "log2FC SjD", "log2FC SE SjD", "Difference (SLE - SjD)", "Difference SE", "z"))
log_msg("SLE vs SjD table:", nrow(it), "rows;", paste(names(table(it$Class[it$Class != ""])), table(it$Class[it$Class != ""]), collapse = ", "))

# --- GSEA: three contrasts ---
write_supp_table(data.frame(Contrast = gsea$contrast, `Cell type` = gsea$cell_type, `Gene set` = gsea$gene_set, NES = gsea$score,
                            P = gsea$pvalue, `Adjusted P` = gsea$padj, check.names = FALSE) %>%
                   arrange(match(Contrast, COMP$contrast), ct_rank(`Cell type`), P),
                 "supp_table_gsea.xlsx", "Gene set enrichment",
                 c("SLE vs HC and SjD vs HC: genes ranked by the PyDESeq2 Wald statistic; SLE vs SjD: by the z statistic of the indirect comparison (positive NES: higher in SLE).",
                   "Adjusted P: Benjamini-Hochberg across all cell type and gene set tests of a contrast."),
                 dec_cols = "NES")

# =============================================================================
# Fig 2: atlas (a-c UMAPs, d lineage composition, e scCODA effects)
# =============================================================================
cells <- read_gz("cells.tsv.gz")
cells[, condition := factor(cond_label(condition), levels = names(CONDITION_COLORS))]
cells[, dataset_label := factor(DATASET_LABEL[dataset], levels = DATASET_LABEL)]
cells[, cell_type := factor(cell_type, levels = CT_ORDER)]
set.seed(42); cells <- cells[sample(.N)]  # shuffle to avoid overplotting bias
log_msg("Fig 2:", nrow(cells), "cells")

# a-c: each legend in one column right of its panel; legends have no titles
p2a <- umap_panel(cells, "dataset_label", DATASET_COLORS, "Cells by dataset") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
p2b <- umap_panel(cells, "condition", CONDITION_COLORS, "Cells by disease group") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
p2c <- umap_panel(cells, "cell_type", CT_COLORS, "Cells by cell type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
legend_of <- function(p)  # legend left-aligned and vertically centred in its cell, right of its panel
  ggdraw(get_plot_component(p + theme(legend.justification.bottom = c(0, 0.5)), "guide-box-bottom"))
leg_a <- legend_of(p2a); leg_b <- legend_of(p2b); leg_c <- legend_of(p2c)

# d: cell-type composition per group, the 15 cell types in the colours and order of 2c (Naive CD4 T on top), with a
# legend styled as 2c's in one column right of the bars
# (segments without stroke under one rounded outline per bar)
ctp <- cells[, .N, by = .(condition, cell_type)][, pct := 100 * N / sum(N), by = condition]
ctp[, cell_type := factor(cell_type, levels = rev(CT_ORDER))]
ctp_tot <- as.data.frame(ctp[, .(pct = sum(pct)), by = condition])
p2d <- ggplot(ctp, aes(x = condition, y = pct)) +
  geom_col(aes(fill = cell_type), width = 0.7, color = NA, position = position_stack(reverse = TRUE),
           key_glyph = "point") +  # round keys, as the 2c legend
  stack_outline(ctp_tot, "condition", "pct", width = 0.7, size = 0.3) +
  scale_fill_manual(values = CT_COLORS, breaks = CT_ORDER, name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_nature() +
  theme(axis.text = element_text(size = MIN_PT),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.2, "cm"),
        legend.key.spacing.y = unit(0.5, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = LEG_PT, margin = margin(l = 1.5)),
        plot.margin = margin(12, 12, 12, 12)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(shape = 21, size = 1.5, colour = NA))) +
  panel_title("Cell-type composition")

# e: scCODA effects, three panels side by side, SLE vs HC, SjD vs HC and SLE vs SjD, sharing a y-axis of cell types (the
# reference cell type marked), scCODA's log2 fold change on the x-axis (the default of pertpy's effects bar plot, as the
# sc-best-practices chapter shows it; breaks -1, 0, 1; an axis line under each panel) with a line at 0; a point per effect
# credible at FDR 0.2 with every reference cell type of the sensitivity analysis (Methods; the log2FC of a
# non-credible cell type is only the compositional shift relative to the reference). Points coloured by sign, the same two
# colours in every panel (those of Fig 4's NES scale).
COMP_COLORS <- c("log2FC > 0" = "#B2182B", "log2FC < 0" = "#2166AC")
fig2e <- comp[, .(contrast = factor(contrast, levels = COMP$contrast), cell_type = factor(cell_type, levels = rev(CT_ORDER)),
                  log2fc, credible = credible & !is.na(all_references) & all_references)]
log_msg("Fig 2e credible effects:", paste(names(table(fig2e$contrast[fig2e$credible])), table(fig2e$contrast[fig2e$credible]), collapse = ", "))
ct_axis_label <- function(x) ifelse(x == reference, paste0(x, " (reference)"), x)
p2e <- ggplot(fig2e[credible == TRUE][, sign := factor(ifelse(log2fc > 0, "log2FC > 0", "log2FC < 0"), levels = names(COMP_COLORS))],
              aes(log2fc, cell_type, colour = sign)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(size = 1.6) +
  facet_wrap(~contrast, nrow = 1) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_discrete(drop = FALSE, labels = ct_axis_label) +
  scale_colour_manual(values = COMP_COLORS, drop = FALSE, name = NULL) +
  labs(x = "log2 fold change", y = NULL) + theme_nature() +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25), panel.grid.major.x = element_blank(),
        axis.line.x = element_line(colour = "black", linewidth = 0.4), axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
        panel.spacing.x = unit(3, "mm"), strip.text = element_text(size = MIN_PT, face = "bold"),
        legend.position = "bottom", legend.justification = "right", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0), legend.key.size = unit(3, "mm"), legend.text = element_text(size = MIN_PT),
        plot.margin = margin(4, 4, 4, 12)) +
  panel_title("Composition changes")  # methods (scCODA, reference, FDR) in the caption

# Compose in two rows: a, b, c, each with its legend in one column to its right, then d, e; row 2 shorter than row 1.
# The three UMAP cells share the width left by the legends equally (legend widths in mm, from the longest label).
ROW2_MM <- 66; LEG_W_MM <- c(a = 13, b = 7, c = 25)
UMAP_W <- (FIG_WIDTH * 25.4 - sum(LEG_W_MM)) / 3
row1 <- plot_grid(p2a + no_legend, leg_a, p2b + no_legend, leg_b, p2c + no_legend, leg_c, ncol = 6,
                  rel_widths = c(UMAP_W, LEG_W_MM["a"], UMAP_W, LEG_W_MM["b"], UMAP_W, LEG_W_MM["c"]),
                  labels = c("a", "", "b", "", "c", ""), label_size = LABEL_SIZE, label_fontface = "bold")
row2 <- plot_grid(p2d, p2e, ncol = 2, rel_widths = c(0.40, 0.60), labels = c("d", "e"),
                  label_size = LABEL_SIZE, label_fontface = "bold")
save_fig(plot_grid(row1, row2, ncol = 1, rel_heights = c(UMAP_MM, ROW2_MM)), FIG_OUT, "fig2_atlas.pdf", UMAP_MM + ROW2_MM)

# =============================================================================
# Fig 3: differential expression (a DEG counts, b SLE vs SjD for the genes of a)
# =============================================================================
deg_counts <- deg[, .N, by = .(cell_type, disease, direction = ifelse(log_fc > 0, "Up", "Down"))]
# 3b: genes changed vs HC in either disease (the DEG rule of panel a), by their class (07_figure_data.py: Wald test of the
# difference and two one-sided tests of equivalence, both adjusted within each cell type over all genes tested in both
# diseases; Methods). Equivalent: the colour of Fig 4b's header.
FIG3B_CATS <- c("Equivalent" = "#CC79A7", "SLE > SjD" = "#D55E00", "SjD > SLE" = "#0072B2", "Inconclusive" = "grey80")
FIG3B_OF_CLASS <- c("Equivalent" = "Equivalent", "Higher in SLE" = "SLE > SjD", "Higher in SjD" = "SjD > SLE", "Inconclusive" = "Inconclusive")
changed <- ind[class != ""][, category := unname(FIG3B_OF_CLASS[class])]
eq_counts <- changed[, .N, by = .(cell_type, category)]
log_msg("Fig 3b genes changed vs HC:", nrow(changed), "|",
        paste(names(FIG3B_CATS), sapply(names(FIG3B_CATS), function(k) sum(changed$category == k)), collapse = ", "))

# Diverging counts: up above zero, down below, at one x position per cell type. BAR_GROUP_W is the width a cell type's
# bars take (3a: SLE and SjD side by side, touching, SLE first; 3b: one bar), the rest is the gap between cell types.
# The y-axis is one scale through zero, its range trimmed to the data (not symmetric).
BAR_GROUP_W <- 0.7
diverging <- function(df, fill_col, fills, ylab, title, up_lab, down_lab, down, group = NULL, x_labels = TRUE) {
  df <- df %>% mutate(cell_type = factor(cell_type, levels = TESTED), n_plot = ifelse(direction %in% down, -N, N))
  lim <- c(min(df$n_plot, 0), max(df$n_plot, 0)) * 1.12
  bars <- if (is.null(group))
    ggchicklet::geom_chicklet(color = "black", size = 0.15, width = BAR_GROUP_W, radius = BAR_RADIUS)
  else
    ggchicklet::geom_chicklet(position = position_dodge(width = BAR_GROUP_W, preserve = "single"),
                              color = "black", size = 0.15, width = BAR_GROUP_W, radius = BAR_RADIUS)
  ggplot(df, aes(x = cell_type, y = n_plot, fill = .data[[fill_col]],
                 group = if (is.null(group)) NULL else .data[[group]])) +
    bars +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = fills, breaks = names(fills), name = NULL,
                      guide = guide_legend(nrow = 1, override.aes = list(color = NA))) +
    scale_y_continuous(limits = lim, labels = function(y) abs(y), expand = expansion(mult = c(0.02, 0.02))) +
    annotate("text", x = Inf, y = lim[2] * 0.9, label = up_lab, hjust = 1.1, vjust = 1, size = MIN_PT / .pt, fontface = "bold") +
    annotate("text", x = Inf, y = lim[1] * 0.9, label = down_lab, hjust = 1.1, vjust = 0, size = MIN_PT / .pt, fontface = "bold") +
    labs(x = NULL, y = ylab) +
    theme_nature() +
    theme(axis.text.x = if (x_labels) element_text(size = MIN_PT, angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = MIN_PT),
          plot.margin = margin(5, 5, 5, 28),
          legend.position = "top",
          legend.justification = "right",
          legend.text = element_text(size = MIN_PT),
          legend.key.size = unit(0.3, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major.x = element_blank()) +
    panel_title(title)
}
deg_counts[, fill := factor(paste(disease, direction), levels = c("SLE Down", "SLE Up", "SjD Down", "SjD Up"))]
deg_counts[, disease := factor(disease, levels = c("SLE", "SjD"))]  # SLE bar first within each cell type
p3a <- diverging(deg_counts, "fill", c("SLE Down" = "#E69F00", "SLE Up" = "#D55E00", "SjD Down" = "#56B4E9", "SjD Up" = "#0072B2"),
                 "Number of DEGs", "DEG counts", "Up", "Down", down = "Down", group = "disease", x_labels = FALSE)
# 3b: stacked bars per cell type (Equivalent at the bottom, Inconclusive on top), same x positions, bar width, outline, theme
# and margins as 3a
p3b <- ggplot(eq_counts %>% mutate(cell_type = factor(cell_type, levels = TESTED),
                                   category = factor(category, levels = names(FIG3B_CATS))),
              aes(x = cell_type, y = N, fill = category)) +
  ggchicklet::geom_chicklet(position = position_stack(reverse = TRUE), color = "black", size = 0.15,
                            width = BAR_GROUP_W, radius = BAR_RADIUS) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = FIG3B_CATS, breaks = names(FIG3B_CATS), name = NULL,
                    guide = guide_legend(nrow = 1, override.aes = list(color = NA))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Genes changed vs HC") +
  theme_nature() +
  theme(axis.text.x = element_text(size = MIN_PT, angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = MIN_PT),
        plot.margin = margin(5, 5, 5, 28),
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = MIN_PT),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.x = element_blank()) +
  panel_title("SLE vs SjD")
# 3a has no cell-type labels (3b's, same order, serve both); panels aligned so the bars line up
save_fig(plot_grid(p3a, p3b, ncol = 1, rel_heights = c(0.75, 1), align = "v", axis = "lr", labels = c("a", "b"),
                   label_size = LABEL_SIZE, label_fontface = "bold"), FIG_OUT, "fig3_de.pdf", 150)

# =============================================================================
# Fig 4: pathways (a NES vs healthy controls, SLE and SjD; b fgsea on the indirect comparison). Heatmap: white cell borders, 7 mm blocks, 45-degree column labels, bodies sized to fill the page width. Column labels at
# 5.5 pt and gaps of 3 and 5 mm: at the 175 mm page width a column is about 2.8 mm, and 7 pt labels at
# 45 degrees need about 3.5 mm, so they overlapped. The sign convention of b is in the caption only.
# =============================================================================
gsea[, pathway := format_pathway_label(gene_set)]
cts <- intersect(CT_ORDER, unique(gsea$cell_type))
# rows: pathways significant in at least one cell type of any block; ordered by total absolute mean NES
# (|mean over all cell types, SLE| + |mean, SjD|), descending
sig_paths <- unique(gsea[padj < 0.05, pathway])
ord <- gsea[block %in% c("SLE", "SjD") & pathway %in% sig_paths, .(m = mean(score, na.rm = TRUE)), by = .(pathway, block)][
  , .(tot = sum(abs(m))), by = pathway][order(-tot), pathway]
ord <- c(ord, setdiff(sig_paths, ord))
mat_of <- function(b) {
  m <- matrix(NA_real_, length(ord), length(cts), dimnames = list(ord, cts))
  s <- gsea[block == b & pathway %in% ord & cell_type %in% cts & padj < 0.05]  # white: adjusted p >= 0.05
  m[cbind(s$pathway, s$cell_type)] <- s$score
  m
}
m_sle <- mat_of("SLE"); m_sjd <- mat_of("SjD"); m_ind <- mat_of("SLE vs SjD")
max_abs <- min(4, max(abs(c(m_sle, m_sjd, m_ind)), na.rm = TRUE))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#2166AC", "white", "#B2182B"))
label_gp <- gpar(fontsize = MIN_PT, fontfamily = FONT_FAMILY, col = TEXT_COL)
col_label_gp <- gpar(fontsize = 5.5, fontfamily = FONT_FAMILY, col = TEXT_COL)
ROW_MM <- 2.5; GAP_MM <- 3; HT_GAP_MM <- 5; ROT <- 45
PAD_MM <- c(2, 2, if (DRAFT) 6 else 2, 2)  # bottom, left, top, right; the DRAFT mark sits in the top padding
block_anno <- function(name, labels, fills) do.call(HeatmapAnnotation, c(
  setNames(list(anno_block(gp = gpar(fill = fills, col = "black", lwd = FRAME_LWD), labels = labels,
                           labels_gp = gpar(col = "white", fontsize = TITLE_PT, fontface = "bold"),
                           height = unit(BLOCK_MM, "mm"))), paste0("block_", name)),
  list(show_annotation_name = FALSE, show_legend = FALSE)))
ht_common <- list(col = col_fun, na_col = "white", rect_gp = gpar(col = "white", lwd = 0.3), border = TRUE,
                  border_gp = gpar(col = "black", lwd = FRAME_LWD), cluster_rows = FALSE, cluster_columns = FALSE,
                  column_names_side = "bottom", column_names_rot = ROT, column_names_gp = col_label_gp,
                  column_title_gp = gpar(fontsize = TITLE_PT, fontface = "bold", col = TEXT_COL),
                  height = unit(length(ord) * ROW_MM, "mm"))
make_hts <- function(body_mm) {  # body_mm: total width of the three blocks' bodies
  w_col <- body_mm / (3 * length(cts))
  ht_a <- do.call(Heatmap, c(list(cbind(m_sle, m_sjd), name = "NES",
    column_split = factor(rep(c("SLE vs HC", "SjD vs HC"), each = length(cts)), levels = c("SLE vs HC", "SjD vs HC")),
    column_labels = rep(cts, 2), column_gap = unit(GAP_MM, "mm"),
    column_title = "Pathway enrichment",  # one title over both blocks, centred on them
    top_annotation = block_anno("a", c("SLE vs HC", "SjD vs HC"), DISEASE_COLORS),
    row_names_side = "left", row_names_gp = label_gp,
    width = unit(2 * length(cts) * w_col, "mm"),
    heatmap_legend_param = list(title_gp = label_gp, labels_gp = label_gp, direction = "horizontal",
                                title_position = "topcenter", legend_width = unit(2.5, "cm"))), ht_common))
  ht_b <- do.call(Heatmap, c(list(m_ind, name = "NES_b", show_heatmap_legend = FALSE, column_labels = cts,
    column_title = "Pathway differences", show_row_names = FALSE,  # title in the style of a (column_title_gp of ht_common)
    top_annotation = block_anno("b", "SLE vs SjD", "#CC79A7"),
    width = unit(length(cts) * w_col, "mm")), ht_common))
  ht_a + ht_b
}
draw_args <- list(ht_gap = unit(HT_GAP_MM, "mm"), heatmap_legend_side = "bottom",
                  align_heatmap_legend = "heatmap_center", padding = unit(PAD_MM, "mm"))
# Absolute sizes: estimate the bodies' width from the page width, then measure the drawn width and take
# the excess off the bodies so the page is exactly FIG_WIDTH
rn_mm <- convertWidth(max_text_width(ord, gp = label_gp), "mm", valueOnly = TRUE)
body_mm <- FIG_WIDTH * 25.4 - PAD_MM[2] - PAD_MM[4] - rn_mm -
  convertWidth(ht_opt("DIMNAME_PADDING"), "mm", valueOnly = TRUE) - GAP_MM - HT_GAP_MM
pdf(NULL)
drawn <- do.call(draw, c(list(make_hts(body_mm)), draw_args))
body_mm <- body_mm - (convertWidth(ComplexHeatmap:::width(drawn), "mm", valueOnly = TRUE) - FIG_WIDTH * 25.4)
drawn <- do.call(draw, c(list(make_hts(body_mm)), draw_args))
w_mm <- convertWidth(ComplexHeatmap:::width(drawn), "mm", valueOnly = TRUE)
h_mm <- convertHeight(ComplexHeatmap:::height(drawn), "mm", valueOnly = TRUE)
invisible(dev.off())
stopifnot(h_mm <= FIG_MAX_HEIGHT * 25.4)
f4 <- file.path(FIG_OUT, "fig4_pathways.pdf")
pdf_device(f4, width = w_mm / 25.4, height = h_mm / 25.4)
do.call(draw, c(list(make_hts(body_mm)), draw_args))
# panel labels level with the column titles: a at the page's left edge, b just left of its block
decorate_column_title("NES", slice = 1, {
  grid.text("a", x = unit(0, "npc") - unit(rn_mm, "mm") - ht_opt("DIMNAME_PADDING"), y = unit(1, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = LABEL_SIZE, fontface = "bold", fontfamily = FONT_FAMILY))
})
decorate_column_title("NES_b", slice = 1, {
  grid.text("b", x = unit(-2, "mm"), y = unit(1, "npc"), just = c("right", "top"),
            gp = gpar(fontsize = LABEL_SIZE, fontface = "bold", fontfamily = FONT_FAMILY))
})
if (DRAFT) grid.text(DRAFT_TEXT, x = unit(1, "npc") - unit(1, "mm"), y = unit(1, "npc") - unit(1, "mm"), just = c("right", "top"),
                     gp = gpar(col = "red", fontsize = MIN_PT, fontface = "bold", fontfamily = FONT_FAMILY))
invisible(dev.off())
log_msg("saved", f4, "|", round(w_mm, 1), "x", round(h_mm, 1), "mm |", length(ord), "pathways,", length(cts), "cell types")

# =============================================================================
# Supplementary figures
# =============================================================================
# --- per-cohort composition: per dataset and cell type, log2 of the mean proportion in its cases over the mean
# proportion in its own controls (donor proportions); a point marks an effect of the joint scCODA model for that
# dataset's disease (SLE vs HC or SjD vs HC) credible with every reference cell type (step 05b), the rule of the text and
# Fig 2e ---
cc <- read_ct(file.path(DATA, "composition_counts.tsv"))
cc[, prop := n / sum(n), by = donor_id]
pc <- cc[, .(mean_prop = mean(prop), donors = uniqueN(donor_id)), by = .(dataset, condition, cell_type)]
pc <- merge(pc[condition != "control"], pc[condition == "control"], by = c("dataset", "cell_type"), suffixes = c("", "_hc"))
pc[, log2_ratio := log2(mean_prop / mean_prop_hc)]
pc[, disease := condition]
pc <- merge(pc, comp[contrast %in% c("SLE vs HC", "SjD vs HC"), .(disease = sub(" vs HC", "", contrast), cell_type, credible = all_references)],
            by = c("disease", "cell_type"), all.x = TRUE)
stopifnot(!anyNA(pc$credible))
n_lab <- unique(pc[, .(dataset, disease, donors, donors_hc)])
n_lab[, label := sprintf("%s\n%d %s / %d HC", DATASET_LABEL[dataset], donors, disease, donors_hc)]
pc[, col := factor(n_lab$label[match(dataset, n_lab$dataset)], levels = n_lab$label[match(c(SLE_DS, SJD_DS), n_lab$dataset)])]
pc[, cell_type := factor(ifelse(cell_type == reference, paste(cell_type, "(reference)"), cell_type),
                         levels = rev(ifelse(CT_ORDER == reference, paste(CT_ORDER, "(reference)"), CT_ORDER)))]
pc[, disease := factor(disease, levels = c("SLE", "SjD"))]
lim <- min(3, ceiling(max(abs(pc$log2_ratio[is.finite(pc$log2_ratio)]))))
pcomp <- ggplot(pc, aes(col, cell_type, fill = pmax(pmin(log2_ratio, lim), -lim))) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = pc[credible == TRUE], shape = 16, size = 1.2, colour = "black") +
  facet_grid(. ~ disease, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim),
                       name = "log2 ratio vs HC") +  # definition in the figure legend
  labs(x = NULL, y = NULL) +
  theme_nature() +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(),  # theme_nature() sets the major grid explicitly
        axis.text.x = element_text(size = MIN_PT, lineheight = 0.9),
        strip.text = element_text(size = MIN_PT, face = "bold"), legend.position = "right",
        legend.key.height = unit(6, "mm"), legend.key.width = unit(3, "mm")) +
  panel_title("Composition per cohort")
save_fig(pcomp, SUPP_OUT, "supp_fig_composition_cohorts.pdf", 95)
log_msg("per-cohort composition: points mark joint-model effects credible with every reference:",
        paste(unique(pc[credible == TRUE, paste(disease, cell_type)]), collapse = "; "), "| log2 ratios capped at +-", lim)

# --- per-cohort DE: a, per-cohort log2FC compared between cohorts of the same disease, over the DEGs of the joint model
# tested in both cohorts; b, fixed-effect meta-analysis vs the joint model, every gene tested in both (joint DEGs in
# colour); c, I^2 of the joint DEGs tested in more than one cohort ---
meta <- read_gz("percohort_meta.tsv.gz")
meta[, joint_deg := is_deg(joint_adj_p_value, joint_log_fc)]
PAIRS <- c(combn(SLE_DS, 2, simplify = FALSE), combn(SJD_DS, 2, simplify = FALSE))
r_lab <- function(x, y) sprintf("r = %.2f\nn = %s", cor(x, y), format(length(x), big.mark = ","))
pair_plot <- function(p) {
  d <- if (p[1] %in% SLE_DS) "SLE" else "SjD"
  x <- meta[disease == d & joint_deg][, .(a = get(paste0("log_fc_", p[1])), b = get(paste0("log_fc_", p[2])))][!is.na(a) & !is.na(b)]
  m <- max(abs(c(x$a, x$b)))
  ggplot(x, aes(a, b)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_abline(linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.6, colour = DISEASE_COLORS[[d]]) +
    annotate("text", x = -m, y = m, label = r_lab(x$a, x$b), hjust = 0, vjust = 1, size = MIN_PT / .pt, family = FONT_FAMILY) +
    coord_equal(xlim = c(-m, m), ylim = c(-m, m)) +
    labs(x = paste("log2FC,", DATASET_LABEL[p[1]]), y = paste("log2FC,", DATASET_LABEL[p[2]])) +
    theme_nature() + panel_title(paste(d, "DEGs"))
}
meta_plot <- function(d) {
  x <- meta[disease == d & !is.na(joint_log_fc)]
  m <- max(abs(c(x$joint_log_fc, x$meta_log_fc)))
  ggplot(x, aes(joint_log_fc, meta_log_fc)) +
    geom_abline(linetype = "dashed", linewidth = 0.3) +
    ggrastr::geom_point_rast(data = x[joint_deg == FALSE], size = 0.2, colour = "grey75", raster.dpi = 300) +
    ggrastr::geom_point_rast(data = x[joint_deg == TRUE], size = 0.4, colour = DISEASE_COLORS[[d]], raster.dpi = 300) +
    annotate("text", x = -m, y = m, hjust = 0, vjust = 1, size = MIN_PT / .pt, family = FONT_FAMILY,
             label = sprintf("all genes: %s\nDEGs: %s", sub("\n", ", ", r_lab(x$joint_log_fc, x$meta_log_fc)),
                             sub("\n", ", ", r_lab(x$joint_log_fc[x$joint_deg], x$meta_log_fc[x$joint_deg])))) +
    coord_equal(xlim = c(-m, m), ylim = c(-m, m)) +
    labs(x = "log2FC, joint model", y = "log2FC, meta-analysis") +
    theme_nature() + panel_title(paste(d, "vs HC"))
}
i2 <- meta[joint_deg & k > 1][, disease := factor(disease, levels = c("SLE", "SjD"))]
pi2 <- ggplot(i2, aes(I2, fill = disease)) +
  geom_histogram(breaks = seq(0, 1, 0.1), position = position_dodge(), colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = DISEASE_COLORS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(I^2), y = "Joint-model DEGs") +
  theme_nature() + theme(legend.position = c(0.75, 0.85), legend.key.size = unit(3, "mm")) +
  panel_title("Heterogeneity between cohorts")
spaced <- theme(plot.margin = margin(6, 9, 6, 9))  # space between the panels
row_a <- plot_grid(plotlist = lapply(lapply(PAIRS, pair_plot), `+`, spaced), nrow = 1, labels = letters[seq_along(PAIRS)],
                   label_size = LABEL_SIZE, label_fontface = "bold")
row_b <- plot_grid(meta_plot("SLE") + spaced, meta_plot("SjD") + spaced, pi2 + spaced, nrow = 1,
                   labels = letters[length(PAIRS) + 1:3], label_size = LABEL_SIZE, label_fontface = "bold")
save_fig(plot_grid(row_a, NULL, row_b, ncol = 1, rel_heights = c(1, 0.08, 1.25)), SUPP_OUT, "supp_fig_de_cohorts.pdf", 130)
for (d in c("SLE", "SjD")) {
  x <- meta[disease == d & joint_deg]
  log_msg(d, "joint DEGs:", nrow(x), "| tested in > 1 cohort:", sum(x$k > 1), "| I2 > 0.5:", sum(x$I2 > 0.5, na.rm = TRUE),
          "| same sign in the meta-analysis:", sum(sign(x$meta_log_fc) == sign(x$joint_log_fc)))
}

# --- atlas: a, UMAP before integration by dataset; b, integrated UMAP of the controls by dataset, an equal
# random sample of control cells per dataset (the smallest dataset's number); c, integrated UMAP by annotation cluster
# (numbered in the cell-type order, excluded clusters last and black); a-c in one row as Fig 2a-c, the dataset legend in
# one row under a and b, the cluster key right of c; d, dot plot of config.SUPP_DOTPLOT_MARKERS by final cell type,
# full width below ---
ac <- read_gz("atlas_cells.tsv.gz")
setnames(ac, "label", "cell_type"); ac[, cell_type := ct_display(cell_type)]
ac[, dataset_label := factor(DATASET_LABEL[dataset], levels = DATASET_LABEL)]
set.seed(42); ac <- ac[sample(.N)]  # random drawing order
ctrl <- ac[condition == "control"]
n_b <- min(ctrl[, .N, by = dataset]$N)
set.seed(42); ctrl <- ctrl[ctrl[, .I[sample(.N, n_b)], by = dataset]$V1][sample(.N)]  # n_b cells per dataset, random order
EXCL <- c("low quality", "platelets")  # config.EXCLUDED_LABELS
cl <- unique(ac[, .(cluster, cell_type)])
stopifnot(!anyDuplicated(cl$cluster), all(cl$cell_type %in% c(CT_ORDER, EXCL)))
cl[, res := as.integer(sub("leiden_(\\d):.*", "\\1", cluster))][, id := as.integer(sub(".*:", "", cluster))]
setorder(cl[, rank := match(cell_type, c(CT_ORDER, EXCL))], rank, res, id)
cl[, number := seq_len(.N)]
ac[, number := factor(cl$number[match(cluster, cl$cluster)], levels = cl$number)]
kept <- cl[!cell_type %in% EXCL]
CL_COLORS <- c(setNames(scales::hue_pal()(nrow(kept)), kept$number), setNames(rep("black", nrow(cl) - nrow(kept)), cl[cell_type %in% EXCL, number]))
# a-c: the UMAP panels of Fig 2; the dataset legend in one row under a and b, serving both; the cluster key right of c
pa <- umap_panel(ac, "dataset_label", DATASET_COLORS, "Before integration", "umap_unint")
pb <- umap_panel(ctrl, "dataset_label", DATASET_COLORS, "Only HC donors, integrated", "umap") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), nrow = 1))  # one row, centred under a and b
cen <- ac[, .(x = median(umap_1), y = median(umap_2)), by = number]
pc_ <- umap_panel(ac, "number", CL_COLORS, "Annotation clusters", "umap") +
  ggrepel::geom_text_repel(data = cen, aes(x, y, label = number), inherit.aes = FALSE, size = LEG_PT / .pt, family = FONT_FAMILY,
                           fontface = "bold", seed = 42, box.padding = 0.08, point.padding = 0, min.segment.length = 0.1,
                           segment.size = 0.2, max.overlaps = Inf, bg.colour = "white", bg.r = 0.12)  # Fig 2 labels no clusters: ggrepel, short segments; a white outline keeps numbers on the black clusters readable
leg_ab <- ggdraw(get_plot_component(pb, "guide-box-bottom"))  # centred under a and b
# cluster key right of c: one entry per line at the line pitch of the legends (key 0.2 cm + 0.5 pt), then a blank line
# and the note, vertically centred
key <- cl[, .(n = paste(number, collapse = ", ")), by = cell_type][
  , txt := sprintf("%s: %s", ifelse(cell_type %in% EXCL, sub("^(.)", "\\U\\1", cell_type, perl = TRUE), cell_type), n)]  # "Low quality", "Platelets"
key_lines <- c(key$txt, "", "Black: excluded at annotation")
pitch <- (2 + 0.5 / 72.27 * 25.4) / UMAP_MM  # line pitch as a fraction of the row height
key_y <- 0.5 + (rev(seq_along(key_lines)) - (length(key_lines) + 1) / 2) * pitch
key_c <- ggdraw() + draw_text(key_lines, x = 0, y = key_y, hjust = 0, size = LEG_PT, family = FONT_FAMILY)
# widest line, measured on the figures' device (the null device's metrics are narrower), plus 3 mm to the page edge
pdf_device(tf <- tempfile(fileext = ".pdf"), 7, 7)
key_w_mm <- max(sapply(key_lines, function(x) convertWidth(grobWidth(textGrob(x, gp = gpar(fontsize = LEG_PT, fontfamily = FONT_FAMILY))), "mm", valueOnly = TRUE))) + 3
invisible(dev.off()); unlink(tf)
dot <- read_ct(file.path(DATA, "atlas_dotplot.tsv"))
dot[, cell_type := ct_display(cell_type)]
mk <- unique(dot$gene)  # in the order of config.SUPP_DOTPLOT_MARKERS (07_figure_data.py writes them so)
stopifnot(setequal(dot$cell_type, CT_ORDER), dot[, .N, by = gene]$N == length(CT_ORDER))
dot[, label := ifelse(n_datasets < n_datasets_all, sprintf("%s (%d/%d)", gene, n_datasets, n_datasets_all), gene)]
dot[, scaled := (mean - min(mean)) / (max(mean) - min(mean)), by = gene]  # per gene, 0-1 over the cell types (sc.pl.dotplot standard_scale = "var")
dot[, label := factor(label, levels = unique(dot$label[match(mk, dot$gene)]))][, cell_type := factor(cell_type, levels = rev(CT_ORDER))]
pd_ <- ggplot(dot, aes(label, cell_type, size = 100 * fraction, colour = scaled)) +
  geom_point(stroke = 0) +
  scale_size_area(max_size = 2.5, name = "Cells with\ncounts (%)", breaks = c(25, 50, 75, 100)) +
  scale_colour_gradient(low = "grey90", high = "#08306B", name = "Mean expression\n(scaled per gene)", breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = NULL) + theme_nature() +
  theme(axis.text.x = element_text(size = LEG_PT, angle = 45, hjust = 1, vjust = 1, face = "italic"), axis.text.y = element_text(size = LEG_PT),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
        legend.position = "right", legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
        plot.margin = margin(4, 4, 4, 12)) +  # as Fig 2e: room for the panel label d
  panel_title("Annotation markers")
# Compose as Fig 2's row of UMAPs: a, b, c the same size; the dataset legend in one row centred under a and b, the key
# right of c; d below at full width
DOT_MM <- 60; LEG_ROW_MM <- 6
UMAP_W <- (FIG_WIDTH * 25.4 - key_w_mm) / 3
ab <- plot_grid(plot_grid(pa + no_legend, pb + no_legend, nrow = 1, labels = c("a", "b"), label_size = LABEL_SIZE, label_fontface = "bold"),
                leg_ab, ncol = 1, rel_heights = c(UMAP_MM, LEG_ROW_MM))
c_key <- plot_grid(plot_grid(pc_ + no_legend, key_c, nrow = 1, rel_widths = c(UMAP_W, key_w_mm), labels = c("c", ""),
                             label_size = LABEL_SIZE, label_fontface = "bold"),
                   NULL, ncol = 1, rel_heights = c(UMAP_MM, LEG_ROW_MM))  # an empty strip level with the legend: c as tall as a and b
row1 <- plot_grid(ab, c_key, nrow = 1, rel_widths = c(2 * UMAP_W, UMAP_W + key_w_mm))
save_fig(plot_grid(row1, plot_grid(pd_, labels = "d", label_size = LABEL_SIZE, label_fontface = "bold"),
                   ncol = 1, rel_heights = c(UMAP_MM + LEG_ROW_MM, DOT_MM)), SUPP_OUT, "supp_fig_atlas.pdf", UMAP_MM + LEG_ROW_MM + DOT_MM)
log_msg("atlas figure: panel b", n_b, "control cells per dataset; key", round(key_w_mm, 1), "mm wide; UMAP cells", round(UMAP_W, 1), "mm")
log_msg("atlas figure:", nrow(ac), "cells,", nrow(ctrl), "control cells drawn,", nrow(cl), "clusters (",
        sum(cl$cell_type %in% EXCL), "excluded ),", length(mk), "markers:",
        paste(unique(dot[n_datasets < n_datasets_all, as.character(label)]), collapse = ", "), "outside the shared genes")
log_msg("done")
