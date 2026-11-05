#!/usr/bin/env Rscript
# Manuscript figures 2-4 and Table 1. Not a pipeline step and not a notebook port. Shared style in figure_style.R
# (theme_nature(), panel_title(), the bar conventions BAR_RADIUS and stack_outline(), the coloured facet strips
# class_strips() and the heatmap layout), sized for the Springer Nature artwork guide (175 mm wide, text 5-7 pt). Fig 1 is
# TikZ (fig1_pipeline.tex).
# Usage: run by make_figures.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration; FIG_DRAFT=1 adds a
#        DRAFT mark)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/: 07_figures/data/ (07_figure_data.py), 04_differential_expression/,
#        05_composition/, 06_pathways/
# Output, overwritten on each run: manuscript/figures/ ($FIG_OUT) fig2_atlas.pdf, fig3_de.pdf, fig4_pathways.pdf;
#         manuscript/tables/ ($TABLE_OUT) table1.xlsx
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot); library(grid)
  library(ComplexHeatmap); library(circlize); library(data.table)
})

ATLAS_RUN <- Sys.getenv("ATLAS_RUN", "02_integration")
ATLAS_ROOT <- Sys.getenv("ATLAS_ROOT", ".")  # set by run.sh
RUN  <- file.path(Sys.getenv("PIPELINE_OUT", file.path(ATLAS_ROOT, "out")), "py", "atlas", ATLAS_RUN)
DATA <- file.path(RUN, "07_figures", "data")
OUT <- Sys.getenv("FIG_OUT", file.path(ATLAS_ROOT, "manuscript", "figures"))
TABLE_OUT <- Sys.getenv("TABLE_OUT", file.path(ATLAS_ROOT, "manuscript", "tables"))
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", file.path(ATLAS_ROOT, "pipeline_py"))  # datasets.tsv
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(TABLE_OUT, recursive = TRUE, showWarnings = FALSE)
DRAFT <- Sys.getenv("FIG_DRAFT", "0") == "1"  # FIG_DRAFT=1 adds a red DRAFT | run | date mark
DRAFT_TEXT <- paste0("DRAFT | ", ATLAS_RUN, " | ", format(Sys.Date()))
log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
pdf(NULL)  # null default device: plot measurements otherwise write Rplots.pdf into the working directory (the repo)

# Style, colours, datasets and cell types: figure_style.R (shared with 09_supplementary.R)
source(file.path(PIPELINE_DIR, "figure_style.R"))

# One cell-type order in every figure: the lineage order of the Fig 2c legend (CT_ORDER)
de <- bind_rows(lapply(c("SLE", "SjD"), function(d)
  read_ct(file.path(RUN, "04_differential_expression", paste0(d, "_vs_control"), "all", "de_all_cell_types.tsv"))[, disease := d]))
de <- as.data.table(de)
deg <- de[!is.na(adj_p_value) & adj_p_value < 0.05 & abs(log_fc) >= 1]
CT_FIG_ORDER <- CT_ORDER
TESTED <- intersect(CT_FIG_ORDER, unique(de$cell_type))  # cell types with a DE fit (HSPC has too few controls)

# =============================================================================
# Table 1: datasets, one Excel workbook with a title row, flat values, no formulas
# =============================================================================
t1 <- read_ct(file.path(DATA, "table1.tsv"))
ACCESSION <- c("Perez" = "GEO GSE174188; CELLxGENE 436154da", "GSE135779_adult" = "GEO GSE135779 (adult arm)",
               "GSE253568" = "GEO GSE253568", "GSE316482" = "GEO GSE316482*", "HRA003613" = "GSA-Human HRA003613")
# Reference: the publication's citation number (Vancouver), in order of first citation in manuscript/main.tex, read at
# build time so the numbers follow the text; GSE316482 has no publication
BIBKEY <- c("Perez" = "perez2022", "GSE135779_adult" = "neharbelaid2020", "GSE253568" = "mcdermott2025",
            "GSE316482" = NA, "HRA003613" = "xiang2023")
tex <- paste(readLines(file.path(ATLAS_ROOT, "manuscript", "main.tex"), warn = FALSE), collapse = "\n")
tex <- sub("^.*?\\\\begin\\{document\\}", "", tex)
cited <- unique(trimws(unlist(strsplit(regmatches(tex, gregexpr("\\\\cite[pt]?\\*?(?:\\[[^]]*\\])?\\{\\K[^}]*", tex, perl = TRUE))[[1]], ","))))
cite_no <- function(k) if (is.na(k)) "GEO deposit only" else if (k %in% cited) as.character(match(k, cited)) else
  stop("Table 1: ", k, " is not cited in main.tex")
tab <- data.frame(`Repository/accession` = ACCESSION[t1$dataset], Disease = t1$disease, Cases = t1$cases,
                  Controls = t1$controls, Cells = t1$cells_analysed,
                  `Median genes per cell` = t1$median_genes_per_cell,
                  Reference = vapply(BIBKEY[t1$dataset], cite_no, ""), check.names = FALSE)
log_msg("Table 1 citation numbers:", paste(t1$dataset, tab$Reference, sep = " = ", collapse = ", "))
tab <- rbind(tab, data.frame(`Repository/accession` = "Total", Disease = NA, Cases = sum(t1$cases), Controls = sum(t1$controls),
                             Cells = sum(t1$cells_analysed), `Median genes per cell` = NA, Reference = NA,
                             check.names = FALSE))
notes <- c("Cases and controls are donors; cells are those of the atlas, after quality control and the removal of low-quality and platelet clusters.",
           "* Deposited counts include intronic reads; the other datasets count exonic reads only.")
r0 <- if (DRAFT) 3 else 2                                   # header row
wb <- openxlsx2::wb_workbook()
wb <- openxlsx2::wb_add_worksheet(wb, "Table 1")
wb <- openxlsx2::wb_add_data(wb, x = "Table 1. scRNA-seq datasets of the atlas", start_row = 1, start_col = 1)
wb <- openxlsx2::wb_add_font(wb, dims = "A1", bold = TRUE)
if (DRAFT) {
  wb <- openxlsx2::wb_add_data(wb, x = DRAFT_TEXT, start_row = 2, start_col = 1)
  wb <- openxlsx2::wb_add_font(wb, dims = "A2", bold = TRUE, color = openxlsx2::wb_color(hex = "FFFF0000"))
}
wb <- openxlsx2::wb_add_data(wb, x = tab, start_row = r0, start_col = 1, na.strings = "")
wb <- openxlsx2::wb_add_font(wb, dims = openxlsx2::wb_dims(rows = r0, cols = seq_len(ncol(tab))), bold = TRUE)
wb <- openxlsx2::wb_add_numfmt(wb, dims = openxlsx2::wb_dims(rows = r0 + seq_len(nrow(tab)), cols = 3:6), numfmt = "#,##0")
wb <- openxlsx2::wb_add_data(wb, x = notes, start_row = r0 + nrow(tab) + 2, start_col = 1)
wb <- openxlsx2::wb_set_col_widths(wb, cols = seq_len(ncol(tab)), widths = c(36, 9, 8, 9, 14, 21, 17))
t1_file <- file.path(TABLE_OUT, "table1.xlsx")
openxlsx2::wb_save(wb, t1_file, overwrite = TRUE)
log_msg("saved", t1_file)

# =============================================================================
# Fig 2: atlas (a-c UMAPs, d lineage composition, e scCODA effects)
# =============================================================================
cells <- read_ct(cmd = paste("zcat", shQuote(file.path(DATA, "cells.tsv.gz"))))  # no R.utils in the image
cells[, condition := factor(cond_label(condition), levels = names(CONDITION_COLORS))]
cells[, dataset_label := factor(DATASET_LABEL[dataset], levels = DATASET_LABEL)]
cells[, cell_type := factor(cell_type, levels = CT_ORDER)]
set.seed(42); cells <- cells[sample(.N)]  # shuffle to avoid overplotting bias
log_msg("Fig 2:", nrow(cells), "cells")

# a-c UMAP panels, the same size: same cell size, no coord_fixed. Legends at 5 pt
# (LEG_PT), each in one column right of its panel.
LEG_PT <- 5
umap_panel <- function(colour, values, title, legend_name = NULL) {
  ggplot(cells, aes(UMAP_1, UMAP_2, colour = .data[[colour]])) +
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
p2a <- umap_panel("dataset_label", DATASET_COLORS, "Cells by dataset") +   # legends have no titles
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
p2b <- umap_panel("condition", CONDITION_COLORS, "Cells by disease group") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
p2c <- umap_panel("cell_type", CT_COLORS, "Cells by cell type") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), ncol = 1))
legend_of <- function(p)  # legend left-aligned and vertically centred in its cell, right of its panel
  ggdraw(get_plot_component(p + theme(legend.justification.bottom = c(0, 0.5)), "guide-box-bottom"))
leg_a <- legend_of(p2a); leg_b <- legend_of(p2b); leg_c <- legend_of(p2c)
no_legend <- theme(legend.position = "none")

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

# e: scCODA effects of the one composition model (step 05): three panels side by side, SLE vs HC, SjD vs HC
# and SLE vs SjD, sharing a y-axis of cell types (order of the other figures; the reference cell type marked), scCODA's
# log2 fold change on the x-axis (the default of pertpy's effects bar plot, as the sc-best-practices chapter shows it;
# breaks -1, 0, 1; an axis line under each panel) with a line at 0; a point per credible effect at FDR 0.2 only (the
# log2FC of a non-credible cell type is only the compositional shift relative to the reference). SLE vs SjD comes from the
# fit with base level SLE (SjD_vs_SLE) with its sign reversed, so that positive means higher in SLE.
COMP_CONTRASTS <- list("SLE vs HC" = c("SLE_vs_control", 1), "SjD vs HC" = c("SjD_vs_control", 1), "SLE vs SjD" = c("SjD_vs_SLE", -1))
# points coloured by sign, the same two colours in every panel (those of Fig 4's NES scale); in SLE vs SjD positive means
# higher in SLE
COMP_COLORS <- c("log2FC > 0" = "#B2182B", "log2FC < 0" = "#2166AC")
reference <- readLines(file.path(RUN, "05_composition", "reference_cell_type.txt"), n = 1)
comp <- bind_rows(lapply(names(COMP_CONTRASTS), function(k) {
  dir <- file.path(RUN, "05_composition", COMP_CONTRASTS[[k]][1])
  cov <- sub("_vs_", "|", COMP_CONTRASTS[[k]][1]); lev <- sub("\\|.*", "", cov); base <- sub(".*\\|", "", cov)
  eff <- read_ct(file.path(dir, "effects.tsv"))
  cred <- read_ct(file.path(dir, "credible_effects.tsv"))[Covariate == sprintf("C(condition, Treatment('%s'))[T.%s]", base, lev)]
  stopifnot(nrow(cred) == nrow(eff))
  eff %>% transmute(cell_type, log2fc = as.numeric(COMP_CONTRASTS[[k]][2]) * `log2-fold change`) %>%
    left_join(cred %>% transmute(cell_type = `Cell Type`, credible = as.logical(credible)), by = "cell_type") %>%
    mutate(contrast = k)
}))
# Fig 2e draws the effects credible with every reference cell type of the sensitivity analysis (Methods):
# 05b_reference_sensitivity/reference_sensitivity.tsv (contrast, cell_type, all_references yes/no; step 05b --summary).
# Without it (a run whose sweep is not done) the credible effects of the cDC2 fit are drawn and the log says so.
sens_file <- file.path(RUN, "05b_reference_sensitivity", "reference_sensitivity.tsv")
if (file.exists(sens_file)) {
  sens <- read_ct(sens_file)[, .(contrast, cell_type, all_references = all_references == "yes")]
  comp <- comp %>% left_join(sens, by = c("contrast", "cell_type")) %>%
    mutate(credible = credible & !is.na(all_references) & all_references)
  log_msg("Fig 2e: effects credible with every reference cell type (", sens_file, ")")
} else {
  log_msg("Fig 2e: no step 05b reference_sensitivity.tsv; credible effects of the cDC2 fit only")
}
comp <- comp %>% mutate(cell_type = factor(cell_type, levels = rev(CT_FIG_ORDER)),
                        contrast = factor(contrast, levels = names(COMP_CONTRASTS)))
log_msg("Fig 2e credible effects:", paste(names(table(comp$contrast[comp$credible])), table(comp$contrast[comp$credible]), collapse = ", "))
ct_axis_label <- function(x) ifelse(x == reference, paste0(x, " (reference)"), x)
p2e <- ggplot(filter(comp, credible) %>% mutate(sign = factor(ifelse(log2fc > 0, "log2FC > 0", "log2FC < 0"), levels = names(COMP_COLORS))),
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
UMAP_MM <- 46; ROW2_MM <- 66; LEG_W_MM <- c(a = 13, b = 7, c = 25)
UMAP_W <- (FIG_WIDTH * 25.4 - sum(LEG_W_MM)) / 3
row1 <- plot_grid(p2a + no_legend, leg_a, p2b + no_legend, leg_b, p2c + no_legend, leg_c, ncol = 6,
                  rel_widths = c(UMAP_W, LEG_W_MM["a"], UMAP_W, LEG_W_MM["b"], UMAP_W, LEG_W_MM["c"]),
                  labels = c("a", "", "b", "", "c", ""), label_size = LABEL_SIZE, label_fontface = "bold")
row2 <- plot_grid(p2d, p2e, ncol = 2, rel_widths = c(0.40, 0.60), labels = c("d", "e"),
                  label_size = LABEL_SIZE, label_fontface = "bold")
save_fig(plot_grid(row1, row2, ncol = 1, rel_heights = c(UMAP_MM, ROW2_MM)), "fig2_atlas.pdf", UMAP_MM + ROW2_MM)

# =============================================================================
# Fig 3: differential expression (a DEG counts, b SLE vs SjD for the genes of a)
# =============================================================================
deg_counts <- deg[, .N, by = .(cell_type, disease, direction = ifelse(log_fc > 0, "Up", "Down"))]
indirect <- read_ct(cmd = paste("zcat", shQuote(file.path(DATA, "indirect.tsv.gz"))))
# 3b: genes changed vs HC in either disease (the DEG rule of panel a: adjusted p < 0.05 and |log2FC| >= 1), classified
# per cell type. d = log2FC_SLE - log2FC_SjD, SE = sqrt(SE_SLE^2 + SE_SjD^2) (Bucher). One multiple-testing rule, as in
# the rest of the paper: both tests are adjusted with Benjamini-Hochberg within each cell type over all genes tested in
# both diseases, and only then restricted to the genes drawn. Different: Wald test of the indirect comparison (padj of
# 07_figure_data.py), adjusted p < 0.05 and |d| >= 1, split by sign. Equivalent: two one-sided tests (Schuirmann 1987;
# Lakens et al. 2018), margin +-1 log2FC, p = the larger of the two one-sided p, adjusted p < 0.05. Inconclusive
# (Lakens et al. 2018): neither.
EQ_MARGIN <- 1
FIG3B_CATS <- c("Equivalent" = "#CC79A7", "SLE > SjD" = "#D55E00", "SjD > SLE" = "#0072B2", "Inconclusive" = "grey80")  # Equivalent: the colour of Fig 4b's header
indirect[, p_tost := pmax(pnorm((diff + EQ_MARGIN) / se, lower.tail = FALSE), pnorm((diff - EQ_MARGIN) / se))]
indirect[, padj_tost := p.adjust(p_tost, "BH"), by = cell_type]
changed <- indirect[(adj_p_value_sle < 0.05 & abs(log_fc_sle) >= 1) | (adj_p_value_sjd < 0.05 & abs(log_fc_sjd) >= 1)]
changed[, category := fifelse(padj < 0.05 & diff >= 1, "SLE > SjD", fifelse(padj < 0.05 & diff <= -1, "SjD > SLE",
                      fifelse(padj_tost < 0.05, "Equivalent", "Inconclusive")))]
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
                   label_size = LABEL_SIZE, label_fontface = "bold"), "fig3_de.pdf", 150)

# =============================================================================
# Fig 4: pathways (a NES vs healthy controls, SLE and SjD; b fgsea on the indirect comparison)
# heatmap: white cell borders, 7 mm blocks, 45-degree column labels, bodies sized to fill the page width. Column labels at
# 5.5 pt and narrower gaps (3 and 5 mm, J: 3 and 8): at the 175 mm page width a column is about 2.8 mm, and 7 pt labels at
# 45 degrees need about 3.5 mm, so they overlapped. The sign convention of b is in the caption only.
# =============================================================================
gsea <- bind_rows(lapply(c("SLE", "SjD"), function(d)
  read_ct(file.path(RUN, "06_pathways", paste0(d, "_vs_control"), "gsea.tsv"))[, block := d]))
gsea <- bind_rows(gsea, read_ct(file.path(DATA, "indirect_gsea.tsv"))[, block := "SLE vs SjD"])
gsea <- as.data.table(gsea); gsea[, pathway := format_pathway_label(gene_set)]
cts <- intersect(CT_FIG_ORDER, unique(gsea$cell_type))
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
f4 <- file.path(OUT, "fig4_pathways.pdf")
pdf_device(f4, width = w_mm / 25.4, height = h_mm / 25.4)  # figure_style.R
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
log_msg("done")
