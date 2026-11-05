#!/usr/bin/env Rscript
# Supplementary tables and figures that do not depend on the GSE253568 intersect run. Not a pipeline step
# and not a notebook port; style from figure_style.R (as 08_figures.R).
# Usage: run by make_supplementary.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration; FIG_DRAFT=1 adds
#        a DRAFT mark)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/: 04_differential_expression/, 04b_percohort_de/meta.tsv.gz, 05_composition/,
#        05b_reference_sensitivity/reference_sensitivity.tsv,
#        06_pathways/, 07_figures/data/indirect_gsea.tsv, 09_supplementary/data/ (composition_counts.tsv, qc_counts.tsv, atlas_cells.tsv.gz, atlas_dotplot.tsv)
# Output, overwritten on each run, in manuscript/supplementary/ ($SUPP_OUT); one Excel workbook per table, flat values
# with a title row:
#   supp_table_composition.xlsx  scCODA joint model, the three contrasts, with the references credible with
#   supp_table_de_sle.xlsx       PyDESeq2, SLE vs HC, every gene tested
#   supp_table_de_sjd.xlsx       PyDESeq2, SjD vs HC, every gene tested
#   supp_table_sle_vs_sjd.xlsx   SLE vs SjD indirect comparison per gene and cell type, with the Fig 3b classes
#   supp_table_gsea.xlsx         fgsea, Hallmark 2025.1, the three contrasts
#   (numbers from inputs/supplementary_numbering.tsv, in the workbook titles)
#   supp_table_qc.xlsx           cells per dataset after each quality-control filter
#   supp_fig_composition_cohorts.pdf  per-cohort composition heatmap
#   supp_fig_de_cohorts.pdf           per-cohort DE and the fixed-effect meta-analysis vs the joint model
#   supp_fig_atlas.pdf                UMAPs before integration and of the controls, annotation clusters, marker dot plot
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot); library(grid); library(data.table)
})

ATLAS_RUN <- Sys.getenv("ATLAS_RUN", "02_integration")
ATLAS_ROOT <- Sys.getenv("ATLAS_ROOT", ".")  # set by run.sh
RUN  <- file.path(Sys.getenv("PIPELINE_OUT", file.path(ATLAS_ROOT, "out")), "py", "atlas", ATLAS_RUN)
OUT <- Sys.getenv("SUPP_OUT", file.path(ATLAS_ROOT, "manuscript", "supplementary"))
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", file.path(ATLAS_ROOT, "pipeline_py"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DRAFT <- Sys.getenv("FIG_DRAFT", "0") == "1"
DRAFT_TEXT <- paste0("DRAFT | ", ATLAS_RUN, " | ", format(Sys.Date()))
log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
pdf(NULL)  # null default device: plot measurements otherwise write Rplots.pdf into the working directory
source(file.path(PIPELINE_DIR, "figure_style.R"))

# =============================================================================
# Tables: one flat workbook each, title in row 1 (and a DRAFT mark in row 2 with FIG_DRAFT=1), header, data, notes
# =============================================================================
# Supplementary numbers, in order of first citation in the paper: inputs/supplementary_numbering.tsv; each workbook
# title carries its number
SUPP_NUM <- fread(file.path(PIPELINE_DIR, "inputs", "supplementary_numbering.tsv"))
write_table <- function(df, file, title, notes = character(), int_cols = character(), dec_cols = character()) {
  num <- SUPP_NUM$number[SUPP_NUM$kind == "Table" & SUPP_NUM$file == file]
  stopifnot(length(num) == 1, startsWith(title, "Supplementary Table. "))
  title <- sub("^Supplementary Table\\. ", sprintf("Supplementary Table %d. ", num), title)
  r0 <- if (DRAFT) 3 else 2
  wb <- openxlsx2::wb_workbook()
  wb <- openxlsx2::wb_add_worksheet(wb, "Table")
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
  wb <- openxlsx2::wb_set_col_widths(wb, cols = seq_len(ncol(df)), widths = "auto")
  openxlsx2::wb_save(wb, file.path(OUT, file), overwrite = TRUE)
  log_msg("saved", file.path(OUT, file), "|", nrow(df), "rows")
}
ct_rank <- function(x) match(x, CT_ORDER)

# --- composition: scCODA joint model, three contrasts. SLE vs SjD is the fit with base level SLE (SjD vs SLE) with
# every estimate negated and the interval bounds swapped, so that positive means higher in SLE, as in Fig 2e. ---
COMP <- list("SLE vs HC" = c("SLE_vs_control", "control", "SLE", 1), "SjD vs HC" = c("SjD_vs_control", "control", "SjD", 1),
             "SLE vs SjD" = c("SjD_vs_SLE", "SLE", "SjD", -1))
reference <- readLines(file.path(RUN, "05_composition", "reference_cell_type.txt"))[1]
comp <- bind_rows(lapply(names(COMP), function(k) {
  x <- COMP[[k]]; s <- as.numeric(x[4]); dir <- file.path(RUN, "05_composition", x[1])
  ef <- read_ct(file.path(dir, "effects.tsv"))
  cr <- read_ct(file.path(dir, "credible_effects.tsv"))[Covariate == sprintf("C(condition, Treatment('%s'))[T.%s]", x[2], x[3])]
  lo <- ef[["ETI 89% lower"]]; hi <- ef[["ETI 89% upper"]]
  data.frame(Contrast = k, `Cell type` = ef$cell_type,
             `Reference cell type` = ifelse(ef$cell_type == reference, "yes", ""),
             `Final parameter` = s * ef[["Final Parameter"]],
             `ETI 89% lower` = if (s > 0) lo else -hi, `ETI 89% upper` = if (s > 0) hi else -lo,
             SD = ef$SD, `Inclusion probability` = ef[["Inclusion probability"]],
             `log2 fold change` = s * ef[["log2-fold change"]],
             `Credible (FDR 0.2)` = ifelse(as.logical(cr$credible[match(ef$cell_type, cr[["Cell Type"]])]), "yes", "no"),
             check.names = FALSE)
})) %>% arrange(match(Contrast, names(COMP)), ct_rank(`Cell type`))
# references with which each effect is credible (step 05b, the reference sensitivity)
sens <- read_ct(file.path(RUN, "05b_reference_sensitivity", "reference_sensitivity.tsv"))
comp <- comp %>% left_join(sens[, .(Contrast = contrast, `Cell type` = cell_type, `Credible with references` = credible_with)],
                           by = c("Contrast", "Cell type")) %>%
  mutate(`Credible with references` = ifelse(`Reference cell type` == "yes", "", `Credible with references`))
write_table(comp, "supp_table_composition.xlsx",
            "Supplementary Table. Cell-type composition",
            c("SLE vs HC and SjD vs HC: fit with healthy controls (HC) as the base level; SLE vs SjD: fit with SLE as the base level, estimates negated so that positive means higher in SLE.",
              "Final parameter: scCODA effect on the log scale, 0 when the effect is not included; ETI: equal-tailed interval; credible: inclusion at a false discovery rate of 0.2.",
              "log2 fold change: scCODA's log2 fold change of the expected proportion, which includes the shift relative to the reference cell type.",
              "Credible with references: references (cDC2, regulatory T, non-classical monocytes, NK, CD8 effector memory T, pDC; the cell type itself excluded) with which the effect is credible; effects credible with every reference are reported."),
            dec_cols = c("Final parameter", "ETI 89% lower", "ETI 89% upper", "SD", "Inclusion probability", "log2 fold change"))

# --- differential expression: every gene tested, per disease ---
for (d in c("SLE", "SjD")) {
  x <- read_ct(file.path(RUN, "04_differential_expression", paste0(d, "_vs_control"), "all", "de_all_cell_types.tsv"))
  t <- data.frame(`Cell type` = x$cell_type, Gene = x$variable, `Base mean` = x$baseMean, log2FC = x$log_fc,
                  `log2FC SE` = x$lfcSE, `Wald statistic` = x$stat, P = x$p_value, `Adjusted P` = x$adj_p_value,
                  DEG = ifelse(!is.na(x$adj_p_value) & x$adj_p_value < 0.05 & abs(x$log_fc) >= 1, "yes", "no"),
                  Cases = x$n_case, Controls = x$n_control, Design = x$design, check.names = FALSE) %>%
    arrange(ct_rank(`Cell type`), P)
  write_table(t, sprintf("supp_table_de_%s.xlsx", tolower(d)),
              sprintf("Supplementary Table. Differential expression, %s vs healthy controls", d),
              c("Every gene tested in each cell type, after the edgeR filterByExpr filter and the filter keeping genes expressed in at least 10% of the cells in at least two pseudobulk samples; Adjusted P: Benjamini-Hochberg within cell type after independent filtering (empty when the gene was filtered out).",
                "DEG: adjusted P < 0.05 and |log2FC| >= 1. Cases and Controls: donors with a pseudobulk sample in the cell type. Design: the model of the cell type, terms with one level dropped."),
              int_cols = c("Cases", "Controls"))
}

# --- SLE vs SjD: the indirect comparison per gene and cell type (Fig 3b; Methods), every gene tested in both diseases;
# equivalence (TOST) and the classes of the genes changed vs HC as in 08_figures.R ---
ind <- read_ct(cmd = paste("zcat", shQuote(file.path(RUN, "07_figures", "data", "indirect.tsv.gz"))))
EQ_MARGIN <- 1  # as 08_figures.R
ind[, p_tost := pmax(pnorm((diff + EQ_MARGIN) / se, lower.tail = FALSE), pnorm((diff - EQ_MARGIN) / se))]
ind[, padj_tost := p.adjust(p_tost, "BH"), by = cell_type]
ind[, `:=`(deg_sle = !is.na(adj_p_value_sle) & adj_p_value_sle < 0.05 & abs(log_fc_sle) >= 1,
           deg_sjd = !is.na(adj_p_value_sjd) & adj_p_value_sjd < 0.05 & abs(log_fc_sjd) >= 1)]
ind[, class := fifelse(!(deg_sle | deg_sjd), "", fifelse(padj < 0.05 & diff >= 1, "Higher in SLE",
               fifelse(padj < 0.05 & diff <= -1, "Higher in SjD", fifelse(padj_tost < 0.05, "Equivalent", "Inconclusive"))))]
stopifnot(!anyNA(ind$class))
it <- data.frame(`Cell type` = ind$cell_type, Gene = ind$variable, `log2FC SLE` = ind$log_fc_sle, `log2FC SE SLE` = ind$lfcSE_sle,
                 `log2FC SjD` = ind$log_fc_sjd, `log2FC SE SjD` = ind$lfcSE_sjd, `Difference (SLE - SjD)` = ind$diff,
                 `Difference SE` = ind$se, z = ind$z, P = ind$pvalue, `Adjusted P` = ind$padj, `Equivalence P` = ind$p_tost,
                 `Adjusted equivalence P` = ind$padj_tost, `DEG in SLE` = ifelse(ind$deg_sle, "yes", "no"),
                 `DEG in SjD` = ifelse(ind$deg_sjd, "yes", "no"), Class = ind$class, check.names = FALSE) %>%
  arrange(ct_rank(`Cell type`), P)
write_table(it, "supp_table_sle_vs_sjd.xlsx",
            "Supplementary Table. SLE vs SjD comparison of gene expression",
            c("Every gene tested in both diseases within a cell type. Difference: log2FC SLE - log2FC SjD, each against the healthy controls of its own datasets; Difference SE: square root of the sum of the two squared SE; z: Difference / Difference SE; P: two-sided, standard normal; Adjusted P: Benjamini-Hochberg within cell type.",
              "Equivalence P: two one-sided tests with an equivalence margin of 1 log2FC (|Difference| < 1), the larger of the two one-sided P; Adjusted equivalence P: Benjamini-Hochberg within cell type.",
              "DEG in SLE, DEG in SjD: adjusted P < 0.05 and |log2FC| >= 1 against healthy controls. Class, for genes differentially expressed in either disease: Higher in SLE or Higher in SjD, adjusted P < 0.05 and |Difference| >= 1; Equivalent, adjusted equivalence P < 0.05; Inconclusive, neither."),
            dec_cols = c("log2FC SLE", "log2FC SE SLE", "log2FC SjD", "log2FC SE SjD", "Difference (SLE - SjD)", "Difference SE", "z"))
log_msg("SLE vs SjD table:", nrow(it), "rows;", paste(names(table(it$Class[it$Class != ""])), table(it$Class[it$Class != ""]), collapse = ", "))

# --- GSEA: three contrasts ---
gs <- bind_rows(read_ct(file.path(RUN, "06_pathways", "SLE_vs_control", "gsea.tsv"))[, Contrast := "SLE vs HC"],
                read_ct(file.path(RUN, "06_pathways", "SjD_vs_control", "gsea.tsv"))[, Contrast := "SjD vs HC"],
                read_ct(file.path(RUN, "07_figures", "data", "indirect_gsea.tsv"))[, Contrast := "SLE vs SjD"])
gt <- data.frame(Contrast = gs$Contrast, `Cell type` = gs$cell_type, `Gene set` = gs$gene_set, NES = gs$score,
                 P = gs$pvalue, `Adjusted P` = gs$padj, check.names = FALSE) %>%
  arrange(match(Contrast, names(COMP)), ct_rank(`Cell type`), P)
write_table(gt, "supp_table_gsea.xlsx",
            "Supplementary Table. Gene set enrichment",
            c("SLE vs HC and SjD vs HC: genes ranked by the PyDESeq2 Wald statistic; SLE vs SjD: by the z statistic of the indirect comparison (positive NES: higher in SLE).",
              "Adjusted P: Benjamini-Hochberg across all cell type and gene set tests of a contrast."),
            dec_cols = "NES")

# --- QC per dataset: the recorded counts of step 01, the atlas and the annotation exclusion (09_supplementary_data.py);
# the per-sample thresholds are not recorded, so the table has none. Doublet rate: of the cells that entered doublet
# detection; annotation exclusion: of the atlas cells. ---
qc <- read_ct(file.path(RUN, "09_supplementary", "data", "qc_counts.tsv"))
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
n_atlas <- sum(qc$final)
write_table(qt, "supp_table_qc.xlsx",
            "Supplementary Table. Quality control: cells per dataset after each filter",
            c(paste0("Outlier filter: cells beyond 5 median absolute deviations (MADs) of the sample median in log total counts, log detected genes or the share of counts in the 20 most expressed genes. ",
                     "Mitochondrial filter: cells above 3 MADs of the sample median in the logit mitochondrial fraction and with more than 8% mitochondrial counts. ",
                     "Doublets: scDblFinder per library; doublet rate of the cells that passed both filters. Excluded at annotation: clusters labelled low quality or platelets."),
              "HRA003613: doublet detection not run; the authors removed doublets before deposit.",
              sprintf("GSE174188: %d samples passed quality control; the %d later-visit samples of donors with more than one visit were removed at integration, leaving the first visit of each donor.",
                      qc[dataset == "Perez", samples_qc], qc[dataset == "Perez", samples_qc - samples]),
              "HRA003613: 1 cell with no counts in the genes of the atlas was removed at integration.",
              "GSE316482: 1 cell empty after ambient RNA correction was removed before doublet detection."),
            int_cols = setdiff(names(qt), c("Dataset", "Doublet rate (%)", "Excluded at annotation (%)")))
log_msg("QC table: final cells", n_atlas, "| samples", sum(qc$samples))

# =============================================================================
# Figures
# =============================================================================
draft_mark <- function(p) {
  if (!DRAFT) return(p)
  ggdraw(p) + draw_label(DRAFT_TEXT, x = 0.995, y = 0.998, hjust = 1, vjust = 1, colour = "red",
                         fontface = "bold", size = MIN_PT, fontfamily = FONT_FAMILY)
}
PREVIEW <- file.path(RUN, "09_supplementary")
save_supp <- function(p, name, height_mm) {
  f <- file.path(OUT, name)
  ggsave(f, draft_mark(p), width = FIG_WIDTH, height = height_mm / 25.4, device = pdf_device, bg = "white")  # figure_style.R
  log_msg("saved", f)
}
SLE_DS <- datasets$dataset[datasets$disease == "SLE"]; SJD_DS <- datasets$dataset[datasets$disease == "SjD"]

# --- per-cohort composition: per dataset and cell type, log2 of the mean proportion in its cases over the mean
# proportion in its own controls (donor proportions); a point marks an effect of the joint scCODA model for that
# dataset's disease (SLE vs HC or SjD vs HC) credible with every reference cell type (step 05b), the rule of the text and
# Fig 2e ---
cc <- read_ct(file.path(RUN, "09_supplementary", "data", "composition_counts.tsv"))
cc[, prop := n / sum(n), by = donor_id]
pc <- cc[, .(mean_prop = mean(prop), donors = uniqueN(donor_id)), by = .(dataset, condition, cell_type)]
pc <- merge(pc[condition != "control"], pc[condition == "control"], by = c("dataset", "cell_type"), suffixes = c("", "_hc"))
pc[, log2_ratio := log2(mean_prop / mean_prop_hc)]
cred <- read_ct(file.path(RUN, "05b_reference_sensitivity", "reference_sensitivity.tsv"))[contrast %in% c("SLE vs HC", "SjD vs HC")]
cred[, `:=`(disease = sub(" vs HC", "", contrast), credible = all_references == "yes")]
pc[, disease := condition]
pc <- merge(pc, cred[, .(disease, cell_type, credible)], by = c("disease", "cell_type"), all.x = TRUE)
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
save_supp(pcomp, "supp_fig_composition_cohorts.pdf", 95)
log_msg("per-cohort composition: points mark joint-model effects credible with every reference:",
        paste(unique(pc[credible == TRUE, paste(disease, cell_type)]), collapse = "; "), "| log2 ratios capped at +-", lim)

# --- per-cohort DE: a, per-cohort log2FC compared between cohorts of the same disease, over the DEGs of the joint model
# tested in both cohorts; b, fixed-effect meta-analysis vs the joint model, every gene tested in both (joint DEGs in
# colour); c, I^2 of the joint DEGs tested in more than one cohort ---
meta <- read_ct(cmd = paste("zcat", shQuote(file.path(RUN, "04b_percohort_de", "meta.tsv.gz"))))
meta[, joint_deg := !is.na(joint_adj_p_value) & joint_adj_p_value < 0.05 & abs(joint_log_fc) >= 1]
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
save_supp(plot_grid(row_a, NULL, row_b, ncol = 1, rel_heights = c(1, 0.08, 1.25)), "supp_fig_de_cohorts.pdf", 130)
for (d in c("SLE", "SjD")) {
  x <- meta[disease == d & joint_deg]
  log_msg(d, "joint DEGs:", nrow(x), "| tested in > 1 cohort:", sum(x$k > 1), "| I2 > 0.5:", sum(x$I2 > 0.5, na.rm = TRUE),
          "| same sign in the meta-analysis:", sum(sign(x$meta_log_fc) == sign(x$joint_log_fc)))
}

# --- atlas: a, UMAP before integration by dataset; b, integrated UMAP of the controls by dataset, an equal
# random sample of control cells per dataset (the smallest dataset's number); c, integrated UMAP by annotation cluster
# (numbered in the cell-type order, excluded clusters last and black); a-c in one row as Fig 2a-c, the dataset legend in
# one row under a and b, the cluster key right of c; d, dot plot of config.SUPP_DOTPLOT_MARKERS by final cell type (09_supplementary_data.py),
# full width below ---
ac <- read_ct(cmd = paste("zcat", shQuote(file.path(RUN, "09_supplementary", "data", "atlas_cells.tsv.gz"))))
setnames(ac, "label", "cell_type"); ac[, cell_type := ct_display(cell_type)]
ac[, dataset_label := factor(DATASET_LABEL[dataset], levels = DATASET_LABEL)]
set.seed(42); ac <- ac[sample(.N)]  # random drawing order
ctrl <- ac[condition == "control"]
n_b <- min(ctrl[, .N, by = dataset]$N)
set.seed(42); ctrl <- ctrl[ctrl[, .I[sample(.N, n_b)], by = dataset]$V1][sample(.N)]  # n_b cells per dataset, random order
EXCL <- c("low quality", "platelets")
cl <- unique(ac[, .(cluster, cell_type)])
stopifnot(!anyDuplicated(cl$cluster), all(cl$cell_type %in% c(CT_ORDER, EXCL)))
cl[, res := as.integer(sub("leiden_(\\d):.*", "\\1", cluster))][, id := as.integer(sub(".*:", "", cluster))]
setorder(cl[, rank := match(cell_type, c(CT_ORDER, EXCL))], rank, res, id)
cl[, number := seq_len(.N)]
ac[, number := factor(cl$number[match(cluster, cl$cluster)], levels = cl$number)]
kept <- cl[!cell_type %in% EXCL]
CL_COLORS <- c(setNames(scales::hue_pal()(nrow(kept)), kept$number), setNames(rep("black", nrow(cl) - nrow(kept)), cl[cell_type %in% EXCL, number]))
# a-c: the UMAP panels of Fig 2 (08_figures.R umap_panel, copied; the data and coordinates as arguments), legends at 5 pt
# (LEG_PT, as 08_figures.R): the dataset legend in one row under a and b, serving both; the cluster key right of c
LEG_PT <- 5
umap_supp <- function(d, colour, values, title, axes, legend_name = NULL) {  # CHANGED from 08_figures.R umap_panel: d, axes
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
pa <- umap_supp(ac, "dataset_label", DATASET_COLORS, "Before integration", "umap_unint")
pb <- umap_supp(ctrl, "dataset_label", DATASET_COLORS, "Only HC donors, integrated", "umap") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), nrow = 1))  # one row, centred under a and b
cen <- ac[, .(x = median(umap_1), y = median(umap_2)), by = number]
pc_ <- umap_supp(ac, "number", CL_COLORS, "Annotation clusters", "umap") +
  ggrepel::geom_text_repel(data = cen, aes(x, y, label = number), inherit.aes = FALSE, size = LEG_PT / .pt, family = FONT_FAMILY,
                           fontface = "bold", seed = 42, box.padding = 0.08, point.padding = 0, min.segment.length = 0.1,
                           segment.size = 0.2, max.overlaps = Inf, bg.colour = "white", bg.r = 0.12)  # Fig 2 labels no clusters: ggrepel, short segments; a white outline keeps numbers on the black clusters readable
leg_ab <- ggdraw(get_plot_component(pb, "guide-box-bottom"))  # centred under a and b
no_legend <- theme(legend.position = "none")
# cluster key right of c: one entry per line at the line pitch of the legends (key 0.2 cm + 0.5 pt), then a blank line
# and the note, vertically centred
key <- cl[, .(n = paste(number, collapse = ", ")), by = cell_type][
  , txt := sprintf("%s: %s", ifelse(cell_type %in% EXCL, sub("^(.)", "\\U\\1", cell_type, perl = TRUE), cell_type), n)]  # "Low quality", "Platelets"
key_lines <- c(key$txt, "", "Black: excluded at annotation")
UMAP_MM <- 46  # as Fig 2
pitch <- (2 + 0.5 / 72.27 * 25.4) / UMAP_MM  # line pitch as a fraction of the row height
key_y <- 0.5 + (rev(seq_along(key_lines)) - (length(key_lines) + 1) / 2) * pitch
key_c <- ggdraw() + draw_text(key_lines, x = 0, y = key_y, hjust = 0, size = LEG_PT, family = FONT_FAMILY)
# widest line, measured on the figures' device (the null device's metrics are narrower), plus 3 mm to the page edge
pdf_device(tf <- tempfile(fileext = ".pdf"), 7, 7)
key_w_mm <- max(sapply(key_lines, function(x) convertWidth(grobWidth(textGrob(x, gp = gpar(fontsize = LEG_PT, fontfamily = FONT_FAMILY))), "mm", valueOnly = TRUE))) + 3
invisible(dev.off()); unlink(tf)
dot <- read_ct(file.path(RUN, "09_supplementary", "data", "atlas_dotplot.tsv"))
dot[, cell_type := ct_display(cell_type)]
mk <- unique(dot$gene)  # in the order of config.SUPP_DOTPLOT_MARKERS (09_supplementary_data.py writes them so)
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
save_supp(plot_grid(row1, plot_grid(pd_, labels = "d", label_size = LABEL_SIZE, label_fontface = "bold"),
                    ncol = 1, rel_heights = c(UMAP_MM + LEG_ROW_MM, DOT_MM)), "supp_fig_atlas.pdf", UMAP_MM + LEG_ROW_MM + DOT_MM)
log_msg("atlas figure: panel b", n_b, "control cells per dataset; key", round(key_w_mm, 1), "mm wide; UMAP cells", round(UMAP_W, 1), "mm")
log_msg("atlas figure:", nrow(ac), "cells,", nrow(ctrl), "control cells drawn,", nrow(cl), "clusters (",
        sum(cl$cell_type %in% EXCL), "excluded ),", length(mk), "markers:",
        paste(unique(dot[n_datasets < n_datasets_all, as.character(label)]), collapse = ", "), "outside the shared genes")
log_msg("done")
