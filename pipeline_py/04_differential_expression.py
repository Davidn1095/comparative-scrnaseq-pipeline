# Differential gene expression between conditions: pseudobulk per donor and cell type, PyDESeq2 per cell type.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19):
#   jupyter-book/conditions/differential_gene_expression.ipynb, setup cells 9-10 and cells 26-66
#   (pseudobulking, variability exploration, feature selection, PyDESeq2 for one cell type). In feature selection (cells
#   49-51) decoupler's filter_by_expr is replaced by edgeR::filterByExpr with the model's design.
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason. Plots
# are written to files instead of displayed.
# Two contrasts, each on its own datasets: SLE vs controls (SLE datasets) and SjD vs controls (SjD datasets). Design
# ~ batch + sex + condition (inferred sex; batch = dataset, or Perez's processing cohort), batch dropped when only one is
# present. All donors only (the
# women-only check was dropped on 2026-11-02).
# Usage: python 04_differential_expression.py   (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/atlas.h5ad (layer "counts"; obs donor_id, dataset, condition, inferred_sex)
#        and 03_annotation/cells_manual_annotation.tsv.gz (manual cell-type labels)
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/04_differential_expression/<contrast>/all/
#   de_<cell type>.tsv (PyDESeq2 results), de_all_cell_types.tsv, design.txt, figures/*.png

# --- cell 9 ---
import warnings

warnings.filterwarnings("ignore")

# --- cell 10 ---
import os  # CHANGED: run folder

import decoupler as dc
import matplotlib  # CHANGED: plots are written to files

matplotlib.use("Agg")  # CHANGED: as above
import matplotlib.pyplot as plt  # CHANGED: as above
import numpy as np
import pandas as pd
import pertpy as pt
import scanpy as sc

from config import DATASETS, OUT_ROOT, batch_of, excluded_labels, log  # CHANGED: pipeline paths; labels left out; batch
import rpy2.robjects as ro  # CHANGED: edgeR::filterByExpr with the design (cell 51)
from rpy2.robjects import numpy2ri, pandas2ri  # CHANGED: as above
from rpy2.robjects.conversion import localconverter  # CHANGED: as above

ro.r("suppressMessages(library(edgeR))")  # CHANGED: as above

# CHANGED: lamindb (ln.connect, ln.track) removed

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline input
out_root = run / "04_differential_expression"  # CHANGED: pipeline output


def save(fig_dir, name, fig=None):  # CHANGED: plots are written to files instead of displayed
    fig_dir.mkdir(parents=True, exist_ok=True)
    (fig or plt.gcf()).savefig(fig_dir / name, dpi=150, bbox_inches="tight")
    plt.close("all")


def safe(name):  # CHANGED: cell-type names as file names
    return "".join(c if c.isalnum() or c in "-." else "_" for c in str(name))


# --- cells 14-23 (data loading and cell/gene filtering of the notebook's dataset) replaced ---
# CHANGED: the input is the atlas of 02_integration.py, whose layer "counts" holds the ambient-corrected counts, with the
# manual cell-type labels of 03_annotation.py; cells without a label are not used. QC and gene filtering were done in
# 01_quality_control.py and 02_integration.py.
adata = sc.read_h5ad(run / "atlas.h5ad")  # CHANGED: as above
for key in list(adata.layers):  # CHANGED: memory; only the counts are used
    if key != "counts":
        del adata.layers[key]
adata.X = None  # CHANGED: as above
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)  # CHANGED: as above
adata.obs["cell_type"] = labels["manual_celltype_annotation"].reindex(adata.obs_names)  # CHANGED: as above
adata = adata[adata.obs["cell_type"].notna() & ~excluded_labels(adata.obs)].copy()  # CHANGED: as above; cells labelled low quality or platelets are left out
adata.obs["sex"] = adata.obs["inferred_sex"].astype(str)  # CHANGED: inferred sex as the design covariate
adata.obs["batch"] = batch_of(adata.obs).values  # CHANGED: batch = dataset, or Perez's processing cohort (Perez_1-4), for the design
log(f"{adata.n_obs} labelled cells, {adata.obs['cell_type'].nunique()} cell types")  # CHANGED: log

contrasts = {  # CHANGED: two contrasts, each on the datasets of its disease
    f"{disease}_vs_control": (disease, DATASETS.loc[DATASETS.disease == disease, "dataset"].tolist())
    for disease in ["SLE", "SjD"]
}

for contrast, (case, datasets) in contrasts.items():  # CHANGED: see above
    for subset in ["all"]:  # CHANGED: all donors, written to all/ (the women-only check was dropped)
        sub = adata[adata.obs["dataset"].isin(datasets) & adata.obs["condition"].isin([case, "control"])]  # CHANGED: as above
        sub = sub.copy()  # CHANGED: as above
        out = out_root / contrast / subset  # CHANGED: pipeline output
        fig = out / "figures"  # CHANGED: as above
        out.mkdir(parents=True, exist_ok=True)  # CHANGED: as above
        log(f"{contrast} / {subset}: {sub.n_obs} cells, {sub.obs['donor_id'].nunique()} donors")  # CHANGED: log

        # --- cell 26 ---
        # CHANGED: one pseudobulk sample per donor (the notebook's sample is replicate x condition, since each of its
        # replicates was measured in both conditions; here each donor has one condition)
        sub.obs["sample"] = pd.Categorical(sub.obs["donor_id"].astype(str))

        # --- cell 28 ---
        adata_pb = dc.pp.pseudobulk(
            sub, sample_col="sample",  # CHANGED: the subset of this contrast
            groups_col="cell_type",
            layer="counts",
            mode="sum"
        )

        # --- cell 30 ---
        dc.pl.filter_samples(
            adata=adata_pb,
            groupby=["condition", "cell_type", "dataset"],  # CHANGED: this atlas's covariates (label/replicate in the notebook)
            min_cells=10,
            min_counts=1000,
            figsize=(5, 8),
        )
        save(fig, "filter_samples.png")  # CHANGED: saved to file

        # --- cell 32 ---
        dc.pp.filter_samples(adata_pb, min_cells=10,min_counts=1000)

        # --- cell 34 ---
        dc.pl.obsbar(adata=adata_pb, y="cell_type", hue="condition", figsize=(6, 3))  # CHANGED: condition (label in the notebook)
        save(fig, "obsbar_cell_type.png")  # CHANGED: saved to file

        # --- cell 38 ---
        adata_pb.layers['counts'] = adata_pb.X.copy()

        sc.pp.normalize_total(adata_pb, target_sum=1e6)
        sc.pp.log1p(adata_pb)
        sc.pp.scale(adata_pb, max_value=10)
        sc.tl.pca(adata_pb)

        dc.pp.swap_layer(adata=adata_pb, key="counts", inplace=True)

        # --- cell 40 ---
        adata_pb.obs["psbulk_counts_log"] = np.log(adata_pb.obs["psbulk_counts"])
        adata_pb.obs["psbulk_cells_log"] = np.log(adata_pb.obs["psbulk_cells"])

        # --- cell 41 ---
        dc.tl.rankby_obsm(adata_pb, key="X_pca",  # CHANGED: obs columns with at least two levels only (exploratory, no
                          obs_keys=[c for c in adata_pb.obs.columns if adata_pb.obs[c].nunique() >= 2])  # effect on the DE): within a contrast "disease" has one level, and ANOVA needs two groups
        sc.pl.pca_variance_ratio(adata_pb, show=False)  # CHANGED: saved below
        save(fig, "pca_variance_ratio.png")  # CHANGED: saved to file
        f = dc.pl.obsm(
            adata=adata_pb,
            return_fig=True,
            nvar=5,
            titles=["PC scores", "Adjusted p-values"],
            figsize=(10, 5)
        )
        save(fig, "pca_covariates.png", f)  # CHANGED: saved to file

        # --- cell 43 ---
        adata_pb.obs = adata_pb.obs.sort_index(axis=1)
        covariates = [c for c in ["condition", "dataset", "batch", "sex", "cell_type", "psbulk_counts_log", "psbulk_cells_log"] if c in adata_pb.obs]  # CHANGED: the design covariates and pseudobulk size (the notebook colours by every obs column)
        sc.pl.pca(adata_pb, color=covariates, ncols=2, size=300, show=False)  # CHANGED: as above; saved below
        save(fig, "pca.png")  # CHANGED: saved to file

        results = []  # CHANGED: one result table per cell type, also written together
        for cell_type in sorted(adata_pb.obs["cell_type"].unique()):  # CHANGED: every cell type (the notebook shows CD14+ monocytes)
            # --- cell 47 ---
            adata_ct = adata_pb[adata_pb.obs["cell_type"] == cell_type].copy()  # CHANGED: as above
            n = adata_ct.obs["condition"].value_counts()  # CHANGED: a contrast needs donors of both conditions
            if n.get(case, 0) < 2 or n.get("control", 0) < 2:  # CHANGED: as above
                log(f"  {cell_type}: skipped ({n.get(case, 0)} {case}, {n.get('control', 0)} control pseudobulk samples)")  # CHANGED: as above
                continue  # CHANGED: as above
            ct_fig = fig / safe(cell_type)  # CHANGED: per-cell-type figures

            # --- cell 49 ---
            # CHANGED: decoupler's filter_by_expr plot left out: the filter is edgeR's own filterByExpr with the design
            # (cell 51 below), which this plot does not show
            dc.pl.filter_by_prop(
                adata=adata_ct,
                min_prop=0.1,
                min_smpls=2,
            )
            save(ct_fig, "filter_by_prop.png")  # CHANGED: saved to file

            # --- cell 51 ---
            # CHANGED: decoupler's filter_by_expr(group="condition") replaced by edgeR::filterByExpr with
            # design = model.matrix(~ batch + sex + condition), the model fitted below, and the chapter's thresholds.
            # edgeR's documentation: MinSampleSize "is the smallest group sample size or, more generally, the minimum
            # inverse leverage computed from the design matrix", and group is only "an alternative way to specify the
            # design matrix when there is only one explanatory factor"; with covariates the minimum sample size comes
            # from the design, not from group. With group="condition" it was about 75 SLE samples
            terms_ct = [c for c in ["batch", "sex"] if adata_ct.obs[c].nunique() > 1] + ["condition"]  # CHANGED: the design of cell 53
            counts_ct = adata_ct.X.toarray() if hasattr(adata_ct.X, "toarray") else np.asarray(adata_ct.X)  # CHANGED: as above
            with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):  # CHANGED: as above
                ro.globalenv["y"] = counts_ct.T.astype(float)  # CHANGED: genes x samples, as edgeR expects
                ro.globalenv["md"] = adata_ct.obs[terms_ct].astype(str).reset_index(drop=True)  # CHANGED: as above
            ro.r(f"design <- model.matrix(~ {' + '.join(terms_ct)}, data = md); "  # CHANGED: as above
                 "keep <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7); "
                 "mss <- 1 / max(hat(design))")
            keep = np.asarray(ro.r("keep"), dtype=bool)  # CHANGED: as above
            log(f"  {cell_type}: edgeR MinSampleSize {float(ro.r('mss')[0]):.1f}, {int(keep.sum())} genes kept by filterByExpr")  # CHANGED: log
            adata_ct = adata_ct[:, keep].copy()  # CHANGED: as above
            dc.pp.filter_by_prop(
                adata=adata_ct,
                min_prop=0.1,
                min_smpls=2,
            )
            log(f"  {cell_type}: {adata_ct.n_vars} genes after filter_by_prop")  # CHANGED: log

            # --- cell 53 ---
            # CHANGED: design ~ batch + sex + condition (the notebook: ~ label); a covariate with one level in this cell
            # type (one batch, or one sex) is dropped, since it cannot be estimated
            terms = [c for c in ["batch", "sex"] if adata_ct.obs[c].nunique() > 1] + ["condition"]  # CHANGED: as above; batch instead of dataset: Perez's processing cohorts are unbalanced between SLE and controls
            design = "~ " + " + ".join(terms)  # CHANGED: as above
            pds2 = pt.tl.PyDESeq2(adata=adata_ct,design=design)  # CHANGED: as above

            # --- cell 54 ---
            pds2.fit()

            # --- cell 56 ---
            pds2.plot_disp_ests()
            save(ct_fig, "dispersion_estimates.png")  # CHANGED: saved to file

            # --- cell 58 ---
            res_df = pds2.test_contrasts(pds2.contrast(
                column="condition",  # CHANGED: condition (label in the notebook)
                baseline="control",  # CHANGED: as above (ctrl)
                group_to_compare=case),  # CHANGED: as above (stim)
                n_cpus=len(os.sched_getaffinity(0)),  # CHANGED: the job's CPUs, as pertpy's fit() uses; by default the Wald test starts one worker per node core (multiprocessing.cpu_count), ~2 min per cell type on a 256-core node; same results
            )

            # --- cell 59 ---
            res_df.to_csv(out / f"de_{safe(cell_type)}.tsv", sep="\t", index=False)  # CHANGED: written out instead of res_df.head(10)
            results.append(res_df.assign(cell_type=cell_type, design=design, n_case=int(n[case]), n_control=int(n["control"])))  # CHANGED: as above

            # --- cell 61 ---
            pds2.plot_volcano(res_df, log2fc_threshold=0)
            save(ct_fig, "volcano.png")  # CHANGED: saved to file

            # --- cell 64 ---
            pds2.plot_ma(res_df)
            save(ct_fig, "ma.png")  # CHANGED: saved to file

            # --- cell 66 (plot_paired by replicate) left out ---
            # CHANGED: the notebook's replicates were measured in both conditions; here each donor has one condition

            log(f"  {cell_type}: {design}; {int(n[case])} {case} vs {int(n['control'])} control; "
                f"{int((res_df['adj_p_value'] < 0.05).sum())} genes with adjusted p < 0.05")  # CHANGED: log

        if results:  # CHANGED: all cell types of the contrast together
            pd.concat(results).to_csv(out / "de_all_cell_types.tsv", sep="\t", index=False)  # CHANGED: as above
log("done")  # CHANGED: log
