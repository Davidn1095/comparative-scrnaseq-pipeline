# Per-cohort differential expression and its fixed-effect meta-analysis (supplementary).
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19), as 04_differential_expression.py:
#   jupyter-book/conditions/differential_gene_expression.ipynb, setup cells 9-10 and cells 26-58 (pseudobulking, feature
#   selection, PyDESeq2), with the changes of step 04 (edgeR::filterByExpr with the model's design in cell 51).
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason. The
# exploratory plots of step 04 (cells 30, 34, 38-43, 49, 56, 61, 64) are left out: they are step 04's diagnostics.
# Each dataset is fitted alone, each disease against that dataset's own controls, with step 04's design minus the
# terms that have one level in the dataset: ~ batch + sex + condition, where batch remains only in Perez (its processing
# cohorts Perez_1-4). Then, per disease, cell type and gene, a fixed-effect (inverse-variance) meta-analysis of the
# per-cohort log2FC with Cochran's Q and I^2, compared with step 04's joint model.
# Usage: python 04b_percohort_de.py   (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration; after step 04)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/atlas.h5ad (layer "counts"), 03_annotation/cells_manual_annotation.tsv.gz,
#        04_differential_expression/<contrast>/all/de_all_cell_types.tsv (joint model)
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/04b_percohort_de/
#   de_<dataset>.tsv.gz (PyDESeq2 results per cell type), meta.tsv.gz (per disease, cell type and gene: number of
#   cohorts k, fixed-effect log2FC and SE, Wald p, BH-adjusted p within disease and cell type, Q, its p, I^2, the joint
#   model's log2FC, SE and adjusted p, and the log2FC of each cohort)

# --- cell 9 ---
import warnings

warnings.filterwarnings("ignore")

# --- cell 10 ---
import os  # CHANGED: run folder

import decoupler as dc
import h5py  # CHANGED: only the counts layer and the cell table are read
import anndata as ad  # CHANGED: as above
import numpy as np
import pandas as pd
import pertpy as pt
from scipy import stats  # CHANGED: meta-analysis

from config import DATASETS, OUT_ROOT, batch_of, excluded_labels, log  # CHANGED: pipeline paths; labels left out; batch
import rpy2.robjects as ro  # CHANGED: edgeR::filterByExpr with the design (cell 51)
from rpy2.robjects import numpy2ri, pandas2ri  # CHANGED: as above
from rpy2.robjects.conversion import localconverter  # CHANGED: as above

ro.r("suppressMessages(library(edgeR))")  # CHANGED: as above

# CHANGED: lamindb (ln.connect, ln.track) removed

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline input
out = run / "04b_percohort_de"  # CHANGED: pipeline output
out.mkdir(parents=True, exist_ok=True)  # CHANGED: as above
ncpu = len(os.sched_getaffinity(0))  # CHANGED: the job's CPUs (step 04, cell 58)

# --- cells 14-23 (data loading and cell/gene filtering of the notebook's dataset) replaced ---
# CHANGED: as step 04: the counts of the atlas of 02_integration.py with the manual cell-type labels of 03_annotation.py;
# cells without a label, labelled low quality or platelets are not used
with h5py.File(run / "atlas.h5ad", "r") as f:  # CHANGED: only the cell table, the gene table and the counts are read
    adata = ad.AnnData(obs=ad.io.read_elem(f["obs"]), var=ad.io.read_elem(f["var"]))  # CHANGED: as above
    adata.layers["counts"] = ad.io.read_elem(f["layers/counts"])  # CHANGED: as above
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)  # CHANGED: as above
adata.obs["cell_type"] = labels["manual_celltype_annotation"].reindex(adata.obs_names)  # CHANGED: as above
adata = adata[adata.obs["cell_type"].notna() & ~excluded_labels(adata.obs)].copy()  # CHANGED: as above
adata.obs["sex"] = adata.obs["inferred_sex"].astype(str)  # CHANGED: inferred sex as the design covariate (step 04)
adata.obs["batch"] = batch_of(adata.obs).values  # CHANGED: batch as in step 04; one level outside Perez
log(f"{adata.n_obs} labelled cells, {adata.obs['cell_type'].nunique()} cell types")  # CHANGED: log

results = []  # CHANGED: per-cohort results of all datasets, for the meta-analysis
for _, row in DATASETS.iterrows():  # CHANGED: one fit per dataset, its disease against its own controls
    dataset, case = row["dataset"], row["disease"]  # CHANGED: as above
    sub = adata[(adata.obs["dataset"] == dataset) & adata.obs["condition"].isin([case, "control"])].copy()  # CHANGED: as above
    log(f"{dataset} ({case}): {sub.n_obs} cells, {sub.obs['donor_id'].nunique()} donors")  # CHANGED: log

    # --- cell 26 ---
    # CHANGED: one pseudobulk sample per donor (step 04)
    sub.obs["sample"] = pd.Categorical(sub.obs["donor_id"].astype(str))

    # --- cell 28 ---
    adata_pb = dc.pp.pseudobulk(
        sub, sample_col="sample",  # CHANGED: the dataset's donors
        groups_col="cell_type",
        layer="counts",
        mode="sum"
    )
    del sub  # CHANGED: memory

    # --- cell 32 ---
    dc.pp.filter_samples(adata_pb, min_cells=10,min_counts=1000)

    res_ds = []  # CHANGED: one result table per cell type, written together
    for cell_type in sorted(adata_pb.obs["cell_type"].unique()):  # CHANGED: every cell type (step 04)
        # --- cell 47 ---
        adata_ct = adata_pb[adata_pb.obs["cell_type"] == cell_type].copy()  # CHANGED: as above
        n = adata_ct.obs["condition"].value_counts()  # CHANGED: a contrast needs donors of both conditions (step 04)
        if n.get(case, 0) < 2 or n.get("control", 0) < 2:  # CHANGED: as above
            log(f"  {cell_type}: skipped ({n.get(case, 0)} {case}, {n.get('control', 0)} control pseudobulk samples)")  # CHANGED: as above
            continue  # CHANGED: as above

        # --- cell 51 ---
        # CHANGED: edgeR::filterByExpr with the design of cell 53, as in step 04
        terms = [c for c in ["batch", "sex"] if adata_ct.obs[c].nunique() > 1] + ["condition"]  # CHANGED: step 04's design; batch remains only in Perez
        counts_ct = adata_ct.X.toarray() if hasattr(adata_ct.X, "toarray") else np.asarray(adata_ct.X)  # CHANGED: as above
        with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):  # CHANGED: as above
            ro.globalenv["y"] = counts_ct.T.astype(float)  # CHANGED: genes x samples, as edgeR expects
            ro.globalenv["md"] = adata_ct.obs[terms].astype(str).reset_index(drop=True)  # CHANGED: as above
        ro.r(f"design <- model.matrix(~ {' + '.join(terms)}, data = md); "  # CHANGED: as above
             "keep <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)")
        adata_ct = adata_ct[:, np.asarray(ro.r("keep"), dtype=bool)].copy()  # CHANGED: as above
        dc.pp.filter_by_prop(
            adata=adata_ct,
            min_prop=0.1,
            min_smpls=2,
        )

        # --- cell 53 ---
        design = "~ " + " + ".join(terms)  # CHANGED: step 04's design (the notebook: ~ label)
        pds2 = pt.tl.PyDESeq2(adata=adata_ct,design=design)  # CHANGED: as above

        # --- cell 54 ---
        pds2.fit()

        # --- cell 58 ---
        res_df = pds2.test_contrasts(pds2.contrast(
            column="condition",  # CHANGED: condition (label in the notebook)
            baseline="control",  # CHANGED: as above (ctrl)
            group_to_compare=case),  # CHANGED: as above (stim)
            n_cpus=ncpu,  # CHANGED: the job's CPUs (step 04)
        )

        # --- cell 59 ---
        res_ds.append(res_df.assign(dataset=dataset, disease=case, cell_type=cell_type, design=design,  # CHANGED: kept for the output and the meta-analysis
                                    n_case=int(n[case]), n_control=int(n["control"])))
        log(f"  {cell_type}: {design}; {int(n[case])} {case} vs {int(n['control'])} control; {adata_ct.n_vars} genes; "
            f"{int(((res_df['adj_p_value'] < 0.05) & (res_df['log_fc'].abs() >= 1)).sum())} DEGs")  # CHANGED: log
    if res_ds:  # CHANGED: all cell types of the dataset together
        res_ds = pd.concat(res_ds)
        res_ds.to_csv(out / f"de_{dataset}.tsv.gz", sep="\t", index=False)
        results.append(res_ds)

# CHANGED: fixed-effect (inverse-variance) meta-analysis of the per-cohort log2FC, per disease, cell type and gene
#: weight w = 1 / SE^2, estimate sum(w * log2FC) / sum(w), SE 1 / sqrt(sum(w)), Wald p, BH within disease
# and cell type; Cochran's Q = sum(w * (log2FC - estimate)^2) on k - 1 degrees of freedom; I^2 = max(0, (Q - (k - 1)) / Q).
# A gene enters from every cohort where it was tested (passed the filters and has a log2FC and SE).
pc = pd.concat(results).dropna(subset=["log_fc", "lfcSE"])
pc["w"] = 1 / pc["lfcSE"] ** 2
pc["wb"] = pc["w"] * pc["log_fc"]
key = ["disease", "cell_type", "variable"]
meta = pc.groupby(key).agg(k=("log_fc", "size"), sw=("w", "sum"), swb=("wb", "sum")).reset_index()
meta["meta_log_fc"] = meta["swb"] / meta["sw"]
meta["meta_se"] = 1 / np.sqrt(meta["sw"])
pc = pc.merge(meta[key + ["meta_log_fc"]], on=key)
pc["q"] = pc["w"] * (pc["log_fc"] - pc["meta_log_fc"]) ** 2
meta = meta.merge(pc.groupby(key)["q"].sum().rename("Q").reset_index(), on=key)
meta["Q_p_value"] = np.where(meta["k"] > 1, stats.chi2.sf(meta["Q"], np.maximum(meta["k"] - 1, 1)), np.nan)
meta["I2"] = np.where((meta["k"] > 1) & (meta["Q"] > 0), np.clip((meta["Q"] - (meta["k"] - 1)) / meta["Q"], 0, 1), np.nan)
meta["meta_p_value"] = 2 * stats.norm.sf(np.abs(meta["meta_log_fc"] / meta["meta_se"]))
meta["meta_adj_p_value"] = meta.groupby(["disease", "cell_type"])["meta_p_value"].transform(
    lambda p: stats.false_discovery_control(p.to_numpy(), method="bh"))
joint = pd.concat(pd.read_csv(run / "04_differential_expression" / f"{d}_vs_control" / "all" / "de_all_cell_types.tsv",
                              sep="\t").assign(disease=d) for d in ["SLE", "SjD"])
meta = meta.merge(joint[key + ["log_fc", "lfcSE", "adj_p_value"]].rename(
    columns={"log_fc": "joint_log_fc", "lfcSE": "joint_se", "adj_p_value": "joint_adj_p_value"}), on=key, how="left")
wide = pc.pivot_table(index=key, columns="dataset", values="log_fc").add_prefix("log_fc_").reset_index()
meta = meta.drop(columns=["sw", "swb"]).merge(wide, on=key, how="left")
meta.to_csv(out / "meta.tsv.gz", sep="\t", index=False)
for d, x in meta.dropna(subset=["joint_log_fc"]).groupby("disease"):
    jd = (x["joint_adj_p_value"] < 0.05) & (x["joint_log_fc"].abs() >= 1)
    log(f"{d}: {len(x)} genes in both; log2FC r (joint, meta) {np.corrcoef(x['joint_log_fc'], x['meta_log_fc'])[0, 1]:.3f}; "
        f"joint DEGs {int(jd.sum())}, same sign in the meta-analysis {(np.sign(x.loc[jd, 'joint_log_fc']) == np.sign(x.loc[jd, 'meta_log_fc'])).mean():.3f}, "
        f"tested in > 1 cohort {(x.loc[jd, 'k'] > 1).mean():.3f}, I2 > 0.5 {(x.loc[jd, 'I2'] > 0.5).mean():.3f}")  # CHANGED: log
log("done")  # CHANGED: log
