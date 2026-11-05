# Step 01 of the Python pipeline: quality control -> ambient RNA correction -> doublet detection
# (the gene filter runs once on the merged atlas), ported from sc-best-practices
# jupyter-book/preprocessing_visualization/quality_control.ipynb, commit 47edbd5 (2026-09-29; the
# same file at book HEAD dd6eb86, 2026-10-19). The notebook's code cells are copied in order; every
# changed line carries a "# CHANGED:" comment with the reason. Plot-only cells (18, 26, 29) and the
# lamindb cells (9, 66) are left out, as marked where they would be.
# Usage: python 01_quality_control.py <dataset> [--stage qc | library [i] | merge]  (no --stage: all three in one process;
#   stages: see the driver below)
# Output: $PIPELINE_OUT/py/datasets/<dataset>/01_quality_control/: cells.tsv (every labelled cell,
#   its metrics and the outcome of each filter), libraries.tsv, samples/<sample>.h5ad

# --- cell 7 ---
from typing import Literal

import anndata as ad
# CHANGED: lamindb not used; inputs are read from the pipeline's per-library files
import matplotlib.pyplot as plt  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns  # noqa: F401
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import csc_matrix
from scipy.special import logit
from scipy.stats import median_abs_deviation

import sys  # CHANGED: dataset given on the command line
from config import DATASETS, DBL_RATE_PER_1K, INPUT_ROOT, REF_DIR, dataset_dir, log  # CHANGED: pipeline paths

# Suppress verbose logging from Scanpy
sc.settings.verbosity = 0

# Set figure parameters for clean, minimal plots
sc.set_figure_params(dpi=80, facecolor="white", frameon=False)

# CHANGED: ln.connect / ln.track removed (no lamindb)

# --- cell 32 ---
import logging

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

rcb.logger.setLevel(logging.ERROR)

# CHANGED: %load_ext rpy2.ipython and %%R cells run as ro.r() in a script

# --- cell 33 ---
ro.r("suppressMessages({library(decontX); library(SingleCellExperiment)})")  # CHANGED: DecontX instead of SoupX

# --- cell 56 ---
ro.r("suppressMessages({library(scater); library(scDblFinder); library(SingleCellExperiment); library(BiocParallel)})")  # CHANGED: %%R cell as ro.r()


# --- cell 20 ---
def mad_thresholds(values: pd.Series, nmads: float) -> tuple[float, float]:
    """Thresholds `nmads` median absolute deviations below and above the median of `values`."""
    median, mad = np.median(values), median_abs_deviation(values)
    return median - nmads * mad, median + nmads * mad


def is_outlier(
    adata: ad.AnnData,
    metric: str,
    nmads: float,
    side: Literal["both", "upper"] = "both",
) -> pd.Series:
    """Flag the cells beyond the MAD thresholds of `metric`, on both sides or only above."""
    values = adata.obs[metric]
    lower, upper = mad_thresholds(values, nmads)
    return values > upper if side == "upper" else (values < lower) | (values > upper)


def qc_metrics(adata):  # CHANGED: notebook cells 14-28 split in two, so that the per-cell metrics are computed per library; cells 14-22 (per-cell metrics; each cell's values depend only on its own counts)
    # --- cell 14 ---
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    # --- cell 16 ---
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    # CHANGED: bare `adata` (display only) left out
    # --- cell 18 (plots) left out ---
    # --- cell 22 ---
    mt_fraction = (adata.obs["total_counts_mt"] + 1) / (adata.obs["total_counts"] + 2)
    adata.obs["logit_mt_fraction"] = logit(mt_fraction)
    return adata  # CHANGED: end of qc_metrics


def qc_flags(adata):  # CHANGED: notebook cells 24-28 (MAD thresholds per sample; only adata.obs is read)
    # --- cell 24 ---
    qc_filters = {
        "log1p_total_counts": (5, "both"),
        "log1p_n_genes_by_counts": (5, "both"),
        "pct_counts_in_top_20_genes": (5, "both"),
        "logit_mt_fraction": (3, "upper"),
    }

    flags = {
        metric: is_outlier(adata, metric, nmads, side)
        for metric, (nmads, side) in qc_filters.items()
    }
    adata.obs["outlier"] = (
        flags["log1p_total_counts"]
        | flags["log1p_n_genes_by_counts"]
        | flags["pct_counts_in_top_20_genes"]
    )
    adata.obs["mt_outlier"] = flags["logit_mt_fraction"] & (adata.obs["pct_counts_mt"] > 8)
    # CHANGED: adata.obs[["outlier", "mt_outlier"]].value_counts() (display only) left out
    for metric, flag in flags.items():  # CHANGED: keep each flag for the per-filter report
        adata.obs[f"flag_{metric}"] = flag  # CHANGED: as above
    # --- cell 26 (plot) left out ---
    # --- cell 28 ---
    log(f"Total number of cells: {adata.n_obs}")  # CHANGED: print -> log
    record = adata.obs.copy()  # CHANGED: keep the metrics of all cells for the report
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    log(f"Number of cells after filtering of low quality cells: {adata.n_obs}")  # CHANGED: print -> log
    # --- cell 29 (plot) left out ---
    return adata, record  # CHANGED: function return


def ambient(adata, data_tod):  # CHANGED: notebook cells 35-51 as a function, run once per 10x library (ambient RNA is shared by the library's samples)
    # --- cell 35 ---
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp, target_sum=1e4)
    sc.pp.log1p(adata_pp)
    # --- cell 37 ---
    sc.pp.pca(adata_pp, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
    sc.pp.neighbors(adata_pp, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
    sc.tl.leiden(
        adata_pp, key_added="soupx_groups", flavor="igraph", n_iterations=2, directed=False,
        random_state=0,  # CHANGED: scanpy default random_state=0 written out (seed)
    )

    # Preprocess variables for SoupX
    adata.obs["soupx_groups"] = adata_pp.obs["soupx_groups"]
    # --- cell 39 ---
    del adata_pp
    # --- cell 41 ---
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T
    # --- cell 43 ---
    # CHANGED: the raw droplets come from the library's raw_feature_bc_matrix.h5 (data_tod argument)
    # --- cell 46 ---
    data_csc = data.tocsc()
    # CHANGED: data_tod_csc and x_tod/i_tod/p_tod/dims_tod built below, only when raw droplets exist

    # Extract sparse components and cast to correct types
    x = data_csc.data.astype(np.float64)
    i = data_csc.indices.astype(np.int32)
    p = data_csc.indptr.astype(np.int32)
    dims = np.array(data_csc.shape, dtype=np.int32)

    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        ro.globalenv["x"] = x
        ro.globalenv["i"] = i
        ro.globalenv["p"] = p
        ro.globalenv["dims"] = dims

        ro.globalenv["genes"] = np.array(genes)
        ro.globalenv["cells"] = np.array(cells)
        ro.globalenv["soupx_groups"] = adata.obs["soupx_groups"].to_numpy()
        # CHANGED: genes_raw / cells_raw not passed: the raw droplets are already subset to the cells' genes
        ro.globalenv["has_tod"] = data_tod is not None  # CHANGED: deposited matrices have no raw droplets
        if data_tod is not None:  # CHANGED: as above
            data_tod_csc = data_tod.tocsc()
            ro.globalenv["x_tod"] = data_tod_csc.data.astype(np.float64)
            ro.globalenv["i_tod"] = data_tod_csc.indices.astype(np.int32)
            ro.globalenv["p_tod"] = data_tod_csc.indptr.astype(np.int32)
            ro.globalenv["dims_tod"] = np.array(data_tod_csc.shape, dtype=np.int32)
    # --- cell 48 ---
    ro.r("""
library(Matrix)

# Manually coerce types to avoid "array" class errors
x <- as.numeric(x)
i <- as.integer(i)
p <- as.integer(p)
dims <- as.integer(dims)

# CHANGED: x_tod/i_tod/p_tod/dims_tod coerced below, only when raw droplets exist

# Reconstruct sparse matrices
data <- new("dgCMatrix",
            Dim = dims,
            x = x,
            i = i,
            p = p)

# Assign row and column names
rownames(data) <- genes
colnames(data) <- cells
# CHANGED: data_tod built in the DecontX block below; its names are not needed by DecontX

# CHANGED: DecontX instead of SoupX, with the notebook's clusters (z) and seed 12345
bg <- NULL
if (has_tod) {
  bg <- SingleCellExperiment(list(counts = new("dgCMatrix", Dim = as.integer(dims_tod), x = as.numeric(x_tod),
                                               i = as.integer(i_tod), p = as.integer(p_tod))))
}
sce <- decontX(SingleCellExperiment(list(counts = data)), z = soupx_groups, background = bg, seed = 12345)
out <- as(round(decontXcounts(sce)), "CsparseMatrix")  # CHANGED: rounding as SoupX adjustCounts(roundToInt = TRUE)
contamination <- sce$decontX_contamination
""")
    # --- cell 50 ---
    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        out_py = ro.conversion.rpy2py(ro.globalenv["out"])
        contamination = np.array(ro.globalenv["contamination"])  # CHANGED: keep DecontX's per-cell estimate

    x = np.array(out_py.slots["x"])
    i = np.array(out_py.slots["i"])
    p = np.array(out_py.slots["p"])
    shape = tuple(out_py.slots["Dim"])

    out_matrix = csc_matrix((x, i, p), shape=shape)
    # --- cell 51 ---
    adata.layers["counts"] = adata.X.copy()
    adata.layers["decontX_counts"] = out_matrix.T.tocsr()  # CHANGED: DecontX instead of SoupX
    adata.X = adata.layers["decontX_counts"]  # CHANGED: as above
    adata.obs["decontX_contamination"] = contamination  # CHANGED: as above
    return adata  # CHANGED: function return


# --- cell 53 ---
# CHANGED: the gene filter (sc.pp.filter_genes(adata, min_cells=20)) is not applied per sample here but
# once on the merged atlas: joining samples with different gene sets would replace real
# counts with zeros, and scDblFinder selects its own features from the library's cells.


def doublets(adata, dbr):  # CHANGED: notebook cells 57-64 as a function, run once per 10x library
    # --- cell 57 ---
    data_mat = adata.X.T
    # --- cell 60 ---
    data_mat = adata.X.T.tocsc()

    x = data_mat.data.astype(np.float64)
    i = data_mat.indices.astype(np.int32)
    p = data_mat.indptr.astype(np.int32)
    dims = np.array(data_mat.shape, dtype=np.int32)

    with localconverter(ro.default_converter + numpy2ri.converter):
        ro.globalenv["x"] = x
        ro.globalenv["i"] = i
        ro.globalenv["p"] = p
        ro.globalenv["dims"] = dims
        ro.globalenv["dbr"] = dbr  # CHANGED: explicit expected doublet rate
    # --- cell 62 ---
    ro.r("""
x <- as.numeric(x)
i <- as.integer(i)
p <- as.integer(p)
dims <- as.integer(dims)

data_mat <- new("dgCMatrix", Dim = dims, x = x, i = i, p = p)

set.seed(123)
sce <- scDblFinder(SingleCellExperiment(list(counts = data_mat)), dbr = dbr)  # CHANGED: dbr passed

doublet_score <- sce$scDblFinder.score
doublet_class <- sce$scDblFinder.class
""")
    doublet_score = np.array(ro.globalenv["doublet_score"])  # CHANGED: %%R -o outputs read from R
    doublet_class = np.array(ro.r("as.character(doublet_class)"))  # CHANGED: as above
    # --- cell 64 ---
    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    log(adata.obs.scDblFinder_class.value_counts().to_dict())  # CHANGED: display -> log
    return adata  # CHANGED: function return


# CHANGED: everything below loops the notebook's steps over the dataset's samples and libraries, in three stages:
#   python 01_quality_control.py <dataset> --stage qc           A: per-cell metrics, MAD flags per sample
#   python 01_quality_control.py <dataset> --stage library [i]  B: library i (default SLURM_ARRAY_TASK_ID):
#                                                                   DecontX, expected doublet rate and check, scDblFinder
#   python 01_quality_control.py <dataset> --stage merge        C: cells.tsv, libraries.tsv, samples/*.h5ad
#   python 01_quality_control.py <dataset>                      A, B for every library, C, in one process
# QC needs a sample's cells from every library it was sequenced in (a Perez sample spans 2-16 libraries),
# so the per-sample thresholds are stage A. Stage files go to out/stages/ and are removed after a successful merge.
import json  # CHANGED: stage bookkeeping
import os  # CHANGED: as above
import shutil  # CHANGED: as above

ds = sys.argv[1]
stage = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "all"  # CHANGED: staged run
mode = DATASETS.set_index("dataset").doublets[ds]
if mode not in ("single", "pooled_cleared", "removed_by_authors"):  # CHANGED: removed_by_authors, deposited cells after the authors' doublet removal
    sys.exit(f"{ds}: doublet design '{mode}' is not a known doublet design (pipeline_py/datasets.tsv)")
libs = pd.read_csv(INPUT_ROOT / ds / "libraries.tsv", sep="\t", dtype=str, keep_default_na=False)
_ref = pd.read_csv(REF_DIR / "genes" / "genes.gtf.gz", sep="\t", comment="#", header=None, usecols=[2, 8], dtype=str)  # CHANGED:
_ref = _ref[_ref[2] == "gene"][8]  # GRCh38-2024-A genes in GTF order, which is Cell Ranger's feature order
REF_NAMES = pd.Series(ad.utils.make_index_unique(pd.Index(_ref.str.extract(r'gene_name "([^"]+)"')[0])).values,
                      index=_ref.str.extract(r'gene_id "([^"]+)"')[0].values)  # CHANGED: unique name of every reference gene
out = dataset_dir(ds, "01_quality_control")
(out / "samples").mkdir(parents=True, exist_ok=True)
st = out / "stages"  # CHANGED: intermediate files of the staged run
QC_VAR = ["n_cells_by_counts", "mean_counts", "log1p_mean_counts", "pct_dropout_by_counts", "total_counts",
          "log1p_total_counts"]  # CHANGED: gene-level columns of calculate_qc_metrics, which differ between samples


def read_library(L):  # CHANGED: cell 9 replaced: the library's Cell Ranger filtered matrix, cells with a donor label
    a = sc.read_10x_h5(L.filtered_h5)
    # --- cell 11 ---
    a.var_names = REF_NAMES[a.var.gene_ids].values  # CHANGED: a.var_names_make_unique() over the whole reference, so that a
    # dataset holding only some reference genes (deposited route) names them as the other datasets do; equal to
    # var_names_make_unique() for a Cell Ranger library (full reference list)
    labelled = pd.read_csv(L.cells_tsv, sep="\t", dtype=str).set_index("barcode")
    a = a[a.obs_names.isin(labelled.index)].copy()
    a.obs = a.obs.join(labelled[["sample_id", "donor_id"]])
    a.obs["library"] = L.library
    a.obs_names = L.library + "_" + a.obs_names
    return a


def ambient_library(L, a):  # CHANGED: 2. ambient RNA correction of one library, on its QC-passed cells, raw droplets as background
    data_tod = None
    if L.raw_h5:
        raw = sc.read_10x_h5(L.raw_h5)
        raw.var_names_make_unique()
        called = sc.read_10x_h5(L.filtered_h5).obs_names
        raw = raw[(~raw.obs_names.isin(called)) & (np.asarray(raw.X.sum(axis=1)).ravel() > 0)]  # empty droplets
        data_tod = raw[:, a.var_names].X.T
    a = ambient(a, data_tod)
    log(f"{L.library}: DecontX on {a.n_obs} cells, background {'none' if data_tod is None else data_tod.shape[1]}")
    return a


def doublets_library(L, a, lab, record):  # CHANGED: 3. doublets of one library, on its QC-passed, corrected cells with all genes
    empty = np.asarray(a.X.sum(axis=1)).ravel() == 0  # cells emptied by ambient correction
    record.loc[a.obs_names[empty], "empty_after_ambient"] = True
    a = a[~empty].copy()
    n_called = float(L.n_called)
    if mode == "removed_by_authors":  # CHANGED: deposited cells after the authors' doublet removal: no second pass
        record.loc[a.obs_names, "scDblFinder_class"] = "not_run_removed_by_authors"
        row = {"library": L.library, "mode": mode, "n_called": n_called, "n_called_source": L.n_called_source,
               "n_donors": lab.donor_id.nunique(), "sum_sq_donor_prop": np.nan, "dbr_used": np.nan,
               "n_empty_after_ambient": int(empty.sum()), "n_doublets": 0, "unlabelled_fraction": np.nan,
               "predicted_cross_donor": np.nan, "check": "n/a"}
        return a, row
    prop = lab.donor_id.value_counts(normalize=True)
    same_donor = float((prop ** 2).sum())
    total_rate = DBL_RATE_PER_1K * n_called / 1000
    dbr = total_rate * same_donor if mode == "pooled_cleared" else total_rate
    a = doublets(a, dbr)
    record.loc[a.obs_names, "scDblFinder_score"] = a.obs.scDblFinder_score
    record.loc[a.obs_names, "scDblFinder_class"] = a.obs.scDblFinder_class
    unlabelled = 1 - len(lab) / n_called
    cross = total_rate * (1 - same_donor) if mode == "pooled_cleared" else np.nan
    row = {"library": L.library, "mode": mode, "n_called": n_called, "n_called_source": L.n_called_source,
           "n_donors": len(prop), "sum_sq_donor_prop": round(same_donor, 4), "dbr_used": round(dbr, 4),
           "n_empty_after_ambient": int(empty.sum()), "n_doublets": int((a.obs.scDblFinder_class == "doublet").sum()),
           "unlabelled_fraction": round(unlabelled, 4), "predicted_cross_donor": round(cross, 4),
           "check": ("ok" if unlabelled >= cross else "FLAG") if mode == "pooled_cleared" else "n/a"}
    return a[a.obs.scDblFinder_class == "singlet"].copy(), row


def write_outputs(record, lib_rows, final):  # CHANGED: outputs: every labelled cell with each filter's outcome; one h5ad per sample (singlets)
    pd.DataFrame(lib_rows).to_csv(out / "libraries.tsv", sep="\t", index=False)
    record["kept"] = record.scDblFinder_class.isin(["singlet", "not_run_removed_by_authors"])  # CHANGED: singlets, or deposited cells after the authors' doublet removal
    cols = ["library", "sample_id", "donor_id", "total_counts", "n_genes_by_counts", "pct_counts_in_top_20_genes",
            "pct_counts_mt", "logit_mt_fraction", "pct_counts_ribo", "pct_counts_hb", "flag_log1p_total_counts",
            "flag_log1p_n_genes_by_counts", "flag_pct_counts_in_top_20_genes", "flag_logit_mt_fraction", "outlier",
            "mt_outlier", "decontX_contamination", "empty_after_ambient", "scDblFinder_score", "scDblFinder_class", "kept"]
    record[cols].rename_axis("cell_id").to_csv(out / "cells.tsv", sep="\t")
    res = ad.concat(final, merge="same")
    for smp, idx in res.obs.groupby("sample_id", observed=True).indices.items():
        res[idx].copy().write_h5ad(out / "samples" / f"{smp}.h5ad")
    log(f"{ds}: {len(record)} labelled cells; QC removed {int((record.outlier | record.mt_outlier).sum())}, "
        f"empty after ambient {int(record.empty_after_ambient.sum())}, doublets {int((record.scDblFinder_class == 'doublet').sum())}; "
        f"{int(record.kept.sum())} kept")


def new_record_columns(record):  # CHANGED: columns filled by the doublet step
    record["empty_after_ambient"] = False
    record["scDblFinder_score"] = np.nan
    record["scDblFinder_class"] = ""
    return record


def stage_qc():  # CHANGED: A. per-cell metrics library by library, then the MAD flags per sample over all its libraries
    shutil.rmtree(st, ignore_errors=True)  # a new run never mixes with an earlier run's stage files
    st.mkdir()
    obs = []
    for _, L in libs.iterrows():
        obs.append(qc_metrics(read_library(L)).obs)
        log(f"{L.library}: QC metrics of {len(obs[-1])} labelled cells")
    obs = pd.concat(obs)
    log(f"{ds}: {len(obs)} labelled cells in {len(libs)} libraries")
    records = [qc_flags(ad.AnnData(obs=obs.iloc[idx].copy()))[1]
               for smp, idx in obs.groupby("sample_id", observed=True).indices.items()]
    record = pd.concat(records)  # record order: samples in groupby order, then library and barcode order
    record.to_pickle(st / "qc_record.pkl")  # pickle keeps dtypes and exact values
    (st / "qc.done").write_text(json.dumps({"n_cells": len(record), "n_libraries": len(libs)}))


def stage_library(i):  # CHANGED: B. one library: ambient RNA correction and doublets
    L = libs.iloc[i]
    if not (st / "qc.done").exists():
        sys.exit(f"{ds}: stage A (qc) has not finished")
    rec_all = pd.read_pickle(st / "qc_record.pkl")
    record = rec_all[rec_all.library == L.library].copy()  # this library's cells in record order
    a = read_library(L)[record.index].copy()
    a = qc_metrics(a)  # var flags (mt, ribo, hb); per-cell values equal stage A's
    a.var = a.var.drop(columns=QC_VAR)  # gene-level QC columns differ between samples; not kept
    a.obs = record.loc[a.obs_names].copy()
    a = a[(~a.obs.outlier) & (~a.obs.mt_outlier)].copy()
    a = ambient_library(L, a)
    record = record.join(a.obs[["decontX_contamination"]])
    labelled = record.copy()
    record = new_record_columns(record)
    a, row = doublets_library(L, a, labelled, record)
    d = st / "libraries"
    d.mkdir(parents=True, exist_ok=True)
    record.to_pickle(d / f"{L.library}.record.pkl")
    pd.to_pickle(row, d / f"{L.library}.row.pkl")  # the libraries.tsv row as a dict (types kept)
    a.write_h5ad(d / f"{L.library}.singlets.h5ad", convert_strings_to_categoricals=False)  # obs as in memory
    (d / f"{L.library}.done").write_text(json.dumps({"n_record": len(record), "n_singlets": a.n_obs}))  # written last
    log(f"{L.library}: stage B done, {len(record)} labelled cells, {a.n_obs} singlets")


def stage_merge():  # CHANGED: C. refuses to run unless every library's stage B output is present and complete
    d = st / "libraries"
    missing = []
    for _, L in libs.iterrows():
        done = d / f"{L.library}.done"
        files = [d / f"{L.library}.{x}" for x in ("record.pkl", "row.pkl", "singlets.h5ad")]
        absent = [f.name for f in [done] + files if not f.exists()]
        if absent:
            missing.append(f"{L.library}: missing {', '.join(absent)}")
            continue
        n = json.loads(done.read_text())
        if len(pd.read_pickle(files[0])) != n["n_record"] or ad.read_h5ad(files[2], backed="r").n_obs != n["n_singlets"]:
            missing.append(f"{L.library}: outputs do not match its .done counts")
    if missing:
        sys.exit(f"{ds}: merge refused, {len(missing)} problem(s): " + "; ".join(missing[:20]))
    rec_all = pd.read_pickle(st / "qc_record.pkl")
    record = pd.concat([pd.read_pickle(d / f"{L.library}.record.pkl") for _, L in libs.iterrows()])
    if len(record) != len(rec_all) or not record.index.isin(rec_all.index).all():
        sys.exit(f"{ds}: merge refused, library records do not cover stage A's cells")
    record = record.loc[rec_all.index]  # stage A's row order
    lib_rows = [pd.read_pickle(d / f"{L.library}.row.pkl") for _, L in libs.iterrows()]
    final = [ad.read_h5ad(d / f"{L.library}.singlets.h5ad") for _, L in libs.iterrows()]
    write_outputs(record, lib_rows, final)
    shutil.rmtree(st)  # stage files duplicate the outputs; removed once the merge has written them
    log(f"{ds}: stage files removed")


if stage == "all":
    stage_qc()
    for i in range(len(libs)):
        stage_library(i)
    stage_merge()
elif stage == "qc":
    stage_qc()
elif stage == "library":
    k = sys.argv.index("--stage") + 2
    stage_library(int(sys.argv[k]) if len(sys.argv) > k else int(os.environ["SLURM_ARRAY_TASK_ID"]))
elif stage == "merge":
    stage_merge()
else:
    sys.exit(f"unknown stage {stage}")
