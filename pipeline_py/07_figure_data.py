# Data of every figure and table, main and supplementary, drawn by
# 08_figures.R, which reads nothing else. Not a pipeline step and not a notebook port.
# Usage: run by make_figures.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration; FIG_DATA=<folder>
#        writes elsewhere than the default output folder)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/: atlas.h5ad (obs, X_umap, X_umap_unintegrated, X), markers_unshared.h5ad,
#        03_annotation/cells_manual_annotation.tsv.gz, 04_differential_expression/, 04b_percohort_de/, 05_composition/,
#        05b_reference_sensitivity/, 06_pathways/; $PIPELINE_OUT/py/datasets/<dataset>/01_quality_control/cells.tsv
# Output: $FIG_DATA, default $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/07_figures/data/
#   Computed here:
#   cells.tsv.gz         UMAP coordinates, dataset, condition and cell type of the cells analysed in steps 04-06 (labelled,
#                        not "low quality" or "platelets"; Fig 2)
#   table1.tsv           per dataset, over the atlas cells (low-quality and platelet clusters removed): cases, controls,
#                        cells analysed, median genes per cell; and cells after QC before that removal (Table 1)
#   indirect.tsv.gz      SLE vs SjD indirect comparison per gene and cell type, genes tested in both contrasts (Fig 3b,
#                        Supplementary Table SLE vs SjD): difference of the two log2FC (each disease vs its own controls),
#                        SE = sqrt(SE_SLE^2 + SE_SjD^2), Wald z, BH within cell type; equivalence by two one-sided tests
#                        (margin 1 log2FC; p_tost, padj_tost BH within cell type); deg_sle, deg_sjd (adjusted p < 0.05 and
#                        |log2FC| >= 1 vs HC); class of the genes changed vs HC in either disease
#   indirect_gsea.tsv    fgsea on the z statistics of the indirect comparison (Fig 4b), as 06_pathways.py: Hallmark 2025.1,
#                        fgseaMultilevel defaults, set.seed(42) before each cell type, BH across all cell type x gene set tests
#   composition_counts.tsv  cells per donor and cell type, with dataset and condition: the cells of step 05's scCODA model
#                        (labelled, not "low quality" or "platelets"); a cell type a donor lacks is a row with n = 0
#   qc_counts.tsv        cells per dataset after each step 01 filter, in the atlas and after the annotation
#                        exclusion: recorded counts only, the per-sample thresholds are not recorded
#   atlas_cells.tsv.gz   per atlas cell: dataset, condition, UMAP before integration (obsm X_umap_unintegrated) and the
#                        integrated UMAP (X_umap), annotation cluster and its label (supplementary atlas figure, a-c)
#   atlas_dotplot.tsv    per final cell type and marker of config.SUPP_DOTPLOT_MARKERS: mean expression (X, log1p of
#                        library-size normalised counts) and fraction of cells with counts, as sc.pl.dotplot; markers
#                        outside the shared genes from markers_unshared.h5ad, over the cells of the datasets that measured
#                        them, as 03_annotation.py (supplementary atlas figure, d)
#   Copied unchanged from the step outputs:
#   de_<contrast>.tsv    04_differential_expression/<contrast>/all/de_all_cell_types.tsv (SLE_vs_control, SjD_vs_control)
#   gsea_<contrast>.tsv  06_pathways/<contrast>/gsea.tsv
#   composition_<contrast>_effects.tsv, composition_<contrast>_credible.tsv
#                        05_composition/<contrast>/effects.tsv and credible_effects.tsv (SLE_vs_control, SjD_vs_control,
#                        SjD_vs_SLE); composition_reference.txt: 05_composition/reference_cell_type.txt
#   reference_sensitivity.tsv  05b_reference_sensitivity/reference_sensitivity.tsv
#   percohort_meta.tsv.gz      04b_percohort_de/meta.tsv.gz
import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scipy.sparse as sp
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.stats import false_discovery_control, norm

from config import DATASETS, OUT_ROOT, PIPELINE_DIR, SUPP_DOTPLOT_MARKERS, dataset_dir, excluded_labels, log

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")
out = Path(os.environ.get("FIG_DATA", run / "07_figures" / "data"))
out.mkdir(parents=True, exist_ok=True)
rcb.consolewrite_warnerror = lambda s: None


def to_r(name, df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv[name] = ro.conversion.py2rpy(df)


def from_r(name):
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(ro.globalenv[name])


# --- step outputs read by 08_figures.R, copied unchanged ---
CONTRASTS = ["SLE_vs_control", "SjD_vs_control"]
copies = {**{f"de_{c}.tsv": run / "04_differential_expression" / c / "all" / "de_all_cell_types.tsv" for c in CONTRASTS},
          **{f"gsea_{c}.tsv": run / "06_pathways" / c / "gsea.tsv" for c in CONTRASTS},
          **{f"composition_{c}_{k}.tsv": run / "05_composition" / c / f"{f}.tsv"
             for c in CONTRASTS + ["SjD_vs_SLE"] for k, f in (("effects", "effects"), ("credible", "credible_effects"))},
          "composition_reference.txt": run / "05_composition" / "reference_cell_type.txt",
          "reference_sensitivity.tsv": run / "05b_reference_sensitivity" / "reference_sensitivity.tsv",
          "percohort_meta.tsv.gz": run / "04b_percohort_de" / "meta.tsv.gz"}
for name, src in copies.items():
    shutil.copyfile(src, out / name)
log(f"copied {len(copies)} step outputs")

# --- atlas cells: obs, both UMAPs and the manual labels ---
with h5py.File(run / "atlas.h5ad", "r") as f:
    obs = ad.io.read_elem(f["obs"])[["dataset", "condition", "donor_id", "sample_id", "n_genes_by_counts"]].copy()
    um = {k: f["obsm"][k][:] for k in ("X_umap_unintegrated", "X_umap")}
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)
obs["cell_type"] = labels["manual_celltype_annotation"].reindex(obs.index)
excluded = excluded_labels(obs)
analysed = (obs["cell_type"].notna() & ~excluded).to_numpy()  # the cells of steps 04-06

# --- Table 1 ---
rows = []
# Donors, cells and median genes count the cells of the atlas: after quality control and after the low-quality and
# platelet clusters are removed at annotation (cells_after_qc keeps the count before that removal, for reference)
for _, d in DATASETS.iterrows():
    inds = (obs["dataset"] == d["dataset"]).to_numpy()
    o = obs[inds & analysed]
    rows.append({"dataset": d["dataset"], "label": d["label"], "disease": d["disease"],
                 "cases": o.loc[o.condition == d["disease"], "donor_id"].nunique(),
                 "controls": o.loc[o.condition == "control", "donor_id"].nunique(),
                 "cells_after_qc": int(inds.sum()), "cells_analysed": len(o),
                 "median_genes_per_cell": int(o["n_genes_by_counts"].median())})
pd.DataFrame(rows).to_csv(out / "table1.tsv", sep="\t", index=False)
log("Table 1:\n" + pd.DataFrame(rows).to_string(index=False))

# --- Fig 2: UMAP and labels of the cells analysed ---
cells = obs.loc[analysed, ["dataset", "condition", "cell_type"]].assign(UMAP_1=um["X_umap"][analysed, 0],
                                                                         UMAP_2=um["X_umap"][analysed, 1])
cells = cells[["UMAP_1", "UMAP_2", "dataset", "condition", "cell_type"]]
cells.round({"UMAP_1": 3, "UMAP_2": 3}).to_csv(out / "cells.tsv.gz", sep="\t", index_label="cell_id")
log(f"cells.tsv.gz: {len(cells)} cells analysed of {len(obs)}")

# --- supplementary atlas figure, a-c: the two UMAPs and the annotation clusters ---
ac = pd.DataFrame({"dataset": obs.dataset.astype(str).values, "condition": obs.condition.astype(str).values,
                   "umap_unint_1": um["X_umap_unintegrated"][:, 0], "umap_unint_2": um["X_umap_unintegrated"][:, 1],
                   "umap_1": um["X_umap"][:, 0], "umap_2": um["X_umap"][:, 1],
                   "cluster": labels["cluster"].reindex(obs.index).values, "label": obs.cell_type.values}, index=obs.index)
assert ac.cluster.notna().all(), "cells without an annotation cluster"
ac.rename_axis("cell_id").to_csv(out / "atlas_cells.tsv.gz", sep="\t")
log(f"atlas_cells.tsv.gz: {len(ac)} cells, {ac.cluster.nunique()} clusters, "
    f"{ac.cluster[excluded.values].nunique()} of them excluded at annotation")

# --- QC counts per dataset. Doublets removed: scDblFinder class "doublet"; kept at step 01: column kept (singlets,
# or every cell of a deposit whose doublets the authors removed); between them, the cells empty after ambient
# correction. Between kept and the atlas: the Perez later visits (02_integration.py) and the cells with no
# counts in the atlas genes. ---
rows = []
for ds in DATASETS.dataset:
    c = pd.read_csv(dataset_dir(ds, "01_quality_control") / "cells.tsv", sep="\t",
                    usecols=["sample_id", "outlier", "mt_outlier", "empty_after_ambient", "scDblFinder_class", "kept"])
    in_ds = obs.dataset.astype(str) == ds
    rows.append({"dataset": ds, "samples": obs.sample_id[in_ds].nunique(), "samples_qc": c.sample_id.nunique(), "before_qc": len(c),
                 "after_outlier": int((~c.outlier).sum()), "after_mt": int((~(c.outlier | c.mt_outlier)).sum()),
                 "doublets": int((c.scDblFinder_class == "doublet").sum()),
                 "doublets_run": bool(c.scDblFinder_class.isin(["singlet", "doublet"]).any()),
                 "empty_after_ambient": int(c.empty_after_ambient.sum()), "kept": int(c.kept.sum()),
                 "in_atlas": int(in_ds.sum()), "excluded_annotation": int((in_ds & excluded).sum()),
                 "final": int((in_ds & analysed).sum())})
qc = pd.DataFrame(rows)
qc.to_csv(out / "qc_counts.tsv", sep="\t", index=False)
log(f"qc_counts.tsv: {len(qc)} datasets, {qc.samples.sum()} samples, {qc.final.sum()} final cells")

# --- supplementary atlas figure, d: dot plot of the markers over the final cell types; X read in row blocks (the whole
# matrix does not fit in memory), only the marker columns kept ---
genes = list(dict.fromkeys(g for v in SUPP_DOTPLOT_MARKERS.values() for g in v))
mx = ad.read_h5ad(run / "markers_unshared.h5ad")
assert mx.obs_names.equals(obs.index), "markers_unshared.h5ad: cells differ from the atlas"
with h5py.File(run / "atlas.h5ad", "r") as f:
    var = ad.io.read_elem(f["var"]).index
    shared = [g for g in genes if g in var]
    missing = [g for g in genes if g not in var and g not in mx.var_names]
    assert not missing, f"dot-plot markers in neither the atlas nor markers_unshared.h5ad: {missing}"
    col = var.get_indexer(shared)
    X, ip = f["X"], f["X"]["indptr"][:]
    blocks, step = [], 100_000
    for i in range(0, len(obs), step):
        j = min(i + step, len(obs))
        b = sp.csr_matrix((X["data"][ip[i]:ip[j]], X["indices"][ip[i]:ip[j]], ip[i:j + 1] - ip[i]), shape=(j - i, len(var)))
        blocks.append(b[:, col])
    M = sp.vstack(blocks).tocsc()
ct = obs.cell_type.astype(str).values
rows = []
for g in genes:
    if g in shared:
        v, keep, n_ds = M[:, shared.index(g)].toarray().ravel(), analysed, len(mx.uns["datasets"])
    else:
        meas = [d for d in mx.uns["datasets"] if mx.var.loc[g, f"measured_{d}"]]
        v = mx[:, g].X.toarray().ravel()
        keep, n_ds = analysed & obs.dataset.astype(str).isin(meas).values, len(meas)
    d = pd.DataFrame({"cell_type": ct[keep], "v": v[keep]})
    r = d.groupby("cell_type").v.agg(mean="mean", fraction=lambda x: (x > 0).mean()).reset_index()
    rows.append(r.assign(gene=g, n_datasets=n_ds, n_datasets_all=len(mx.uns["datasets"])))
dot = pd.concat(rows)[["gene", "cell_type", "mean", "fraction", "n_datasets", "n_datasets_all"]]
dot.to_csv(out / "atlas_dotplot.tsv", sep="\t", index=False)
log(f"atlas_dotplot.tsv: {len(genes)} markers ({len(genes) - len(shared)} outside the shared genes: "
    + ", ".join(f"{g} {int(dot[dot.gene == g].n_datasets.iloc[0])}/{len(mx.uns['datasets'])}" for g in genes if g not in shared)
    + f"), {dot.cell_type.nunique()} cell types")

# --- composition counts: cells per donor and cell type, the cells of step 05 ---
o = obs.loc[analysed, ["donor_id", "dataset", "condition", "cell_type"]].astype(str)
counts = o.groupby(["donor_id", "cell_type"]).size().unstack(fill_value=0).stack().rename("n").reset_index()
donors = o.drop_duplicates("donor_id").set_index("donor_id")[["dataset", "condition"]]
counts = counts.join(donors, on="donor_id")[["donor_id", "dataset", "condition", "cell_type", "n"]]
counts.to_csv(out / "composition_counts.tsv", sep="\t", index=False)
log(f"composition_counts.tsv: {counts.donor_id.nunique()} donors, {counts.cell_type.nunique()} cell types, "
    f"{int(counts.n.sum())} cells")

# --- SLE vs SjD indirect comparison: Wald test of the difference, equivalence test and classes (Fig 3b; Methods) ---
cols = ["variable", "cell_type", "log_fc", "lfcSE", "stat", "adj_p_value"]
sle = pd.read_csv(out / "de_SLE_vs_control.tsv", sep="\t", usecols=cols)
sjd = pd.read_csv(out / "de_SjD_vs_control.tsv", sep="\t", usecols=cols)
m = sle.merge(sjd, on=["cell_type", "variable"], suffixes=("_sle", "_sjd"))
m = m.dropna(subset=["log_fc_sle", "lfcSE_sle", "log_fc_sjd", "lfcSE_sjd"])
m["diff"] = m["log_fc_sle"] - m["log_fc_sjd"]
m["se"] = np.sqrt(m["lfcSE_sle"] ** 2 + m["lfcSE_sjd"] ** 2)
m["z"] = m["diff"] / m["se"]
m["pvalue"] = 2 * norm.sf(np.abs(m["z"]))
m["padj"] = m.groupby("cell_type")["pvalue"].transform(lambda p: false_discovery_control(p.to_numpy(), method="bh"))
# Equivalence: two one-sided tests (Schuirmann 1987; Lakens et al. 2018), margin +-1 log2FC, p = the larger of the two
# one-sided p, BH within cell type over all genes tested in both diseases; computed in R (pnorm, p.adjust), as the figures
# were drawn. Classes, for the genes changed vs HC in either disease (the DEG rule: adjusted p < 0.05 and |log2FC| >= 1):
# Higher in SLE or Higher in SjD, adjusted p < 0.05 and |diff| >= 1; Equivalent, adjusted equivalence p < 0.05;
# Inconclusive (Lakens et al. 2018), neither.
EQ_MARGIN = 1
to_r("ind", m[["cell_type", "diff", "se"]])
ro.r(f"""ind$p_tost <- pmax(pnorm((ind$diff + {EQ_MARGIN}) / ind$se, lower.tail = FALSE), pnorm((ind$diff - {EQ_MARGIN}) / ind$se))
ind$padj_tost <- ave(ind$p_tost, ind$cell_type, FUN = function(p) p.adjust(p, "BH"))""")
tost = from_r("ind")
m["p_tost"], m["padj_tost"] = tost["p_tost"].to_numpy(), tost["padj_tost"].to_numpy()
for d in ("sle", "sjd"):
    m[f"deg_{d}"] = (m[f"adj_p_value_{d}"] < 0.05) & (m[f"log_fc_{d}"].abs() >= 1)
m["class"] = np.select([~(m.deg_sle | m.deg_sjd), (m.padj < 0.05) & (m["diff"] >= 1), (m.padj < 0.05) & (m["diff"] <= -1),
                        m.padj_tost < 0.05], ["", "Higher in SLE", "Higher in SjD", "Equivalent"], "Inconclusive")
m.to_csv(out / "indirect.tsv.gz", sep="\t", index=False)
log(f"indirect comparison: {len(m)} gene x cell type pairs, {int((m['padj'] < 0.05).sum())} with padj < 0.05, "
    f"{int(((m['padj'] < 0.05) & (m['diff'].abs() >= 1)).sum())} also with |diff| >= 1; classes of the genes changed vs HC: "
    + ", ".join(f"{k} {v}" for k, v in m.loc[m["class"] != "", "class"].value_counts().items()))

# --- fgsea on the indirect z statistics (settings of 06_pathways.py) ---
ro.r("suppressMessages(library(fgsea))")
hallmark = dc.pp.read_gmt(PIPELINE_DIR / "reference" / "h.all.v2025.1.Hs.symbols.gmt")
ro.globalenv["pathways"] = ro.ListVector({g: ro.StrVector(t.tolist()) for g, t in hallmark.groupby("source")["target"]})
res_all = []
for ct_, g in m.groupby("cell_type"):
    stats = ro.FloatVector(g["z"].to_numpy(dtype=float))
    stats.names = ro.StrVector(g["variable"].astype(str).tolist())
    ro.globalenv["stats"] = stats
    ro.r("set.seed(42); res <- as.data.frame(fgsea(pathways = pathways, stats = stats))[, c('pathway', 'pval', 'NES', 'log2err')]")
    res_all.append(from_r("res").assign(cell_type=ct_))
gs = pd.concat(res_all).rename(columns={"pathway": "gene_set", "NES": "score", "pval": "pvalue"}).dropna(subset=["pvalue"])
gs["padj"] = false_discovery_control(gs["pvalue"].to_numpy(), method="bh")
gs[["cell_type", "gene_set", "score", "pvalue", "padj", "log2err"]].sort_values(["padj", "cell_type"]).to_csv(
    out / "indirect_gsea.tsv", sep="\t", index=False)
log(f"fgsea on the indirect z: {len(gs)} tests, {int((gs['padj'] < 0.05).sum())} with padj < 0.05")
log("done")
