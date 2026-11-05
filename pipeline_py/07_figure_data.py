# Inputs of the manuscript figures, drawn by 08_figures.R. Not a pipeline step and
# not a notebook port.
# Usage: run by make_figures.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/: atlas.h5ad (obs, X_umap), 03_annotation/cells_manual_annotation.tsv.gz,
#        04_differential_expression/<contrast>/all/de_all_cell_types.tsv
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/07_figures/data/
#   cells.tsv.gz      UMAP coordinates, dataset, condition and cell type of the cells analysed in steps 04-06 (labelled,
#                     not "low quality" or "platelets"; Fig 2)
#   table1.tsv        per dataset, over the atlas cells (low-quality and platelet clusters removed): cases, controls,
#                     cells analysed, median genes per cell; and cells after QC before that removal (Table 1)
#   indirect.tsv.gz   SLE vs SjD indirect comparison per gene and cell type (Fig 3b): difference of the two log2FC (each
#                     disease vs its own controls), SE = sqrt(SE_SLE^2 + SE_SjD^2), Wald z, BH within cell type; genes
#                     tested in both contrasts
#   indirect_gsea.tsv fgsea on the z statistics of the indirect comparison (Fig 4b), as 06_pathways.py: Hallmark 2025.1,
#                     fgseaMultilevel defaults, set.seed(42) before each cell type, BH across all cell type x gene set tests
import os

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.stats import false_discovery_control, norm

from config import DATASETS, OUT_ROOT, PIPELINE_DIR, excluded_labels, log

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")
out = run / "07_figures" / "data"
out.mkdir(parents=True, exist_ok=True)

# --- cells: UMAP and labels ---
with h5py.File(run / "atlas.h5ad", "r") as f:
    obs = ad.io.read_elem(f["obs"])[["dataset", "condition", "donor_id", "n_genes_by_counts"]].copy()
    umap = f["obsm/X_umap"][:]
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)
obs["cell_type"] = labels["manual_celltype_annotation"].reindex(obs.index)
obs["UMAP_1"], obs["UMAP_2"] = umap[:, 0], umap[:, 1]
analysed = (obs["cell_type"].notna() & ~excluded_labels(obs)).to_numpy()

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

cells = obs.loc[analysed, ["UMAP_1", "UMAP_2", "dataset", "condition", "cell_type"]]
cells.round({"UMAP_1": 3, "UMAP_2": 3}).to_csv(out / "cells.tsv.gz", sep="\t", index_label="cell_id")
log(f"cells.tsv.gz: {len(cells)} cells analysed of {len(obs)}")

# --- SLE vs SjD indirect comparison ---
de = run / "04_differential_expression"
cols = ["variable", "cell_type", "log_fc", "lfcSE", "stat", "adj_p_value"]
sle = pd.read_csv(de / "SLE_vs_control" / "all" / "de_all_cell_types.tsv", sep="\t", usecols=cols)
sjd = pd.read_csv(de / "SjD_vs_control" / "all" / "de_all_cell_types.tsv", sep="\t", usecols=cols)
m = sle.merge(sjd, on=["cell_type", "variable"], suffixes=("_sle", "_sjd"))
m = m.dropna(subset=["log_fc_sle", "lfcSE_sle", "log_fc_sjd", "lfcSE_sjd"])
m["diff"] = m["log_fc_sle"] - m["log_fc_sjd"]
m["se"] = np.sqrt(m["lfcSE_sle"] ** 2 + m["lfcSE_sjd"] ** 2)
m["z"] = m["diff"] / m["se"]
m["pvalue"] = 2 * norm.sf(np.abs(m["z"]))
m["padj"] = m.groupby("cell_type")["pvalue"].transform(lambda p: false_discovery_control(p.to_numpy(), method="bh"))
m.to_csv(out / "indirect.tsv.gz", sep="\t", index=False)
log(f"indirect comparison: {len(m)} gene x cell type pairs, {int((m['padj'] < 0.05).sum())} with padj < 0.05, "
    f"{int(((m['padj'] < 0.05) & (m['diff'].abs() >= 1)).sum())} also with |diff| >= 1 (all genes tested in both; Fig 3b draws the DEGs vs HC among them)")

# --- fgsea on the indirect z statistics (settings of 06_pathways.py) ---
rcb.consolewrite_warnerror = lambda s: None
ro.r("suppressMessages(library(fgsea))")
hallmark = dc.pp.read_gmt(PIPELINE_DIR / "reference" / "h.all.v2025.1.Hs.symbols.gmt")
ro.globalenv["pathways"] = ro.ListVector({g: ro.StrVector(t.tolist()) for g, t in hallmark.groupby("source")["target"]})
res_all = []
for ct, g in m.groupby("cell_type"):
    stats = ro.FloatVector(g["z"].to_numpy(dtype=float))
    stats.names = ro.StrVector(g["variable"].astype(str).tolist())
    ro.globalenv["stats"] = stats
    ro.r("set.seed(42); res <- as.data.frame(fgsea(pathways = pathways, stats = stats))[, c('pathway', 'pval', 'NES', 'log2err')]")
    with localconverter(ro.default_converter + pandas2ri.converter):
        res_all.append(ro.conversion.rpy2py(ro.globalenv["res"]).assign(cell_type=ct))
gs = pd.concat(res_all).rename(columns={"pathway": "gene_set", "NES": "score", "pval": "pvalue"}).dropna(subset=["pvalue"])
gs["padj"] = false_discovery_control(gs["pvalue"].to_numpy(), method="bh")
gs[["cell_type", "gene_set", "score", "pvalue", "padj", "log2err"]].sort_values(["padj", "cell_type"]).to_csv(
    out / "indirect_gsea.tsv", sep="\t", index=False)
log(f"fgsea on the indirect z: {len(gs)} tests, {int((gs['padj'] < 0.05).sum())} with padj < 0.05")
log("done")
