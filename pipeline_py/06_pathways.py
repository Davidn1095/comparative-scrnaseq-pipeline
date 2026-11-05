# Gene set enrichment analysis (GSEA) of the differential expression between conditions, per cell type.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19):
#   jupyter-book/conditions/gsea_pathway.ipynb, setup cells 15-16, gene-set cells 25 and 27, and cells 45-52
#   ("Population contrast enrichment scoring", GSEA part). AUCell, ULM and PROGENy are not ported.
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason. Plots
# are written to files instead of displayed.
# Gene sets: MSigDB Hallmark, release 2025.1 (pipeline_py/reference/h.all.v2025.1.Hs.symbols.gmt) instead of Reactome.
# GSEA: fgsea (Bioconductor; fgseaMultilevel with its defaults, through rpy2) instead of decoupler's dc.mt.gsea.
# Multiple testing: Benjamini-Hochberg across all cell type x gene set tests of each contrast together, instead of within
# each cell type.
# Usage: python 06_pathways.py   (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/04_differential_expression/<contrast>/all/de_all_cell_types.tsv
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/06_pathways/<contrast>/: gsea.tsv (cell type, gene set, score, p-value,
#   adjusted p-value, log2err), figures/barplot_<cell type>.png; 06_pathways/fgsea_settings.txt (fgsea version, seed)

# --- cell 15 ---
import warnings

warnings.filterwarnings(
    "ignore",
    message=r".*The method obs_vector is deprecated.*",
    category=FutureWarning,
)

warnings.filterwarnings(
    "ignore",
    message=r".*Tight layout not applied.*",
    category=UserWarning,
)

# --- cell 16 ---
import logging  # CHANGED: R messages through rpy2
import os  # CHANGED: run folder

import decoupler as dc
import matplotlib  # CHANGED: plots are written to files

matplotlib.use("Agg")  # CHANGED: as above
import matplotlib.pyplot as plt
import pandas as pd  # CHANGED: tables
import scanpy as sc
import rpy2.rinterface_lib.callbacks as rcb  # CHANGED: fgsea through rpy2, as 01_quality_control.py runs DecontX and scDblFinder
import rpy2.robjects as ro  # CHANGED: as above
from rpy2.robjects import pandas2ri  # CHANGED: as above
from rpy2.robjects.conversion import localconverter  # CHANGED: as above
from scipy.stats import false_discovery_control  # CHANGED: Benjamini-Hochberg across all tests of a contrast

from config import OUT_ROOT, PIPELINE_DIR, log  # CHANGED: pipeline paths

# CHANGED: lamindb (ln.connect, ln.track) and squidpy not imported: squidpy is not used by the GSEA part

sc.set_figure_params(figsize=(3, 3), frameon=False)

rcb.logger.setLevel(logging.ERROR)  # CHANGED: as in 01_quality_control.py
ro.r("suppressMessages(library(fgsea))")  # CHANGED: fgsea instead of dc.mt.gsea, see cell 50
FGSEA_SEED = 42  # CHANGED: seed of R's random number generator, set before every fgsea call (see cell 50)

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline input
out_root = run / "06_pathways"  # CHANGED: pipeline output


def safe(name):  # CHANGED: cell-type names as file names
    return "".join(c if c.isalnum() or c in "-." else "_" for c in str(name))


# --- cell 25 ---
# CHANGED: MSigDB Hallmark gene sets, release 2025.1, from a file in the repository instead of the Reactome set (v7.5.1) in
# lamindb
gmt = PIPELINE_DIR / "reference" / "h.all.v2025.1.Hs.symbols.gmt"  # CHANGED: as above

hallmark = dc.pp.read_gmt(gmt)  # CHANGED: as above (reactome in the notebook)

# --- cell 27 ---
geneset_size = hallmark.groupby("source").size()  # CHANGED: hallmark (reactome in the notebook)
hallmark = hallmark.loc[  # CHANGED: as above
    [15 < geneset_size.loc[gset] < 500 for gset in hallmark["source"]]  # CHANGED: as above
]
log(f"Hallmark 2025.1: {hallmark['source'].nunique()} gene sets of 15-500 genes")  # CHANGED: log
ro.globalenv["pathways"] = ro.ListVector({gset: ro.StrVector(genes.tolist()) for gset, genes in hallmark.groupby("source")["target"]})  # CHANGED: gene sets as an R list for fgsea
fgsea_version = ro.r('as.character(packageVersion("fgsea"))')[0]  # CHANGED: recorded with the seed
out_root.mkdir(parents=True, exist_ok=True)  # CHANGED: as above
(out_root / "fgsea_settings.txt").write_text(f"fgsea {fgsea_version}, fgseaMultilevel with default arguments\nset.seed({FGSEA_SEED}) before each call (one call per contrast x cell type)\np-values: fgsea's raw pval; Benjamini-Hochberg across all cell type x gene set tests of each contrast\n")  # CHANGED: as above
log(f"fgsea {fgsea_version}, seed {FGSEA_SEED}")  # CHANGED: log

for contrast in ["SLE_vs_control", "SjD_vs_control"]:  # CHANGED: the two contrasts of 04_differential_expression.py
    # --- cell 46 ---
    res_df = pd.read_csv(run / "04_differential_expression" / contrast / "all" / "de_all_cell_types.tsv", sep="\t")  # CHANGED: PyDESeq2 results of 04_differential_expression.py, all cell types

    # --- cell 48 ---
    # CHANGED: one row per cell type (the notebook has one cell type, "stim.vs.ctrl"); genes without a statistic are left out
    data = res_df.pivot_table(index="cell_type", columns="variable", values="stat")  # CHANGED: as above
    out = out_root / contrast  # CHANGED: pipeline output
    fig_dir = out / "figures"  # CHANGED: as above
    fig_dir.mkdir(parents=True, exist_ok=True)  # CHANGED: as above

    # --- cell 50 ---
    # CHANGED: fgsea (Bioconductor) instead of decoupler's dc.mt.gsea. The book names fgsea as the fast implementation of
    # GSEA (Korotkevich et al. 2021, bioRxiv 060012). dc.mt.gsea (decoupler 2.2.0) computes p-values from 1,000
    # permutations as b/n against same-sign nulls, with no +1, so the smallest non-zero p-value is about 0.002 and a
    # p-value can be 0: too coarse for Benjamini-Hochberg over about 1,250 tests per contrast. fgsea's default,
    # fgseaMultilevel, estimates small p-values to arbitrary precision (down to eps = 1e-50).
    # fgsea is run with its defaults, once per cell type with all gene sets, ranked by the PyDESeq2 Wald statistic of that
    # cell type's own genes (those kept by the expression filters of 04_differential_expression.py). Its raw p-values are
    # used; its padj (within one call, i.e. within the cell type) is not. Benjamini-Hochberg is applied once to all p-values
    # of the contrast. R's seed is set before each call, so each cell type's result does not depend on the order of calls.
    scores, pvalues, log2errs = {}, {}, {}  # CHANGED: as above
    for cell_type in data.index:  # CHANGED: as above
        row = data.loc[cell_type].dropna()  # CHANGED: as above
        stats = ro.FloatVector(row.to_numpy(dtype=float))  # CHANGED: named numeric vector for fgsea
        stats.names = ro.StrVector(row.index.astype(str).tolist())  # CHANGED: as above
        ro.globalenv["stats"] = stats  # CHANGED: as above
        ro.r(f"set.seed({FGSEA_SEED}); res <- as.data.frame(fgsea(pathways = pathways, stats = stats))"
             "[, c('pathway', 'pval', 'NES', 'log2err')]")  # CHANGED: as above
        with localconverter(ro.default_converter + pandas2ri.converter):  # CHANGED: as above
            res = ro.conversion.rpy2py(ro.globalenv["res"])  # CHANGED: as above
        for r in res.itertuples(index=False):  # CHANGED: as above
            scores[(cell_type, r.pathway)] = r.NES  # CHANGED: normalised enrichment score, as dc.mt.gsea's score
            pvalues[(cell_type, r.pathway)] = r.pval  # CHANGED: as above
            log2errs[(cell_type, r.pathway)] = r.log2err  # CHANGED: fgsea's estimate of the p-value's standard deviation (log2)
    score = pd.Series(scores).unstack()  # CHANGED: cell type x gene set
    pvalue = pd.Series(pvalues).unstack()  # CHANGED: as above
    flat = pvalue.stack().dropna()  # CHANGED: as above (gene sets fgsea leaves out for a cell type are absent)
    padj = pd.Series(false_discovery_control(flat.values, method="bh"), index=flat.index).unstack()  # CHANGED: as above
    table = pd.concat({"score": score.stack().dropna(), "pvalue": flat, "padj": padj.stack().dropna(),
                       "log2err": pd.Series(log2errs)}, axis=1)  # CHANGED: one table of all tests
    table.index.names = ["cell_type", "gene_set"]  # CHANGED: as above
    table.reset_index().sort_values(["padj", "cell_type"]).to_csv(out / "gsea.tsv", sep="\t", index=False)  # CHANGED: as above

    for cell_type in score.index:  # CHANGED: every cell type (one contrast row in the notebook)
        mask = (padj.loc[[cell_type]].T < 0.05).iloc[:, 0].fillna(False)  # CHANGED: this cell type's row
        if not mask.any():  # CHANGED: nothing to plot
            continue  # CHANGED: as above
        score_ct = score.loc[[cell_type], mask]  # CHANGED: as above

        # --- cell 52 ---
        fig = dc.pl.barplot(data=score_ct, name=cell_type, top=15, return_fig=True)  # CHANGED: this cell type

        ax = fig.axes[0]
        ax.tick_params(axis="y", labelsize=8)
        fig.savefig(fig_dir / f"barplot_{safe(cell_type)}.png", dpi=150, bbox_inches="tight")  # CHANGED: saved to file (plt.show() in the notebook)
        plt.close("all")  # CHANGED: as above
    log(f"{contrast}: {len(table)} cell type x gene set tests; adjusted p < 0.05 in {int((table['padj'] < 0.05).sum())}")  # CHANGED: log
log("done")  # CHANGED: log
