# Cell-type annotation: marker genes and differentially expressed genes per Leiden cluster, and CellTypist.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19):
#   jupyter-book/cellular_structure/annotation.ipynb, sections "From markers to cluster annotation",
#   "From cluster differentially expressed genes to cluster annotation" and "Classifiers based on a wider set of genes"
#   (cells 5-110 and 215). "Annotation by mapping to a reference" (a bone-marrow reference; CellTypist covers automated
#   annotation here) and the large-language-model section are not ported.
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason. Plots
# are written to files instead of displayed.
# Usage: python 03_annotation.py            (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
#        python 03_annotation.py --apply    (only applies the filled-in manual annotation; no recomputation)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/atlas.h5ad (02_integration.py: scVI neighbours and UMAP, X = log1p of
#   counts normalised by library size, layer "counts")
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/03_annotation/
#   cells.tsv.gz                  per cell: leiden_1, leiden_2, CellTypist labels and confidence (coarse, fine), UMAP
#   dea_leiden_{1,2}.tsv.gz, dea_leiden_{1,2}_filtered.tsv.gz   Wilcoxon results per cluster
#   figures/*.png                 UMAPs, marker dotplots, DE dotplots, dendrograms
#   manual_annotation.tsv         template (cluster, label, reason) to fill in; never overwritten
#   cells_manual_annotation.tsv.gz  written by --apply once manual_annotation.tsv has labels

# --- cell 5 ---
import os  # CHANGED: run folder, model folder and the --apply mode
import sys  # CHANGED: as above (the notebook imports it for the reference-mapping section)
from pathlib import Path  # CHANGED: kept from the notebook; used for output paths

from config import DEV_ROOT, OUT_ROOT, PIPELINE_DIR, REF_DIR, curated_markers, log  # CHANGED: pipeline paths

os.environ["CELLTYPIST_FOLDER"] = str(DEV_ROOT / "reference" / "celltypist")  # CHANGED: models stored once, read at import

import celltypist

# CHANGED: lamindb, shutil, pandas_indexes_base, scvi, cellmapper and csr_matrix not imported: used only by the
# reference-mapping section, which is not ported
import matplotlib

matplotlib.use("Agg")  # CHANGED: plots are written to files
import matplotlib.pyplot as plt
import anndata as ad  # CHANGED: make_index_unique for the CellTypist gene names; the plot object below
import scipy.sparse as sp  # CHANGED: markers kept apart
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns  # noqa: F401  (imported by the notebook; used only by the reference-mapping section)
from celltypist import models

# CHANGED: ln.connect and ln.track removed (no lamindb)

# --- cell 6 ---
# filter out some deprecation and performance warnings that do not affect our code

import warnings

# celltypist switches off all warnings when it is imported, so we restore Python's defaults
warnings.resetwarnings()
for category in (
    DeprecationWarning,
    PendingDeprecationWarning,
    ImportWarning,
    ResourceWarning,
):
    warnings.simplefilter("ignore", category=category)

from numba.core.errors import NumbaDeprecationWarning

warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# --- cell 9 ---
sc.set_figure_params(figsize=(5, 5))

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: input is the output of 02_integration.py
out = run / "03_annotation"  # CHANGED: tables and figures instead of an annotated h5ad
fig = out / "figures"  # CHANGED: as above
fig.mkdir(parents=True, exist_ok=True)  # CHANGED: as above
TEMPLATE = out / "manual_annotation.tsv"  # CHANGED: manual labels come from a file filled in by hand
# --apply reads the labels of the published atlas (MANUAL_ANNOTATION overrides it, e.g. with a filled-in TEMPLATE)
ANNOTATION = Path(os.environ.get("MANUAL_ANNOTATION", PIPELINE_DIR / "inputs" / "manual_annotation.tsv"))


def save(name):  # CHANGED: plots are written to files instead of displayed
    plt.savefig(fig / name, dpi=150, bbox_inches="tight")
    plt.close("all")


def apply_manual_annotation():
    # CHANGED: replaces cells 54, 57, 74 and 75 (a cluster -> label dictionary typed into the notebook): the labels come
    # from manual_annotation.tsv (cluster, label, reason), with clusters written as leiden_1:<id> or leiden_2:<id>. A
    # leiden_2 entry takes precedence over the leiden_1 entry of the same cells, so a resolution-1 cluster can be split.
    cells = pd.read_csv(out / "cells.tsv.gz", sep="\t", dtype={"leiden_1": str, "leiden_2": str}, index_col=0)
    t = pd.read_csv(ANNOTATION, sep="\t", dtype=str, keep_default_na=False)
    t = t[t.label.str.strip() != ""]
    if t.empty:
        log(f"{ANNOTATION} has no labels yet; nothing applied")
        return
    res = t.cluster.str.split(":", n=1).str[0]
    bad = sorted(set(res) - {"leiden_1", "leiden_2"})
    if bad:
        sys.exit(f"{ANNOTATION}: clusters must be written as leiden_1:<id> or leiden_2:<id>; found {bad}")
    cl_annotation = {k: dict(zip(t.cluster[res == k].str.split(":", n=1).str[1], t.label[res == k])) for k in ("leiden_1", "leiden_2")}
    # --- cell 57 / 75 ---
    manual = cells.leiden_1.map(cl_annotation["leiden_1"])  # CHANGED: cells table instead of adata.obs; see above
    fine = cells.leiden_2.map(cl_annotation["leiden_2"])  # CHANGED: resolution-2 entries override, see above
    manual = fine.where(fine.notna(), manual)  # CHANGED: as above
    src = np.where(fine.notna(), "leiden_2:" + cells.leiden_2, np.where(manual.notna(), "leiden_1:" + cells.leiden_1, ""))
    o = pd.DataFrame({"manual_celltype_annotation": manual, "cluster": src}, index=cells.index)
    o.to_csv(out / "cells_manual_annotation.tsv.gz", sep="\t")
    log(f"manual annotation applied to {manual.notna().sum()} of {len(o)} cells; unlabelled {manual.isna().sum()}")
    # --- cell 58 ---
    a = sc.AnnData(obs=o.assign(manual_celltype_annotation=o.manual_celltype_annotation.astype("category")))  # CHANGED: plot from the cells table
    a.obsm["X_umap"] = cells[["umap_1", "umap_2"]].to_numpy()  # CHANGED: as above
    sc.pl.umap(a, color=["manual_celltype_annotation"], show=False)  # CHANGED: show=False, saved below
    save("umap_manual_celltype_annotation.png")  # CHANGED: saved to file


if "--apply" in sys.argv:  # CHANGED: apply the filled-in manual annotation without recomputing anything
    apply_manual_annotation()
    sys.exit(0)

# --- cell 12 ---
adata = sc.read_h5ad(run / "atlas.h5ad")  # CHANGED: the output of 02_integration.py instead of the lamindb artifact
del adata.layers["logcounts"]  # CHANGED: identical to X (memory); X holds log1p of library-size-normalised counts
log(f"atlas: {adata.n_obs} cells, {adata.n_vars} genes")  # CHANGED: log

# --- cell 19 ---
# keys: cell types or populations, values: lists of marker genes
# CHANGED: the notebook's bone-marrow marker dictionary is replaced by the curated markers of the CellTypist
# Immune_All_Low cell types (CellTypist Pan Immune Atlas v2), grouped below by their Immune_All_High cell type
markers_table, marker_genes = curated_markers()  # CHANGED: as above (config.curated_markers; two symbols corrected)
high_of_low = dict(zip(markers_table["Low-hierarchy cell types"], markers_table["High-hierarchy cell types"]))  # CHANGED: as above

# CHANGED: curated markers not measured by every dataset (written by 02_integration.py); their values
# come only from cells of the datasets that measured them, and their labels say how many datasets did
mx = sc.read_h5ad(run / "markers_unshared.h5ad")
if not mx.obs_names.equals(adata.obs_names):
    sys.exit("markers_unshared.h5ad: cells differ from the atlas")
n_ds = len(mx.uns["datasets"])
mx_label = {g: f"{g} ({int(n)}/{n_ds})" for g, n in mx.var.n_datasets.items()}
mx_cells = {g: adata.obs["dataset"].isin([d for d in mx.uns["datasets"] if mx.var.loc[g, f"measured_{d}"]]).values for g in mx.var_names}
log(f"markers outside the shared genes: {', '.join(mx_label.values())}")


_shared_m = sorted({m for v in marker_genes.values() for m in v} & set(adata.var_names))
_mcols = {"atlas": (sp.csc_matrix(adata[:, _shared_m].X), {g: i for i, g in enumerate(_shared_m)}),  # one slice of all markers
          "apart": (sp.csc_matrix(mx.X), {g: i for i, g in enumerate(mx.var_names)})}


def dot_tables(res, genes):  # CHANGED: mean (of X, log1p normalised) and fraction expressing per cluster, as sc.pl.dotplot
    """computes them, for a shared gene over all cells and for a marker kept apart over the cells that measured it."""
    groups = adata.obs[res].cat.categories
    color, size = pd.DataFrame(index=groups, dtype=float), pd.DataFrame(index=groups, dtype=float)
    for g in genes:
        src, keep, lab = ("atlas", np.ones(adata.n_obs, bool), g) if g in adata.var_names else ("apart", mx_cells[g], mx_label[g])
        M, col = _mcols[src]
        v = M[:, col[g]].toarray().ravel()[keep]
        c = adata.obs[res].values[keep]
        d = pd.DataFrame({"v": v, "c": c})
        color[lab] = d.groupby("c", observed=False).v.mean().reindex(groups)
        size[lab] = d.assign(e=d.v > 0).groupby("c", observed=False).e.mean().reindex(groups)
    return color, size


# --- cell 22 ---
marker_genes_in_data = {}
for ct, markers in marker_genes.items():
    markers_found = []
    for marker in markers:
        if marker in adata.var.index or marker in mx.var_names:  # CHANGED: or measured by some datasets (kept apart)
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
st = pd.read_csv(run / "curated_markers_status.tsv", sep="\t", index_col=0).status  # CHANGED: written by 02_integration.py
for k in ["excluded_decision_2", "under_20_cells", "not_measured"]:  # CHANGED: log, by reason
    log(f"curated markers not used, {k.replace('_', ' ')}: " + ", ".join(st.index[st == k]))

# --- cell 25 ---
# CHANGED: adata.X = adata.layers["scran_normalization"] not needed: X already holds the log1p of library-size-normalised
# counts computed by 02_integration.py, used for the marker plots and the differential expression
# --- cells 27, 29, 31 and 33 (HVG, PCA, neighbours, UMAP) left out ---
# CHANGED: the neighbours (on the scVI latent space) and the UMAP of 02_integration.py are used instead of a PCA of
# scran-normalised data

# --- cell 36 ---
# CHANGED: the notebook picks the B and plasma cell types; here every Immune_All_High cell type is a group, and each group's
# Immune_All_Low cell types are plotted together (cells 51 and 52, below)
groups = {}  # CHANGED: as above
for ct, hi in high_of_low.items():  # CHANGED: as above
    groups.setdefault(hi, []).append(ct)  # CHANGED: as above

# --- cell 38 (one UMAP per marker gene of each cell type) left out ---
# CHANGED: with about 100 cell types this would be several hundred UMAPs; the per-cluster dotplots and marker scores below
# carry the same information

# --- cell 42 ---
sc.tl.leiden(adata, resolution=1, key_added="leiden_1", flavor="igraph", n_iterations=2)
# --- cell 43 ---
sc.pl.umap(adata, color="leiden_1", show=False)  # CHANGED: show=False, saved below
save("umap_leiden_1.png")  # CHANGED: saved to file
# --- cell 45 ---
sc.tl.leiden(adata, resolution=2, key_added="leiden_2", flavor="igraph", n_iterations=2)
# --- cell 47 ---
sc.pl.umap(adata, color="leiden_2", legend_loc="on data", show=False)  # CHANGED: show=False, saved below
save("umap_leiden_2.png")  # CHANGED: saved to file
# --- cell 49 ---
sc.pl.umap(adata, color="leiden_1", legend_loc="on data", show=False)  # CHANGED: show=False, saved below
save("umap_leiden_1_on_data.png")  # CHANGED: saved to file
log(f"Leiden: {adata.obs.leiden_1.nunique()} clusters at resolution 1, {adata.obs.leiden_2.nunique()} at resolution 2")  # CHANGED: log

for res in ["leiden_1", "leiden_2"]:  # CHANGED: dotplots per cluster at both resolutions
    for hi, cts in groups.items():  # CHANGED: one dotplot per Immune_All_High cell type, see cell 36
        # --- cell 51 ---
        group_markers = {
            ct: [m for m in ct_markers if m in adata.var.index or m in mx.var_names]  # CHANGED: or kept apart
            for ct, ct_markers in marker_genes.items()
            if ct in cts  # CHANGED: the group's cell types instead of B_plasma_cts
        }
        group_markers = {ct: m for ct, m in group_markers.items() if m}  # CHANGED: drop cell types without markers in the data
        if not group_markers:  # CHANGED: as above
            continue  # CHANGED: as above
        # --- cell 52 ---
        # CHANGED: sc.pl.dotplot(adata, groupby=res, var_names=group_markers, standard_scale="var") is drawn from the tables
        # of dot_tables, so that a marker kept apart uses only the cells that measured it; standard_scale="var" (each gene
        # scaled to 0-1 across clusters) is applied to the colour table as scanpy does
        genes = list(dict.fromkeys(m for ms in group_markers.values() for m in ms))
        color, size = dot_tables(res, genes)
        color = (color - color.min()) / (color.max() - color.min()).replace(0, 1)
        labelled = {ct: [m if m in adata.var_names else mx_label[m] for m in ms] for ct, ms in group_markers.items()}
        plot = ad.AnnData(X=sp.csr_matrix((adata.n_obs, len(color.columns))), obs=adata.obs[[res]],
                          var=pd.DataFrame(index=color.columns))  # carries only the clusters and gene labels for the plot
        cols = [m for ms in labelled.values() for m in ms]  # a marker of several cell types is drawn once per cell type, as
        sc.pl.DotPlot(plot, var_names=labelled, groupby=res, dot_color_df=color[cols], dot_size_df=size[cols]).make_figure()  # sc.pl.dotplot does
        save(f"dotplot_markers_{res}_{hi.replace(' ', '_').replace('/', '_')}.png")  # CHANGED: saved to file
        tab = pd.concat({"mean_scaled": color, "fraction": size}, axis=1).stack(level=1, future_stack=True)  # CHANGED: values kept
        tab.rename_axis(["cluster", "gene"]).to_csv(out / f"markers_{res}_{hi.replace(' ', '_').replace('/', '_')}.tsv", sep="\t")

# --- cells 54, 57 and 58 (manual cluster -> label dictionary): replaced by manual_annotation.tsv, see apply_manual_annotation ---

# CHANGED: X stored as CSC instead of CSR, storage format only (same values): the Wilcoxon test takes
# a few genes at a time over all cells, a column slice that scans the whole CSR matrix each time; nothing after the
# differential expression reads X (CellTypist uses the counts layer, the dendrograms X_pca)
adata.X = sp.csc_matrix(adata.X)

for res in ["leiden_1", "leiden_2"]:  # CHANGED: differential expression at both resolutions (the notebook uses leiden_1)
    # --- cell 62 ---
    sc.tl.rank_genes_groups(
        adata, groupby=res, method="wilcoxon", key_added=f"dea_{res}"  # CHANGED: both resolutions
    )
    # --- cell 64 ---
    sc.tl.dendrogram(
        adata,
        groupby=res,  # CHANGED: both resolutions
    )

    sc.pl.rank_genes_groups_dotplot(
        adata, groupby=res, standard_scale="var", n_genes=5, key=f"dea_{res}", show=False  # CHANGED: both resolutions; saved below
    )
    save(f"dotplot_dea_{res}.png")  # CHANGED: saved to file
    # --- cell 66 ---
    sc.tl.filter_rank_genes_groups(
        adata,
        min_in_group_fraction=0.2,
        max_out_group_fraction=0.2,
        key=f"dea_{res}",  # CHANGED: both resolutions
        key_added=f"dea_{res}_filtered",  # CHANGED: as above
    )
    # --- cell 68 ---
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=res,  # CHANGED: both resolutions
        standard_scale="var",
        n_genes=5,
        key=f"dea_{res}_filtered",  # CHANGED: as above
        show=False,  # CHANGED: saved below
    )
    save(f"dotplot_dea_{res}_filtered.png")  # CHANGED: saved to file
    for key in (f"dea_{res}", f"dea_{res}_filtered"):  # CHANGED: DE tables written out (the notebook keeps them in adata.uns)
        sc.get.rank_genes_groups_df(adata, group=None, key=key).to_csv(out / f"{key}.tsv.gz", sep="\t", index=False)  # CHANGED: as above
    log(f"{res}: differential expression done")  # CHANGED: log

# --- cell 70 (UMAP of MYOM2, KLRD1, KLRF1: the notebook's NK example) left out ---
# CHANGED: example specific to the notebook's dataset
# --- cells 74 and 75: replaced by manual_annotation.tsv, see apply_manual_annotation ---

# --- cell 87 ---
adata_celltypist = sc.AnnData(  # CHANGED: counts, cell metadata and the neighbourhood graph only, instead of adata.copy():
    X=adata.layers["counts"].copy(),  # a full copy of the atlas would double its memory. The graph is kept as in the
    obs=adata.obs[[]].copy(),  # notebook's copy, so CellTypist's majority voting over-clusters on it (here the scVI
    var=adata.var[[]].copy(),  # neighbours) instead of building a new one
    obsp={"connectivities": adata.obsp["connectivities"], "distances": adata.obsp["distances"]},
    uns={"neighbors": adata.uns["neighbors"]},
)  # set adata.X to raw counts
sc.pp.normalize_total(
    adata_celltypist, target_sum=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# CHANGED: .X is kept sparse: dense it would take about 126 GB for this atlas, and CellTypist 1.7.1 accepts sparse input

# --- cell 89 ---
# CHANGED: models.download_models(force_update=True, ...) removed: compute nodes have no internet; Immune_All_Low and
# Immune_All_High (v2, 2022-07-16) were downloaded once on the login node into CELLTYPIST_FOLDER
missing = [m for m in ["Immune_All_Low.pkl", "Immune_All_High.pkl"] if not (Path(os.environ["CELLTYPIST_FOLDER"]) / "data" / "models" / m).exists()]  # CHANGED: as above
if missing:  # CHANGED: as above
    sys.exit(f"CellTypist models missing in {os.environ['CELLTYPIST_FOLDER']}: {missing}; download them on the login node")  # CHANGED: as above

# --- cell 91 ---
model_low = models.Model.load(model="Immune_All_Low.pkl")
model_high = models.Model.load(model="Immune_All_High.pkl")


def map_model_genes(model):  # CHANGED: the models name genes by older symbols; mapped with the rule of the data
    """Rename the model's genes to the atlas gene names: symbol -> its GRCh38-2020-A Ensembl ID -> the reference gene with
    that ID; otherwise the reference gene whose symbol is the same and unique, if no other model gene reached it."""
    gtf = pd.read_csv(REF_DIR / "genes" / "genes.gtf.gz", sep="\t", comment="#", header=None, usecols=[2, 8], dtype=str)
    gtf = gtf[gtf[2] == "gene"][8]
    ref = pd.DataFrame({"gene_id": gtf.str.extract(r'gene_id "([^"]+)"')[0].values,
                        "gene_name": gtf.str.extract(r'gene_name "([^"]+)"')[0].values})
    atlas_name = pd.Series(ad.utils.make_index_unique(pd.Index(ref.gene_name)).values, index=ref.gene_id.values)  # as step 01
    f20 = pd.read_csv(PIPELINE_DIR / "reference" / "GRCh38-2020-A_features.tsv.gz", sep="\t", header=None, dtype=str)
    s20 = f20[f20[1].map(f20[1].value_counts()) == 1].set_index(1)[0]
    s24 = ref[ref.gene_name.map(ref.gene_name.value_counts()) == 1].set_index("gene_name").gene_id
    feat = pd.Series(model.features)
    t = feat.map(s20).where(lambda x: x.isin(atlas_name.index))
    sym = feat.map(s24)
    cand = t.isna() & sym.notna() & ~sym.isin(t.dropna()) & ~sym.duplicated(keep=False)
    t[cand] = sym[cand]
    new = np.where(t.notna(), atlas_name.reindex(t.fillna("")).values, "unmapped:" + feat.values)  # unmapped: matches no gene
    model.classifier.features = new
    log(f"CellTypist genes mapped: {int(t.notna().sum() - cand.sum())} by 2020-A ID, {int(cand.sum())} by symbol, "
        f"{int(t.isna().sum())} unmapped, of {len(feat)}")
    return model


model_low, model_high = map_model_genes(model_low), map_model_genes(model_high)  # CHANGED: see map_model_genes

# --- cell 93 ---
log("Immune_All_High cell types:", len(model_high.cell_types))  # CHANGED: logged instead of displayed
# --- cell 94 ---
log("Immune_All_Low cell types:", len(model_low.cell_types))  # CHANGED: logged instead of displayed

# --- cell 97 ---
predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
)

# --- cell 99 ---
predictions_high_adata = predictions_high.to_adata()

# --- cell 101 ---
adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "conf_score"
]
del predictions_high, predictions_high_adata  # CHANGED: memory

# --- cell 103 ---
predictions_low = celltypist.annotate(
    adata_celltypist, model=model_low, majority_voting=True
)
# --- cell 104 ---
predictions_low_adata = predictions_low.to_adata()
# --- cell 105 ---
adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "conf_score"
]
del predictions_low, predictions_low_adata, adata_celltypist  # CHANGED: memory
log("CellTypist done")  # CHANGED: log

# --- cell 107 ---
sc.pl.umap(
    adata,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
    show=False,  # CHANGED: saved below
)
save("umap_celltypist_coarse.png")  # CHANGED: saved to file

# --- cell 108 ---
sc.pl.umap(
    adata,
    color=["celltypist_cell_label_fine", "celltypist_conf_score_fine"],
    frameon=False,
    sort_order=False,
    wspace=1,
    show=False,  # CHANGED: saved below
)
save("umap_celltypist_fine.png")  # CHANGED: saved to file

# --- cell 110 ---
sc.tl.dendrogram(
    adata,
    groupby="celltypist_cell_label_fine",
)

sc.pl.dendrogram(adata, groupby="celltypist_cell_label_fine", show=False)  # CHANGED: saved below
save("dendrogram_celltypist_fine.png")  # CHANGED: saved to file

# --- cell 215 ---
# CHANGED: the notebook formats the DE names for h5ad and saves the annotated object to lamindb. The atlas is not written
# again (about 40 GB); the cell-level annotation is written as a table, with the UMAP for plotting
cells = adata.obs[["leiden_1", "leiden_2", "celltypist_cell_label_coarse", "celltypist_conf_score_coarse",
                   "celltypist_cell_label_fine", "celltypist_conf_score_fine"]].copy()  # CHANGED: as above
cells[["umap_1", "umap_2"]] = adata.obsm["X_umap"]  # CHANGED: as above
cells.to_csv(out / "cells.tsv.gz", sep="\t")  # CHANGED: as above

if not TEMPLATE.exists():  # CHANGED: template for the manual labels (cells 54-75); never overwritten
    pd.DataFrame({"cluster": [f"leiden_1:{c}" for c in sorted(adata.obs.leiden_1.cat.categories, key=int)],
                  "label": "", "reason": ""}).to_csv(TEMPLATE, sep="\t", index=False)
    log(f"template for manual annotation written: {TEMPLATE}")
apply_manual_annotation()
log("done")
