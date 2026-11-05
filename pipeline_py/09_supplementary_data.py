# Inputs of the supplementary figures that the step outputs do not hold, drawn by 10_supplementary.R. Not a
# pipeline step and not a notebook port.
# Usage: run by make_supplementary.sh (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/: atlas.h5ad (obs), 03_annotation/cells_manual_annotation.tsv.gz;
#   $PIPELINE_OUT/py/datasets/<dataset>/01_quality_control/cells.tsv (step 01's outcome of each filter, per cell)
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/09_supplementary/data/
#   composition_counts.tsv  cells per donor and cell type, with dataset and condition: the cells of step 05's scCODA model
#                           (labelled, not "low quality" or "platelets"); a cell type a donor lacks is a row with n = 0
#   qc_counts.tsv           cells per dataset after each step 01 filter, in the atlas and after the annotation exclusion
#                          : recorded counts only, the per-sample thresholds are not recorded
#   atlas_cells.tsv.gz      per atlas cell: dataset, condition, UMAP before integration (obsm X_umap_unintegrated) and the
#                           integrated UMAP (X_umap), annotation cluster and its label (supplementary atlas figure, a-c)
#   atlas_dotplot.tsv       per final cell type and marker of config.SUPP_DOTPLOT_MARKERS: mean expression (X, log1p of
#                           library-size normalised counts) and fraction of cells with counts, as sc.pl.dotplot; markers
#                           outside the shared genes from markers_unshared.h5ad, over the cells of the datasets that
#                           measured them, as 03_annotation.py (figure panel d)
import os

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from config import DATASETS, OUT_ROOT, SUPP_DOTPLOT_MARKERS, dataset_dir, excluded_labels, log

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")
out = run / "09_supplementary" / "data"
out.mkdir(parents=True, exist_ok=True)

with h5py.File(run / "atlas.h5ad", "r") as f:
    obs = ad.io.read_elem(f["obs"])[["donor_id", "sample_id", "dataset", "condition"]].copy()
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)
obs["cell_type"] = labels["manual_celltype_annotation"].reindex(obs.index)
excluded = excluded_labels(obs)

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
                 "final": int((in_ds & ~excluded & obs.cell_type.notna()).sum())})
qc = pd.DataFrame(rows)
qc.to_csv(out / "qc_counts.tsv", sep="\t", index=False)
log(f"qc_counts.tsv: {len(qc)} datasets, {qc.samples.sum()} samples, {qc.final.sum()} final cells")

# --- atlas figure, a-c: the two UMAPs and the annotation clusters ---
with h5py.File(run / "atlas.h5ad", "r") as f:
    um = {k: f["obsm"][k][:] for k in ("X_umap_unintegrated", "X_umap")}
cells = pd.DataFrame({"dataset": obs.dataset.astype(str).values, "condition": obs.condition.astype(str).values,
                      "umap_unint_1": um["X_umap_unintegrated"][:, 0], "umap_unint_2": um["X_umap_unintegrated"][:, 1],
                      "umap_1": um["X_umap"][:, 0], "umap_2": um["X_umap"][:, 1],
                      "cluster": labels["cluster"].reindex(obs.index).values, "label": obs.cell_type.values}, index=obs.index)
assert cells.cluster.notna().all(), "cells without an annotation cluster"
cells.rename_axis("cell_id").to_csv(out / "atlas_cells.tsv.gz", sep="\t")
log(f"atlas_cells.tsv.gz: {len(cells)} cells, {cells.cluster.nunique()} clusters, "
    f"{cells.cluster[excluded.values].nunique()} of them excluded at annotation")

# --- atlas figure, d: dot plot of the markers over the final cell types; X read in row blocks (the whole matrix does
# not fit in memory), only the marker columns kept ---
final = (~excluded & obs.cell_type.notna()).values
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
        v, keep, n_ds = M[:, shared.index(g)].toarray().ravel(), final, len(mx.uns["datasets"])
    else:
        meas = [d for d in mx.uns["datasets"] if mx.var.loc[g, f"measured_{d}"]]
        v = mx[:, g].X.toarray().ravel()
        keep, n_ds = final & obs.dataset.astype(str).isin(meas).values, len(meas)
    d = pd.DataFrame({"cell_type": ct[keep], "v": v[keep]})
    r = d.groupby("cell_type").v.agg(mean="mean", fraction=lambda x: (x > 0).mean()).reset_index()
    rows.append(r.assign(gene=g, n_datasets=n_ds, n_datasets_all=len(mx.uns["datasets"])))
dot = pd.concat(rows)[["gene", "cell_type", "mean", "fraction", "n_datasets", "n_datasets_all"]]
dot.to_csv(out / "atlas_dotplot.tsv", sep="\t", index=False)
log(f"atlas_dotplot.tsv: {len(genes)} markers ({len(genes) - len(shared)} outside the shared genes: "
    + ", ".join(f"{g} {int(dot[dot.gene == g].n_datasets.iloc[0])}/{len(mx.uns['datasets'])}" for g in genes if g not in shared)
    + f"), {dot.cell_type.nunique()} cell types")

obs = obs[obs["cell_type"].notna() & ~excluded]
obs = obs.drop(columns="sample_id")
for c in obs.columns:
    obs[c] = obs[c].astype(str)
counts = obs.groupby(["donor_id", "cell_type"]).size().unstack(fill_value=0).stack().rename("n").reset_index()
donors = obs.drop_duplicates("donor_id").set_index("donor_id")[["dataset", "condition"]]
counts = counts.join(donors, on="donor_id")[["donor_id", "dataset", "condition", "cell_type", "n"]]
counts.to_csv(out / "composition_counts.tsv", sep="\t", index=False)
log(f"composition_counts.tsv: {counts.donor_id.nunique()} donors, {counts.cell_type.nunique()} cell types, "
    f"{int(counts.n.sum())} cells")
