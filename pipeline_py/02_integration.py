# Step 02 of the Python pipeline: sex inference and gene exclusion -> merge samples -> gene filter on the
# merged atlas -> normalisation -> HVG + scVI with the 10x library as batch, neighbours
# -> Leiden.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19):
#   jupyter-book/cellular_structure/integration.ipynb (last changed d3c6e8a, 2026-09-29), cells 11-67 (scVI part);
#   jupyter-book/cellular_structure/clustering.ipynb (last changed 893f34f, 2026-09-29), cells 12 and 15.
# The notebooks' code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason.
# Plot and display cells are left out, as marked. The book has no code for sex inference or gene exclusion: that part is
# new (donor sex from XIST and chrY genes; mitochondrial and red-cell genes removed) and marked as such.
# Usage: python 02_integration.py   (ATLAS_DATASETS=a,b,c and ATLAS_RUN=<folder> for a subset run, e.g. a code test)
# Output: $PIPELINE_OUT/py/atlas/02_integration/: atlas.h5ad (all kept genes; layers counts and logcounts; obsm X_pca,
#   X_umap_unintegrated, X_scVI, X_umap; leiden clusterings), scvi_model/, donor_sex.tsv, confounding.tsv,
#   excluded_genes.tsv, genes.tsv, markers_unshared.h5ad (curated markers not measured by every dataset, for 03),
#   curated_markers_status.tsv (whether each curated marker is in the atlas, kept apart, excluded, under 20 cells or not measured)

# --- integration cell 11 ---
import warnings

# This looks for any warning containing this specific text
warnings.filterwarnings("ignore", message=".*encoding metadata.*")
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings(
    "ignore", message=".*The default of observed=False is deprecated.*"
)

# --- integration cell 13 ---
# Python packages
# CHANGED: anndata2ri, bbknn, lamindb and scib not imported: no R calls here, scVI is the method, no lamindb,
# scIB metrics run in a later step (02b_integration_metrics.py)
import anndata as ad  # CHANGED: needed to read and merge the per-sample files of step 01
import matplotlib.pyplot as plt  # noqa: F401
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

import scipy.sparse as sp  # CHANGED: markers kept apart
from config import DATASETS, DEV_ROOT, INPUT_ROOT, OUT_ROOT, curated_markers, dataset_dir, log  # CHANGED: pipeline paths

# CHANGED: %load_ext rpy2.ipython, anndata2ri.set_ipython_converter(), ln.connect, ln.track removed (not used here)
# --- integration cell 15 (%%R library(Seurat)) left out: not used ---

scvi.settings.seed = 0  # CHANGED: the book sets no seed (scvi-tools default None); seed 0 for reproducibility
import os  # CHANGED: optional subset and output folder, for code tests
RUN_DATASETS = [d for d in os.environ.get("ATLAS_DATASETS", ",".join(DATASETS.dataset)).split(",") if d]  # CHANGED: as above
out = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline output
out.mkdir(parents=True, exist_ok=True)  # CHANGED: as above

# ===== not in the book: sex inference, then gene exclusion, and the merge of the samples =====
# Reference genes of GRCh38-2024-A with chromosome and position (written by the pipeline from the Cell Ranger GTF)
genes = pd.read_csv(DEV_ROOT / "reference" / "genes.tsv", sep="\t")
PAR_Y = [(10001, 2781479), (56887903, 57217415)]  # GRCh38 pseudoautosomal regions of chrY
in_par = np.zeros(len(genes), bool)
for s, e in PAR_Y:
    in_par |= (genes.start >= s).values & (genes.end <= e).values
RED_CELL = ["HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ", "HBM", "HBQ1", "SLC4A1", "ALAS2", "AHSP", "CA1"]
# sex-chromosome genes (XIST, TSIX, non-PAR chrY) stay in the matrix (provisional, 2026-10-27): sex enters the
# differential-expression design instead (~ dataset + sex + condition)
group = pd.Series(np.select(
    [genes.chrom.eq("chrM").values, genes.gene_name.isin(RED_CELL).values],
    ["mitochondrial", "red_cell"], default=""), index=genes.gene_id.values)
y_ids = set(genes.gene_id[(genes.chrom == "chrY").values & ~in_par])
xist_id = genes.gene_id[genes.gene_name == "XIST"].iloc[0]

parts, sex_rows = [], []
for ds in RUN_DATASETS:
    d = dataset_dir(ds, "01_quality_control") / "samples"
    if not d.exists() or not any(d.glob("*.h5ad")):
        log(f"{ds}: no step-01 output, skipped")
        continue
    a = ad.concat([ad.read_h5ad(f) for f in sorted(d.glob("*.h5ad"))], merge="same")
    a.obs["dataset"] = ds
    a.obs_names = ds + ":" + a.obs_names
    lib = a.obs["library"].astype(str)
    if a.obs["library"].isna().any() or lib.isin(["", "nan", "None"]).any():  # the integration batch
        raise SystemExit(f"{ds}: {int((a.obs['library'].isna() | lib.isin(['', 'nan', 'None'])).sum())} cells without a library ID")
    a.obs["library"] = ds + ":" + lib  # unique across datasets
    donors = pd.read_csv(INPUT_ROOT / ds / "donors.tsv", sep="\t", dtype=str)
    a.obs["condition"] = a.obs.sample_id.map(donors.drop_duplicates("sample_id").set_index("sample_id").condition).values
    a.obs["disease"] = DATASETS.set_index("dataset").disease[ds]
    if ds == "Perez":  # one sample per donor: the first visit
        s = donors[["donor_id", "sample_id", "disease_state", "processing_cohort"]].copy()
        s["rank"] = s.processing_cohort.astype(float) * 10 + s.disease_state.map({"managed": 1, "flare": 2, "treated": 3, "na": 4})
        first = s.loc[s.groupby("donor_id")["rank"].idxmin(), "sample_id"]
        log(f"Perez first visit only: {len(s) - len(first)} later samples dropped")
        a = a[a.obs.sample_id.isin(first).values].copy()
    # sex inference: per donor, XIST and non-PAR chrY genes in counts per million of the corrected counts;
    # female if XIST CPM exceeds chrY CPM
    gid = a.var.gene_ids.values
    X = a.X.tocsr()
    tot = np.asarray(X.sum(1)).ravel()
    xs = np.asarray(X[:, gid == xist_id].sum(1)).ravel()
    ys = np.asarray(X[:, np.isin(gid, list(y_ids))].sum(1)).ravel()
    per = pd.DataFrame({"donor_id": a.obs.donor_id.values, "condition": a.obs.condition.values, "total": tot, "xist": xs, "y": ys})
    per = per.groupby(["donor_id", "condition"], observed=True)[["total", "xist", "y"]].sum().reset_index()
    per["XIST_CPM"] = 1e6 * per.xist / per.total
    per["chrY_CPM"] = 1e6 * per.y / per.total
    per["inferred_sex"] = np.where(per.XIST_CPM > per.chrY_CPM, "female", "male")
    per.insert(0, "dataset", ds)
    sex_rows.append(per)
    a.obs["inferred_sex"] = a.obs.donor_id.map(per.set_index("donor_id").inferred_sex).values
    # gene exclusion: red-cell and mitochondrial genes, after sex inference
    a = a[:, group.reindex(gid).fillna("").values == ""].copy()
    a.X = a.layers["decontX_counts"]  # the ambient-corrected counts are the counts of every later step
    del a.layers["decontX_counts"]
    parts.append(a)
    log(f"{ds}: {a.n_obs} cells, {per.shape[0]} donors ({(per.inferred_sex == 'female').sum()} female)")

sex = pd.concat(sex_rows)
sex.round(2).to_csv(out / "donor_sex.tsv", sep="\t", index=False)
conf = sex.groupby(["dataset", "condition", "inferred_sex"]).size().unstack(fill_value=0).reset_index()
conf.to_csv(out / "confounding.tsv", sep="\t", index=False)  # check that disease is not confounded with donor sex
pd.DataFrame({"gene_id": group.index, "group": group.values}).query("group != ''").merge(
    genes[["gene_id", "gene_name"]]).to_csv(out / "excluded_genes.tsv", sep="\t", index=False)
# the atlas keeps the genes shared by all datasets (inner join below). Curated markers that some
# dataset did not measure (e.g. IGHM, IGKC, TRDC: the GSE135779 deposit is hg19 without immunoglobulin and T-cell receptor
# constant genes) are kept apart for annotation (03), with a record of which dataset measured each; unmeasured is not zero
shared = set.intersection(*[set(a.var_names) for a in parts])
measured_any = set.union(*[set(a.var_names) for a in parts])
extra = sorted({m for v in curated_markers()[1].values() for m in v} & (measured_any - shared))
xparts = []
for a in parts:
    have = [g for g in extra if g in a.var_names]
    idx = {g: i for i, g in enumerate(have)}
    sub = a[:, have].X.tocsc() if have else sp.csc_matrix((a.n_obs, 0))
    cols = [sub[:, idx[g]] if g in idx else sp.csc_matrix((a.n_obs, 1)) for g in extra]
    x = ad.AnnData(X=sp.hstack(cols).tocsr() if extra else sp.csr_matrix((a.n_obs, 0)),
                   obs=pd.DataFrame({"dataset": a.obs["dataset"].values}, index=a.obs_names),
                   var=pd.DataFrame({f"measured_{a.obs['dataset'].iloc[0]}": [g in idx for g in extra]}, index=extra))
    xparts.append(x)
markers_extra = ad.concat(xparts, merge="first")  # one var column per dataset: measured or not
del xparts
log(f"genes: {len(shared)} shared by all datasets, {len(measured_any)} measured by any; curated markers outside the shared "
    f"genes, kept for annotation: {len(extra)} ({', '.join(extra)})")
adata_raw = ad.concat(parts, merge="same")  # one AnnData with all samples (anndata concat, as integration md 30 says)  # CHANGED: inner join = shared genes
if not adata_raw.obs_names.equals(markers_extra.obs_names):  # CHANGED: the two objects must hold the same cells in the same order
    raise SystemExit("marker genes outside the shared set: cell order differs from the atlas")
del parts
adata_raw.layers["counts"] = adata_raw.X.copy()  # the corrected counts
# ===== end of the part not in the book =====

# --- integration cell 17 ---
# CHANGED: the dataset is the merged atlas above, not the lamindb artifact
adata_raw.layers["logcounts"] = adata_raw.X
# --- integration cell 19 ---
# CHANGED: no label_key: cell types are assigned after clustering
batch_key = "library"  # CHANGED: the 10x library (capture) as batch (the book uses the sample, the same unit where a library holds one donor)
# --- integration cell 21 (value_counts display) left out ---
# --- integration cell 23 ---
# CHANGED: all batches kept (the book keeps three to save time)
adata = adata_raw
del adata_raw
# --- integration cell 25 (feature_types display) left out ---
# --- integration cell 27 ---
# CHANGED: gene expression only already; no feature_types subset
sc.pp.filter_genes(adata, min_cells=20)  # CHANGED: min_cells=20 (on the merged atlas; the book uses 1 here)
# CHANGED: cells with no counts in the atlas genes are dropped (no size factor; all their counts are in
# excluded or unshared genes), with the same cells of the markers kept apart; scanpy's filter, in place
empty = np.asarray(adata.X.sum(axis=1)).ravel() == 0
log("cells with no counts in the atlas genes, dropped: " + (", ".join(
    f"{k} {v}" for k, v in adata.obs.dataset[empty].value_counts().items() if v) if empty.any() else "none"))
sc.pp.filter_cells(adata, min_counts=1)
markers_extra = markers_extra[~empty].copy()
if not adata.obs_names.equals(markers_extra.obs_names):
    raise SystemExit("markers kept apart: cell order differs from the atlas after dropping empty cells")
# CHANGED: the same filter for the markers kept apart; a dataset that did not measure a gene has no
# entries for it, so only cells of the datasets that measured it are counted
extra_n = pd.Series(np.asarray((markers_extra.X > 0).sum(axis=0)).ravel(), index=markers_extra.var_names)
markers_extra = markers_extra[:, (extra_n >= 20).values].copy()
# CHANGED: why each curated marker is or is not available for annotation (03 reads this table)
cur = sorted({m for v in curated_markers()[1].values() for m in v})
excl = set(genes.gene_name[(group.reindex(genes.gene_id).fillna("") != "").values])
status = pd.Series(
    ["in_atlas" if m in adata.var_names else
     "kept_apart" if m in markers_extra.var_names else
     "excluded_decision_2" if m in excl else
     "under_20_cells" if m in measured_any else "not_measured" for m in cur], index=cur, name="status")
status.rename_axis("gene").to_csv(out / "curated_markers_status.tsv", sep="\t")
log("curated markers: " + "; ".join(f"{k} {v}" for k, v in status.value_counts().items()) +
    "; under 20 cells: " + ", ".join(status.index[status == "under_20_cells"]) +
    "; excluded: " + ", ".join(status.index[status == "excluded_decision_2"]) +
    "; not measured: " + ", ".join(status.index[status == "not_measured"]))
adata.var.to_csv(out / "genes.tsv", sep="\t")  # CHANGED: record of the genes kept
# --- integration cell 29 ---
adata.X = adata.layers["counts"].copy()
lib_size = np.asarray(adata.X.sum(axis=1)).ravel()  # CHANGED: each cell's counts over the atlas genes, for the markers kept apart
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
# CHANGED: the markers kept apart are normalised with each cell's atlas size factor (normalize_total's default
# target, the median library size), so their values are on the scale of the atlas genes; then log1p
scale = np.median(lib_size[lib_size > 0]) / lib_size
markers_extra.layers["counts"] = markers_extra.X.copy()
markers_extra.X = sp.csr_matrix(markers_extra.X.multiply(scale[:, None]))
markers_extra.X.data = np.log1p(markers_extra.X.data)
measured = markers_extra.var.filter(like="measured_").fillna(False).astype(bool)
markers_extra.var = measured.assign(n_datasets=measured.sum(1).values)
markers_extra.uns["datasets"] = sorted(c.removeprefix("measured_") for c in measured.columns)
markers_extra.write_h5ad(out / "markers_unshared.h5ad")
log(f"markers outside the shared genes: {markers_extra.n_vars} written with the datasets that measured each")
# --- integration cell 32 ---
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
sc.pp.neighbors(adata, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
sc.tl.umap(adata, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"]  # CHANGED: kept, the scVI UMAP below overwrites X_umap
# --- integration cell 34 (UMAP plot) left out ---
# --- integration cell 38 ---
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key
)
adata.var.to_csv(out / "genes.tsv", sep="\t")  # CHANGED: display of adata.var replaced by a file
# --- integration cell 40 (n_batches plot) left out ---
# --- integration cell 42 ---
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
# --- integration cell 45 ---
adata_scvi = adata_hvg.copy()
# --- integration cell 48 ---
scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_key)
# --- integration cell 52 ---
model_scvi = scvi.model.SCVI(adata_scvi)
log(model_scvi)  # CHANGED: display -> log
# --- integration cell 54 (view_anndata_setup display) left out ---
# --- integration cell 58 ---
max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
log(f"max_epochs_scvi {max_epochs_scvi}")  # CHANGED: print -> log
# --- integration cell 60 ---
# CHANGED: the book's model_scvi.train() uses scvi-tools' epoch rule, which gives 6 epochs for 1.36 M
# cells while the default KL warm-up spans 400 epochs (KL weight 0.013 at the end). Instead: per-epoch validation on
# the default 10% split, early stopping on validation ELBO (default patience 45), KL warm-up over 20 epochs and no
# stopping before epoch 40, so the KL weight is 1 well before training stops; max_epochs 400 (the rule's own cap).
model_scvi.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_warmup_epochs=40,
    check_val_every_n_epoch=1,
    plan_kwargs={"n_epochs_kl_warmup": 20},
)
hist = pd.concat(model_scvi.history.values(), axis=1)  # CHANGED: loss curves kept
hist.index.name = "epoch"
hist.to_csv(out / "scvi_history.tsv", sep="\t")
log(f"scVI trained {len(hist)} epochs; final kl_weight {hist.kl_weight.iloc[-1]:.3f}; best elbo_validation "
    f"{hist.elbo_validation.min():.3f} at epoch {int(hist.elbo_validation.astype(float).idxmin()) + 1}")
model_scvi.save(out / "scvi_model", overwrite=True)  # CHANGED: model kept
# --- integration cell 64 ---
adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
# --- integration cell 67 ---
sc.pp.neighbors(adata_scvi, use_rep="X_scVI", random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
sc.tl.umap(adata_scvi, random_state=0)  # CHANGED: scanpy default random_state=0 written out (seed)
# --- integration cell 69 (UMAP plot) left out; scANVI, BBKNN, Seurat and benchmarking cells not used ---

# --- clustering cell 12 ---
sc.tl.leiden(adata_scvi, flavor="igraph", n_iterations=2, random_state=0)  # CHANGED: adata_scvi; default random_state=0 written out (seed)
# --- clustering cell 15 ---
sc.tl.leiden(
    adata_scvi, key_added="leiden_res0_25", resolution=0.25, flavor="igraph", n_iterations=2, random_state=0  # CHANGED: adata_scvi; seed
)
sc.tl.leiden(
    adata_scvi, key_added="leiden_res0_5", resolution=0.5, flavor="igraph", n_iterations=2, random_state=0  # CHANGED: adata_scvi; seed
)
sc.tl.leiden(
    adata_scvi, key_added="leiden_res1", resolution=1.0, flavor="igraph", n_iterations=2, random_state=0  # CHANGED: adata_scvi; seed
)
# --- clustering cell 17 (UMAP plot) left out ---

# CHANGED: the integrated embedding, UMAP and clusterings are carried to the full-gene object, which is written out
adata.obsm["X_scVI"] = adata_scvi.obsm["X_scVI"]
adata.obsm["X_umap"] = adata_scvi.obsm["X_umap"]
adata.obsp["connectivities"] = adata_scvi.obsp["connectivities"]
adata.obsp["distances"] = adata_scvi.obsp["distances"]
adata.uns["neighbors"] = adata_scvi.uns["neighbors"]
for k in ["leiden", "leiden_res0_25", "leiden_res0_5", "leiden_res1"]:
    adata.obs[k] = adata_scvi.obs[k]
adata.write_h5ad(out / "atlas.h5ad")
log(f"atlas: {adata.n_obs} cells, {adata.n_vars} genes, {adata.var.highly_variable.sum()} HVGs; "
    f"Leiden clusters: {adata.obs.leiden.nunique()} (res 1), {adata.obs.leiden_res0_5.nunique()} (0.5), {adata.obs.leiden_res0_25.nunique()} (0.25)")
