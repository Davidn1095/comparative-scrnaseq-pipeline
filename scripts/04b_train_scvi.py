#!/usr/bin/env python3
"""
Train scVI for fig 1 visualisation:
 - Load full atlas counts (~869k cells, 25 common cell types, 2000 HVGs).
 - Stratified subsample (N=SCVI_N_TRAIN, default 10k) within dataset_id x
   condition strata for TRAINING ONLY.
 - Train scVI (n_latent=30, n_layers=2, n_hidden=128, 200 epochs early stop).
 - INFERENCE on the FULL atlas via get_latent_representation(adata).
 - UMAP on the full-atlas latent.

Outputs (results/integration/scvi_visualisation/):
  scvi_model/             Trained model
  training_history.json   Loss curves + epoch count
  latent.npy              (n_full_cells x 30) latent for ALL cells
  scvi_umap.tsv           Per-cell UMAP coords + dataset/condition/cell_type
"""
import json
import os
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scvi
import torch
from anndata import AnnData

# Repository root: ATLAS_ROOT if set, otherwise the directory above scripts/.
ATLAS_ROOT = Path(os.environ.get("ATLAS_ROOT", Path(__file__).resolve().parents[1]))
EXTRACTED = ATLAS_ROOT / "results/integration/scvi_visualisation/extracted"
OUT_NAME = os.environ.get("SCVI_OUT_DIR", "scvi_visualisation")
OUT = ATLAS_ROOT / "results/integration" / OUT_NAME
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR = OUT / "scvi_model"

# SCVI_N_TRAIN: 0 (or unset/<=0) means train on ALL cells; positive N means
# stratified subsample of N cells for training (inference still on full).
N_TRAIN = int(os.environ.get("SCVI_N_TRAIN", "0"))

# Methods §4.3 states seeds were set to 42 for all stochastic operations; scVI
# training was previously unseeded, so the latent was not reproducible.
scvi.settings.seed = 42

print(f"CUDA available: {torch.cuda.is_available()}", flush=True)
print(f"scvi-tools: {scvi.__version__}", flush=True)
print(f"scvi seed:  {scvi.settings.seed}", flush=True)
print(f"scanpy:     {sc.__version__}", flush=True)
print(f"Train subsample size: {N_TRAIN}", flush=True)

# ---- Load full counts ----
print("Loading counts MTX ...", flush=True)
t0 = time.time()
counts = scipy.io.mmread(EXTRACTED / "counts.mtx").tocsr()  # genes x cells
counts = counts.T.tocsr()                                    # cells x genes
barcodes = pd.read_csv(EXTRACTED / "barcodes.tsv", header=None)[0].tolist()
features = pd.read_csv(EXTRACTED / "features.tsv", header=None)[0].tolist()
meta = pd.read_csv(EXTRACTED / "metadata.tsv", sep="\t", dtype=str)
meta = meta.set_index("barcode").reindex(barcodes)
print(f"  Loaded counts {counts.shape} + {len(barcodes)} cells in {time.time()-t0:.1f}s", flush=True)

adata = AnnData(X=counts, obs=meta, var=pd.DataFrame(index=features))
adata.layers["counts"] = adata.X.copy()
adata.obs["dataset_id"] = adata.obs["dataset_id"].astype("category")
adata.obs["condition"]  = adata.obs["condition"].astype("category")
print(adata, flush=True)

# ---- Training data: full atlas (N_TRAIN<=0) or stratified subsample ----
if N_TRAIN <= 0:
    print(f"Training on FULL atlas ({adata.n_obs} cells)", flush=True)
    adata_train = adata
    n_train_cells = adata.n_obs
else:
    rng = np.random.RandomState(42)
    strata = adata.obs.groupby(["dataset_id", "condition"], observed=True).indices
    total_cells = sum(len(v) for v in strata.values())
    train_idx = []
    print("Stratified subsample allocation (proportional within dataset x condition):", flush=True)
    for key, idx in sorted(strata.items()):
        n_k = max(1, int(round(len(idx) * N_TRAIN / total_cells)))
        n_k = min(n_k, len(idx))
        pick = rng.choice(idx, size=n_k, replace=False)
        train_idx.extend(pick.tolist())
        print(f"  {key} : stratum n={len(idx):>7} -> picked {n_k}", flush=True)
    train_idx = np.array(train_idx)
    print(f"Total training cells: {len(train_idx)}", flush=True)
    adata_train = adata[train_idx, :].copy()
    n_train_cells = len(train_idx)
print(f"Training adata: {adata_train.shape}", flush=True)
scvi.model.SCVI.setup_anndata(adata_train, layer="counts", batch_key="dataset_id")

model = scvi.model.SCVI(adata_train, n_layers=2, n_latent=30, n_hidden=128)
print("Model created, starting training ...", flush=True)
t0 = time.time()
model.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=20,
    early_stopping_min_delta=0.001,
    train_size=0.9,
    accelerator="auto",
    devices="auto",
)
train_minutes = (time.time() - t0) / 60.0
n_epochs = len(model.history["elbo_train"])
print(f"Training done in {train_minutes:.1f} min, {n_epochs} epochs", flush=True)

# ---- Save model + history ----
model.save(str(MODEL_DIR), overwrite=True)
print(f"Saved model -> {MODEL_DIR}", flush=True)

history = {k: v.iloc[:, 0].tolist() for k, v in model.history.items() if hasattr(v, "iloc")}
with open(OUT / "training_history.json", "w") as f:
    json.dump({
        "train_minutes":   train_minutes,
        "n_epochs":        n_epochs,
        "n_train_cells":   int(n_train_cells),
        "n_total_cells":   int(adata.n_obs),
        "final_elbo_train": float(model.history["elbo_train"].iloc[-1, 0]),
        "final_elbo_val":   float(model.history["elbo_validation"].iloc[-1, 0]),
        "history":         history,
    }, f, indent=2)
print("Saved training_history.json", flush=True)

# ---- Inference on FULL atlas ----
print(f"Setting up full adata for inference ({adata.n_obs} cells)...", flush=True)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset_id")
t0 = time.time()
latent_full = model.get_latent_representation(adata)
print(f"Inference done in {(time.time()-t0)/60:.1f} min, latent shape {latent_full.shape}", flush=True)
np.save(OUT / "latent.npy", latent_full)
print(f"Saved latent -> {OUT/'latent.npy'}", flush=True)

# ---- UMAP on full-atlas latent ----
print("Computing UMAP on full-atlas latent ...", flush=True)
adata.obsm["X_scvi"] = latent_full
t0 = time.time()
sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3)
print(f"UMAP done in {(time.time()-t0)/60:.1f} min", flush=True)

umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"],
                       index=adata.obs_names)
umap_df = umap_df.join(adata.obs[["dataset_id", "condition", "cell_type"]])
umap_df.to_csv(OUT / "scvi_umap.tsv", sep="\t")
print(f"Saved UMAP -> {OUT/'scvi_umap.tsv'}", flush=True)
print("=== Done ===", flush=True)
