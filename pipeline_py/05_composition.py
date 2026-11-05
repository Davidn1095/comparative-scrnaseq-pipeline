# Compositional analysis of cell-type counts per donor with scCODA.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19):
#   jupyter-book/conditions/compositional.ipynb, setup cells 7-8 and cells 25-40 ("With labeled clusters", scCODA).
#   tascCODA and Milo are not ported.
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason. Plots
# are written to files instead of displayed.
# One model with all donors (2026-11-02): formula batch + sex + condition (batch = dataset, or Perez's processing
# cohort), condition levels control, SLE and SjD,
# reference cell type cDC2; fitted with base level control (SLE vs control, SjD vs control) and with base level SLE
# (SjD vs SLE, the adjusted indirect comparison).
# Usage: python 05_composition.py   (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration;
#   SCCODA_BASES=control or SLE runs one of the two fits, default both)
# Input: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/atlas.h5ad (obs only: donor_id, sample_id, dataset, condition) and
#        03_annotation/cells_manual_annotation.tsv.gz (manual cell-type labels)
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/05_composition/: reference_cell_type.txt; figures/ (boxplots, stacked
#   barplot); <contrast>/ for SLE_vs_control, SjD_vs_control and SjD_vs_SLE: effects.tsv (at FDR 0.2) and
#   credible_effects.tsv (all covariates of its fit); base_<level>/figures/effects_barplot.png

# --- cell 7 ---
import logging
import warnings

logging.getLogger("matplotlib").setLevel(logging.WARNING)
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# --- cell 8 ---
import os  # CHANGED: run folder

import anndata as ad  # CHANGED: the cell table is read into an AnnData without expression
import h5py  # CHANGED: the cell table is read alone from atlas.h5ad
import matplotlib

matplotlib.use("Agg")  # CHANGED: plots are written to files (%matplotlib inline in the notebook)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pertpy as pt

from config import OUT_ROOT, batch_of, excluded_labels, log  # CHANGED: pipeline paths; labels left out; batch

# CHANGED: lamindb, mudata, scanpy, schist, scvi and seaborn not imported: used only by the notebook's data loading,
# tascCODA and Milo sections, which are not ported; ln.connect and ln.track removed

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline input
out_root = run / "05_composition"  # CHANGED: pipeline output


def save(fig_dir, name, fig=None):  # CHANGED: plots are written to files instead of displayed
    fig_dir.mkdir(parents=True, exist_ok=True)
    (fig or plt.gcf()).savefig(fig_dir / name, dpi=150, bbox_inches="tight")
    plt.close("all")


# --- cell 9 (data loading) replaced ---
# CHANGED: cell metadata of the atlas of 02_integration.py with the manual cell-type labels of 03_annotation.py; cells
# without a label are not used. scCODA uses only the cell metadata, so no expression is loaded.
with h5py.File(run / "atlas.h5ad", "r") as f:  # CHANGED: only the cell table is read (backed mode still loads the layers)
    obs = ad.io.read_elem(f["obs"])[["donor_id", "sample_id", "dataset", "condition", "inferred_sex"]].copy()  # CHANGED: as above
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)  # CHANGED: as above
obs["cell_type"] = labels["manual_celltype_annotation"].reindex(obs.index)  # CHANGED: as above
obs["batch"] = batch_of(obs).values  # CHANGED: batch = dataset, or Perez's processing cohort (Perez_1-4): Perez's cohorts are unbalanced between SLE and controls
obs["sex"] = obs.pop("inferred_sex").astype(str)  # CHANGED: inferred donor sex, a covariate as in the DE design (2026-11-03)
obs = obs[obs["cell_type"].notna() & ~excluded_labels(obs)]  # CHANGED: as above; cells labelled low quality or platelets are left out
for c in obs.columns:  # CHANGED: plain strings for the formula and the sample table
    obs[c] = obs[c].astype(str)

# CHANGED: one model with all donors (2026-11-02; replaces two separate models, one per disease):
# formula batch + sex + condition (batch from 2026-11-02, sex from 2026-11-03), condition levels control (HC), SLE and SjD, reference cell type cDC2. Fitted twice with
# the same data, prior and rng_key: with base level control it gives SLE vs HC and SjD vs HC; with base level SLE it gives
# SjD vs SLE, the adjusted indirect comparison (SLE and SjD donors come from different datasets; every dataset has its
# own controls, so the dataset effects are estimated from them).
REFERENCE = "cDC2"  # CHANGED: fixed reference cell type
fits = {"control": ["SLE", "SjD"], "SLE": ["SjD"]}  # CHANGED: base level -> condition levels reported against it
fits = {k: v for k, v in fits.items() if k in os.environ.get("SCCODA_BASES", "control,SLE").split(",")}  # CHANGED: SCCODA_BASES=control or SLE runs one fit per job (default both); each fit writes only its own contrast folders
adata = ad.AnnData(obs=obs)  # CHANGED: all donors
log(f"one model: {adata.n_obs} cells, {obs['donor_id'].nunique()} donors, {obs['cell_type'].nunique()} cell types")  # CHANGED: log

# --- cell 25 ---
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_type",  # CHANGED: this atlas's column (cell_label in the notebook)
    sample_identifier="donor_id",  # CHANGED: one sample per donor (batch in the notebook)
    covariate_obs=["condition", "batch", "sex"],  # CHANGED: batch and sex added for the formula below
)
fig_all = out_root / "figures"  # CHANGED: pipeline output

# --- cell 27 ---
sccoda_model.plot_boxplots(
    sccoda_data,
    modality_key="coda",
    feature_name="condition",
    figsize=(12, 5),
    add_dots=True,
)
save(fig_all, "boxplots.png")  # CHANGED: saved to file (plt.show() in the notebook)

# --- cell 29 ---
sccoda_model.plot_stacked_barplot(
    sccoda_data, modality_key="coda", feature_name="condition", figsize=(4, 2)
)
save(fig_all, "stacked_barplot.png")  # CHANGED: saved to file (plt.show() in the notebook)

automatic = sccoda_model.prepare(sccoda_data.copy(), modality_key="coda", formula="batch + sex + condition",
                                 reference_cell_type="automatic")["coda"].uns["scCODA_params"]["reference_cell_type"]  # CHANGED: pertpy's automatic choice, recorded only
(out_root / "reference_cell_type.txt").write_text(f"{REFERENCE}\nrule: fixed; pertpy's automatic choice for this model: {automatic}\n")  # CHANGED: as above
log(f"reference cell type {REFERENCE} (pertpy's automatic choice: {automatic})")  # CHANGED: log

for base, levels in fits.items():  # CHANGED: one fit per base level (see above)
    formula = f"batch + sex + C(condition, Treatment('{base}'))"  # CHANGED: batch + sex + condition, coded against the base level
    data = sccoda_data.copy()  # CHANGED: each fit starts from the loaded data

    # --- cell 31 ---
    data = sccoda_model.prepare(
        data,
        modality_key="coda",
        formula=formula,  # CHANGED: as above
        reference_cell_type=REFERENCE,  # CHANGED: as above ("Endocrine", a cell type of its gut dataset, in the notebook)
    )
    sccoda_model.run_nuts(data, modality_key="coda", rng_key=42)

    # --- cell 33 ---
    # CHANGED: the effects are written after set_fdr (cell 36), not here: before it the table holds the effects at
    # pertpy's default FDR of 0.05, not at the FDR 0.2 of the credible effects and the barplot

    # --- cell 34 (display of sccoda_data) left out ---

    # --- cell 36 ---
    sccoda_model.set_fdr(data, 0.2)
    params = data["coda"].uns["scCODA_params"]  # CHANGED: pertpy misaligns its SD column (values in cell-type-major
    beta = np.asarray(params["mcmc"]["samples"]["beta"])  # order, rows covariate-major), so the SD of each effect is
    # taken from its MCMC samples (draws x covariates x cell types)
    credible = sccoda_model.credible_effects(data, modality_key="coda")  # CHANGED: cell 38, computed once per fit
    for level in levels:  # CHANGED: one output folder per contrast
        cov = f"C(condition, Treatment('{base}'))[T.{level}]"
        effect_key = f"effect_df_{cov}"
        contrast = f"{level}_vs_{base}"  # SLE_vs_control, SjD_vs_control, SjD_vs_SLE
        out = out_root / contrast
        out.mkdir(parents=True, exist_ok=True)
        effects = data["coda"].varm[effect_key].copy()  # CHANGED: cell 33's display, written out after set_fdr (see above)
        effects["SD"] = beta[:, params["covariate_names"].index(cov), :].std(axis=0, ddof=1)  # CHANGED: as above
        effects.to_csv(out / "effects.tsv", sep="\t")  # CHANGED: as above

        # --- cell 38 ---
        credible.to_frame("credible").to_csv(out / "credible_effects.tsv", sep="\t")  # CHANGED: written out instead of displayed; all covariates of this fit
        n = int(credible.xs(cov, level="Covariate").sum())
        log(f"{contrast}: formula {formula}; reference cell type {REFERENCE}; threshold_prob "
            f"{params['threshold_prob']:.3f}; credible effects at FDR 0.2: {n}")  # CHANGED: log

    # --- cell 40 ---
    sccoda_model.plot_effects_barplot(
        data, modality_key="coda", covariates="condition"
    )
    save(out_root / f"base_{base}" / "figures", "effects_barplot.png")  # CHANGED: saved to file (plt.show() in the notebook)
log("done")  # CHANGED: log
