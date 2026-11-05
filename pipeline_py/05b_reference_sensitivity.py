# Reference sensitivity of the scCODA composition model: step 05's model refitted with each of six
# reference cell types, one fit per reference and base level; then a summary of the credible effects across references.
# Ported from sc-best-practices (book commit dd6eb86, 2026-10-19), as 05_composition.py:
#   jupyter-book/conditions/compositional.ipynb, setup cells 7-8 and cells 25, 31, 36 and 38 (scCODA load, prepare and
#   NUTS, FDR, credible effects); the plots of cells 27-29 and 40 are left out (step 05 draws them).
# The notebook's code cells are copied in order; every changed line carries a "# CHANGED:" comment with the reason.
# Model: as step 05, one model with all donors, formula batch + sex + C(condition, Treatment(<base>)), rng_key 42, FDR 0.2.
# References: cDC2 (step 05's) and the five next cell types whose proportions vary least across donors, Regulatory T,
# Non-classical monocytes, NK, CD8 effector memory T and pDC. Sampler diagnostics per fit: divergences (numpyro's
# "diverging" field, collected on top of pertpy's), acceptance, tree-depth saturation, step size, and per parameter the
# split R-hat and effective sample size of the single chain.
# Usage: python 05b_reference_sensitivity.py <base> <reference>   one fit; base control or SLE, the reference with "_" for
#        spaces (e.g. Regulatory_T); one Slurm job per fit (QOS short)
#        python 05b_reference_sensitivity.py --summary              after the 12 fits
#        (ATLAS_RUN=<folder of 02_integration.py>, default 02_integration)
# Input: as step 05 (atlas.h5ad obs, 03_annotation/cells_manual_annotation.tsv.gz)
# Output: $PIPELINE_OUT/py/atlas/<ATLAS_RUN>/05b_reference_sensitivity/
#   base_<base>__<reference>.tsv   effects of every condition covariate of the fit (pertpy's table, SD from the samples,
#                                  credible at FDR 0.2)
#   diag_base_<base>__<reference>.tsv  sampler diagnostics of the fit
#   reference_sensitivity.tsv      (--summary) per contrast and cell type: references with a credible effect / references
#                                  tested (the cell type itself excluded), sign of the credible effects, the cDC2 fit's
#                                  log2 fold change and the range over references, and all_references (yes/no);
#                                  SLE vs SjD is the base-SLE fit with its sign reversed, as in step 05 and Fig 2e
#   diagnostics.tsv                (--summary) the 12 fits' diagnostics

# --- cell 7 ---
import logging
import warnings

logging.getLogger("matplotlib").setLevel(logging.WARNING)
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# --- cell 8 ---
import os  # CHANGED: run folder
import sys  # CHANGED: base level and reference from the command line

import anndata as ad  # CHANGED: the cell table is read into an AnnData without expression
import h5py  # CHANGED: the cell table is read alone from atlas.h5ad
import numpy as np
import pandas as pd
import pertpy as pt
import numpyro.infer.mcmc as numpyro_mcmc  # CHANGED: divergences are collected (see below)
from numpyro.diagnostics import effective_sample_size, split_gelman_rubin  # CHANGED: sampler diagnostics

from config import OUT_ROOT, batch_of, excluded_labels, log  # CHANGED: pipeline paths; labels left out; batch

# CHANGED: lamindb, mudata, scanpy, schist, scvi and seaborn not imported (see step 05); ln.connect and ln.track removed

run = OUT_ROOT / "atlas" / os.environ.get("ATLAS_RUN", "02_integration")  # CHANGED: pipeline input
out = run / "05b_reference_sensitivity"  # CHANGED: pipeline output
out.mkdir(parents=True, exist_ok=True)  # CHANGED: as above
REFERENCES = ["cDC2", "Regulatory T", "Non-classical monocytes", "NK", "CD8 effector memory T", "pDC"]  # CHANGED: see header
CONTRASTS = {"SLE vs HC": ("control", "SLE", 1), "SjD vs HC": ("control", "SjD", 1), "SLE vs SjD": ("SLE", "SjD", -1)}  # CHANGED: as step 05


def tag(base, ref):  # CHANGED: file names
    return f"base_{base}__{ref.replace(' ', '_')}"


if sys.argv[1:] == ["--summary"]:  # CHANGED: summary of the 12 fits (no notebook counterpart)
    fits = {(b, r): out / f"{tag(b, r)}.tsv" for b in ("control", "SLE") for r in REFERENCES}
    missing = [str(p.name) for p in fits.values() if not p.exists()]
    if missing:
        sys.exit(f"missing fits: {missing}")
    r = pd.concat(pd.read_csv(p, sep="\t") for p in fits.values())
    rows = []
    for k, (base, lev, s) in CONTRASTS.items():
        x = r[(r.base == base) & (r.covariate == f"C(condition, Treatment('{base}'))[T.{lev}]")].copy()
        x["lfc"] = s * x["log2-fold change"]
        for ct in x.cell_type.unique():
            y = x[(x.cell_type == ct) & (x.reference != ct)]
            cr = y[y.credible.astype(str) == "True"]
            sign = "+" if len(cr) and (cr.lfc > 0).all() else "-" if len(cr) and (cr.lfc < 0).all() else ("mixed" if len(cr) else "")
            l0 = x[(x.cell_type == ct) & (x.reference == "cDC2")].lfc
            rows.append(dict(contrast=k, cell_type=ct, credible_with=f"{len(cr)}/{len(y)}", n_credible=len(cr), n_references=len(y),
                             sign_credible=sign, log2FC_cDC2=round(float(l0.iloc[0]), 4) if len(l0) else np.nan,
                             log2FC_min=round(float(y.lfc.min()), 4), log2FC_max=round(float(y.lfc.max()), 4),
                             all_references="yes" if len(y) and len(cr) == len(y) else "no"))
    t = pd.DataFrame(rows)
    t.to_csv(out / "reference_sensitivity.tsv", sep="\t", index=False)
    d = pd.concat([pd.read_csv(out / f"diag_{tag(b, ref)}.tsv", sep="\t", index_col=0)["value"].rename(f"{b}/{ref}")
                   for b, ref in fits], axis=1).T if all((out / f"diag_{tag(b, ref)}.tsv").exists() for b, ref in fits) else None  # one row per fit
    if d is not None:
        d.to_csv(out / "diagnostics.tsv", sep="\t")
    log("effects credible with every reference: " + ", ".join(f"{k} {int((t[t.contrast == k].all_references == 'yes').sum())}" for k in CONTRASTS))
    sys.exit(0)

base, reference = sys.argv[1], sys.argv[2].replace("_", " ")  # CHANGED: one fit per job
assert base in ("control", "SLE") and reference in REFERENCES, (base, reference)

# CHANGED: numpyro's "diverging" field is collected with pertpy's extra fields (pertpy tracks potential energy, number
# of steps, step size and acceptance only); collecting a field does not change the sampling
_run = numpyro_mcmc.MCMC.run


def _run_with_divergences(self, *args, extra_fields=(), **kwargs):
    return _run(self, *args, extra_fields=tuple(extra_fields) + ("diverging",), **kwargs)


numpyro_mcmc.MCMC.run = _run_with_divergences

# --- cell 9 (data loading) replaced ---
# CHANGED: as step 05: the cell table of atlas.h5ad with the manual labels; cells without a label, labelled low quality
# or platelets are not used
with h5py.File(run / "atlas.h5ad", "r") as f:  # CHANGED: only the cell table is read
    obs = ad.io.read_elem(f["obs"])[["donor_id", "sample_id", "dataset", "condition", "inferred_sex"]].copy()  # CHANGED: as above
labels = pd.read_csv(run / "03_annotation" / "cells_manual_annotation.tsv.gz", sep="\t", index_col=0)  # CHANGED: as above
obs["cell_type"] = labels["manual_celltype_annotation"].reindex(obs.index)  # CHANGED: as above
obs["batch"] = batch_of(obs).values  # CHANGED: batch as step 05
obs["sex"] = obs.pop("inferred_sex").astype(str)  # CHANGED: inferred donor sex as step 05
obs = obs[obs["cell_type"].notna() & ~excluded_labels(obs)]  # CHANGED: as above
for c in obs.columns:  # CHANGED: plain strings for the formula and the sample table
    obs[c] = obs[c].astype(str)
adata = ad.AnnData(obs=obs)  # CHANGED: as step 05

# --- cell 25 ---
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_type",  # CHANGED: this atlas's column (cell_label in the notebook)
    sample_identifier="donor_id",  # CHANGED: one sample per donor (batch in the notebook)
    covariate_obs=["condition", "batch", "sex"],  # CHANGED: as step 05
)

# --- cell 31 ---
formula = f"batch + sex + C(condition, Treatment('{base}'))"  # CHANGED: step 05's formula, coded against the base level
sccoda_data = sccoda_model.prepare(
    sccoda_data,
    modality_key="coda",
    formula=formula,  # CHANGED: as above
    reference_cell_type=reference,  # CHANGED: the reference of this fit ("Endocrine" in the notebook)
)
sccoda_model.run_nuts(sccoda_data, modality_key="coda", rng_key=42)

# --- cell 36 ---
sccoda_model.set_fdr(sccoda_data, 0.2)

# --- cell 38 ---
credible = sccoda_model.credible_effects(sccoda_data, modality_key="coda")  # CHANGED: kept for the output below
s = sccoda_data["coda"]  # CHANGED: effects of every condition covariate written out, as step 05
params = s.uns["scCODA_params"]  # CHANGED: as above
beta = np.asarray(params["mcmc"]["samples"]["beta"])  # CHANGED: pertpy misaligns its SD column (step 05)
rows = []  # CHANGED: as above
for key in [k for k in s.varm if k.startswith("effect_df_") and "condition" in k]:  # CHANGED: as above
    cov = key.removeprefix("effect_df_")
    eff = s.varm[key].copy()
    eff["SD"] = beta[:, params["covariate_names"].index(cov), :].std(axis=0, ddof=1)
    eff["credible"] = credible.xs(cov, level="Covariate").reindex(eff.index).values
    eff.insert(0, "covariate", cov)
    rows.append(eff)
res = pd.concat(rows)  # CHANGED: as above
res.insert(0, "reference", reference)  # CHANGED: as above
res.insert(0, "base", base)  # CHANGED: as above
res.to_csv(out / f"{tag(base, reference)}.tsv", sep="\t", index_label="cell_type")  # CHANGED: as above

# CHANGED: sampler diagnostics (no notebook counterpart): one chain, so R-hat is the split R-hat of that chain
ef = sccoda_model.mcmc.get_extra_fields()
samples = sccoda_model.mcmc.get_samples(group_by_chain=True)
steps = np.asarray(ef["num_steps"])
diag = {"base": base, "reference": reference, "num_chains": sccoda_model.mcmc.num_chains, "num_samples": steps.size,
        "divergences": int(np.asarray(ef["diverging"]).sum()), "mean_accept_prob": float(np.asarray(ef["accept_prob"]).mean()),
        "step_size": float(np.asarray(ef["adapt_state.step_size"]).mean()),
        "tree_depth_saturated_frac": float((steps >= 2 ** 10 - 1).mean()),  # numpyro's default maximum tree depth, 10
        "median_num_steps": float(np.median(steps))}
for site, x in samples.items():
    x = np.asarray(x).reshape(np.asarray(x).shape[0], np.asarray(x).shape[1], -1)
    keep = x.std(axis=(0, 1)) > 0
    diag[f"max_rhat_{site}"] = float(np.nanmax(split_gelman_rubin(x[:, :, keep]))) if keep.any() else np.nan
    diag[f"min_ess_{site}"] = float(np.nanmin(effective_sample_size(x[:, :, keep]))) if keep.any() else np.nan
pd.Series(diag).to_frame("value").to_csv(out / f"diag_{tag(base, reference)}.tsv", sep="\t")
log(f"{formula}; reference {reference}; threshold_prob {params['threshold_prob']:.3f}; divergences {diag['divergences']}; "
    f"credible condition effects {int(res.credible.astype(bool).sum())}")
