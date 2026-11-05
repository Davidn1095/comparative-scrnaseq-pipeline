"""Paths and dataset sheet for the pipeline (pipeline_py/), which ports the sc-best-practices notebooks
(theislab/single-cell-best-practices, book commit dd6eb86, 2026-10-19). Per-library inputs are written by
00_deposited_input.py to $PIPELINE_OUT/input/<dataset>/ (layout in its docstring)."""
import os
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parent
REPO = Path(os.environ.get("ATLAS_ROOT", PIPELINE_DIR.parent))
DEV_ROOT = Path(os.environ.get("PIPELINE_OUT", REPO / "out"))  # outputs
INPUT_ROOT = DEV_ROOT / "input"
OUT_ROOT = DEV_ROOT / "py"
REF_DIR = Path(os.environ.get("REF_DIR", REPO / "refdata" / "refdata-gex-GRCh38-2024-A"))  # Cell Ranger reference
DATA_ROOT = Path(os.environ.get("DATA_ROOT", REPO / "data"))  # raw downloads, data/<disease>/<accession>/ (README)
DATASETS = pd.read_csv(PIPELINE_DIR / "datasets.tsv", sep="\t", comment="#", dtype=str)

# scDblFinder's default expected doublet rate per 1,000 cells (scDblFinder 1.26.9,
# dbr.per1k = 0.008), always passed explicitly
DBL_RATE_PER_1K = 0.008


def dataset_dir(ds, step):
    return OUT_ROOT / "datasets" / ds / step


def log(*a):
    import datetime
    print(datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]"), *a, flush=True)


# manual labels of clusters left out of steps 04-06 (06 reads the DE tables of 04)
EXCLUDED_LABELS = ("low quality", "platelets")


# markers of the supplementary atlas dot plot (09_supplementary_data.py, 10_supplementary.R), per final cell
# type in the figures' order: positive markers named in 03_annotation/manual_annotation.tsv. Illustrative only, not used
# for annotation (step 03 uses the curated markers below)
SUPP_DOTPLOT_MARKERS = {
    "Naive CD4 T": ["CCR7", "LEF1", "TCF7"],
    "Memory CD4 T": ["IL7R", "AQP3", "KLRB1"],
    "Regulatory T": ["FOXP3", "IL2RA", "RTKN2"],
    "Naive CD8 T": ["CD8A", "CD8B", "CCR7"],
    "CD8 effector memory T": ["GZMK", "CCL5"],
    "CD8 effector T (TEMRA)": ["GZMH", "FGFBP2"],
    "MAIT and gamma-delta T": ["SLC4A10", "ZBTB16", "TRDC"],
    "NK": ["GNLY", "KLRF1"],
    "Cycling T&NK": ["MKI67", "TOP2A"],
    "Naive B": ["MS4A1", "TCL1A", "IGHD"],
    "Memory B": ["TNFRSF13B", "CD27"],
    "Classical monocytes": ["CD14", "S100A12", "VCAN"],
    "Non-classical monocytes": ["FCGR3A", "CX3CR1"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A"],
    "pDC": ["LILRA4", "IL3RA"],
    "HSPC": ["CD34", "SPINK2"],
}


def excluded_labels(obs, label="cell_type"):
    """Return a mask of the cells whose manual label is in EXCLUDED_LABELS, logging their number per dataset and label."""
    ex = obs[label].isin(EXCLUDED_LABELS)
    n = obs[ex].groupby(["dataset", label], observed=True).size()
    log("cells left out (labels " + ", ".join(EXCLUDED_LABELS) + "): " +
        ("; ".join(f"{d} {lab} {k}" for (d, lab), k in n.items()) if len(n) else "none") + f"; total {int(ex.sum())}")
    return ex


def write_10x_h5(path, X, barcodes, gene_ids, gene_names):
    """Write a genes x cells matrix in the Cell Ranger HDF5 layout. X: cells x genes sparse."""
    import h5py
    import scipy.sparse as sp
    M = sp.csc_matrix(X.T)  # genes x cells, CSC = per-cell columns
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=M.data.astype(np.int32), compression="gzip")
        g.create_dataset("indices", data=M.indices.astype(np.int64), compression="gzip")
        g.create_dataset("indptr", data=M.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(M.shape, dtype=np.int32))
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        ft = g.create_group("features")
        ft.create_dataset("id", data=np.array(gene_ids, dtype="S"))
        ft.create_dataset("name", data=np.array(gene_names, dtype="S"))
        ft.create_dataset("feature_type", data=np.array(["Gene Expression"] * len(gene_ids), dtype="S"))
        ft.create_dataset("genome", data=np.array(["GRCh38"] * len(gene_ids), dtype="S"))
        ft.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


# curated markers of the CellTypist Pan Immune Atlas v2 low-hierarchy cell types, read by 02 (markers kept
# beyond the shared genes) and 03 (annotation). A typo and a retired symbol are corrected when read; the
# downloaded file is kept unchanged.
MARKER_FIX = {"SCL40A1": "SLC40A1", "SEPP1": "SELENOP"}


def curated_markers():
    """Return (markers table, {low-hierarchy cell type: [marker genes]})."""
    t = pd.read_excel(PIPELINE_DIR / "reference" / "celltypist_v2_Basic_celltype_information.xlsx", sheet_name="Main")
    m = {r["Low-hierarchy cell types"]: [MARKER_FIX.get(g.strip(), g.strip()) for g in str(r["Curated markers"]).split(",")
                                         if g.strip()] for _, r in t.iterrows()}
    return t, m


# batch for the DE design and the composition formula. Equal to the dataset, except for Perez,
# where it is the authors' processing cohort (Processing_Cohort, deposited per sample; 00_deposited_input.py writes it
# to input/Perez/donors.tsv as processing_cohort): Perez_1 to Perez_4. Perez's batches are unbalanced between SLE and
# controls (batch 1 controls only, batch 2 mostly SLE), so the batch is modelled like a dataset.
def batch_of(obs):
    """Return the batch of each cell of `obs` (columns dataset and sample_id)."""
    donors = pd.read_csv(INPUT_ROOT / "Perez" / "donors.tsv", sep="\t", dtype={"sample_id": str})
    cohort = donors.set_index("sample_id")["processing_cohort"]
    perez = obs["dataset"].astype(str) == "Perez"
    batch = obs["dataset"].astype(str).copy()
    c = obs.loc[perez, "sample_id"].astype(str).map(cohort)
    if c.isna().any():
        raise ValueError(f"Perez samples without a processing cohort: {sorted(obs.loc[perez][c.isna()].sample_id.unique())[:5]}")
    batch[perez] = "Perez_" + c.astype(float).astype(int).astype(str)
    return batch
