"""Step 00 (deposited-matrix route): write a deposited cell-called matrix in the
standard per-library input layout that Cell Ranger outputs also follow:

  input/<dataset>/libraries.tsv          library, filtered_h5, raw_h5, cells_tsv, n_called, n_called_source
  input/<dataset>/<library>/filtered_feature_bc_matrix.h5   reference gene IDs (GRCh38-2024-A)
  input/<dataset>/<library>/cells.tsv    barcode, sample_id, donor_id (cells with a donor label)
  input/<dataset>/donors.tsv             donor_id, sample_id, condition, sex_reported, age, ...
  input/<dataset>/gene_map.tsv           deposited feature -> reference gene ID and how

Genes are mapped to GRCh38-2024-A by Ensembl ID (version stripped), then, for IDs not in the reference, by a
symbol that names exactly one reference gene not otherwise matched; deposited rows that reach the same
reference gene are summed. A library carries only the reference genes its deposit has (CHANGED 2026-10-29: no
all-zero rows for absent genes, which gene_map.tsv lists); the step 02 merge keeps the genes shared by all
datasets. Counts must be non-negative integers.

Cell Ranger route (CHANGED 2026-10-28): for Perez and GSE253568 the counts are our Cell
Ranger 10.0.0 outputs (home copies under data/<disease>/<dataset>/cellranger/<library>/), and only the
labels come from the deposited objects (Perez: the h5ad's barcode, sample and donor; GSE253568: the
hashtag labels with the hashtag-genotype rule of 00_genotype.py; intersection rule from 2026-11-02). CHANGED 2026-11-01: the
counts are exon-only, rebuilt from Cell Ranger's molecule_info.h5 (molecules with umi_type 1) into
input/<dataset>/<library>/filtered_feature_bc_matrix.h5 and raw_feature_bc_matrix.h5, keeping Cell Ranger's
cell calls; libraries.tsv points to them. Labelled cells that our Cell Ranger does not call are dropped and
counted; n_called is Cell Ranger's called-cell count. Their earlier deposited-matrix test route is
removed so that it cannot be used by mistake; so is GSE157278's (excluded; a public-FASTQ
dataset, which would go through Cell Ranger if it were included again).

Usage: python 00_inputs.py <dataset>
"""
import gzip
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from config import DATA_ROOT, DEV_ROOT, INPUT_ROOT, REF_DIR, log, write_10x_h5

REPO = Path(os.environ.get("ATLAS_ROOT", Path(__file__).resolve().parent.parent))
RAW = {
    "Perez": DATA_ROOT / "SLE/GSE174188",
    "GSE253568": DATA_ROOT / "SjD/GSE253568",
    "GSE135779_adult": DATA_ROOT / "SLE/GSE135779",
    "GSE316482": DATA_ROOT / "SjD/GSE316482",  # CHANGED 2026-10-29: deposited-matrix route
    "HRA003613": DATA_ROOT / "SjD/HRA003613",  # CHANGED 2026-10-29: as above
}
# CHANGED: Cell Ranger outputs of the public-FASTQ route, data/<disease>/<dataset>/cellranger/<library>/; GSE253568 capture1 is GEM well 1 only
CELLRANGER = {"Perez": RAW["Perez"] / "cellranger", "GSE253568": RAW["GSE253568"] / "cellranger"}
CR_DIR = {"GSE253568": {"capture1": "capture1_well1"}}
CONDITION = {"healthy": "control", "normal": "control", "HD": "control", "sle": "SLE", "SLE": "SLE",
             "systemic lupus erythematosus": "SLE", "sjs": "SjD",
             "HC": "control", "pSS": "SjD", "XN": "SjD"}  # CHANGED: GSE316482 and HRA003613 labels


def reference_genes():
    import re
    ids, names = [], []
    with gzip.open(REF_DIR / "genes" / "genes.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            ids.append(re.search(r'gene_id "([^"]+)"', p[8]).group(1))
            names.append(re.search(r'gene_name "([^"]+)"', p[8]).group(1))
    return pd.DataFrame({"gene_id": ids, "gene_name": names})


def map_to_reference(dep_ids, dep_names, ref):
    """Map deposited genes to GRCh38-2024-A. Return (n_kept x n_dep aggregation matrix, the
    reference genes kept, gene_map table). CHANGED 2026-10-29: Ensembl ID first (version stripped); a
    deposited gene whose ID is not in the reference is then mapped by its symbol, where the symbol names
    exactly one reference gene that no deposited ID matched and no other deposited symbol claims. Deposited
    IDs reaching the same reference gene are summed. Only reference genes present in the deposit are kept
    (no all-zero rows for absent genes); the step 02 merge keeps the genes shared by every dataset."""
    ref_idx = pd.Series(np.arange(len(ref)), index=ref.gene_id)
    target = pd.Series(dep_ids).str.replace(r"\.\d+$", "", regex=True).map(ref_idx)
    miss = target.isna().values
    u = ref.gene_name.value_counts()
    uniq = ref[ref.gene_name.map(u) == 1]
    by_sym = pd.Series(dep_names).map(pd.Series(uniq.index.values, index=uniq.gene_name))
    cand = miss & by_sym.notna().values & ~by_sym.isin(target[~miss]).values
    cand &= ~by_sym.duplicated(keep=False).values  # two deposited symbols on one gene: neither is used
    target[cand] = by_sym[cand]
    ok = target.notna().values
    by = np.where(cand, "symbol", np.where(ok, "ensembl", "unmapped"))
    kept = np.unique(target[ok].astype(int).values)  # reference order
    row = pd.Series(np.arange(len(kept)), index=kept)
    A = sp.csr_matrix((np.ones(ok.sum()), (row[target[ok].astype(int).values].values, np.where(ok)[0])),
                      shape=(len(kept), len(dep_ids)))
    gm = pd.DataFrame({"deposited_id": dep_ids, "deposited_name": dep_names, "mapped_by": by,
                       "gene_id": np.where(ok, ref.gene_id.values[target.fillna(0).astype(int)], "")})
    dup = gm.gene_id[(gm.gene_id != "") & gm.gene_id.duplicated(keep=False)]
    gm["summed"] = gm.gene_id.isin(dup)
    absent = ref[~ref.gene_id.isin(gm.gene_id)]  # recorded only; not written to the matrix
    gm = pd.concat([gm, pd.DataFrame({"deposited_id": "", "deposited_name": "", "mapped_by": "absent_from_deposit",
                                      "gene_id": absent.gene_id.values, "summed": False})])
    return A, ref.iloc[kept].reset_index(drop=True), gm


def check_counts(ds, lib, X):
    """Deposited counts must be non-negative integers (raw UMI counts)."""
    d = X.data
    if d.size and (d.min() < 0 or not np.array_equal(d, np.round(d))):
        sys.exit(f"{ds} {lib}: counts are not non-negative integers")


def read_mtx_dir(mtx, barcodes, features):
    X = scipy.io.mmread(mtx).tocsc()  # genes x cells
    bc = pd.read_csv(barcodes, header=None, sep="\t")[0].astype(str).values
    ft = pd.read_csv(features, header=None, sep="\t", dtype=str)
    return X, bc, ft


def write_library(ds, lib, X_genes_cells, barcodes, cells, genes):  # genes: the dataset's mapped reference genes
    if len(set(barcodes)) != len(barcodes):
        sys.exit(f"{ds} {lib}: duplicated barcodes within the library")
    check_counts(ds, lib, X_genes_cells)
    d = INPUT_ROOT / ds / lib
    write_10x_h5(d / "filtered_feature_bc_matrix.h5", X_genes_cells.T, barcodes, genes.gene_id.values,
                 genes.gene_name.values)
    cells.to_csv(d / "cells.tsv", sep="\t", index=False)
    return {"library": lib, "filtered_h5": str(d / "filtered_feature_bc_matrix.h5"), "raw_h5": "",
            "cells_tsv": str(d / "cells.tsv"), "n_called": len(barcodes)}


def exon_only_matrices(ds, lib, cr, d, called, ref):
    """CHANGED 2026-11-01: exon-only counts from Cell Ranger's molecule_info.h5, which holds every
    counted molecule with umi_type 1 (transcriptomic) or 0 (counted only because Cell Ranger 10 includes introns).
    Molecules with umi_type 1 are counted per barcode and gene; the filtered matrix keeps Cell Ranger's called
    barcodes (pass_filter, checked against its filtered matrix), the raw matrix every barcode with an exonic UMI
    or called. Written to input/<dataset>/<library>/."""
    with h5py.File(cr / "molecule_info.h5", "r") as f:
        if not np.array_equal(f["features/id"][:].astype(str), ref.gene_id.values):
            sys.exit(f"{ds} {lib}: molecule_info features differ from the reference gene list")
        if not (f["features/feature_type"][:].astype(str) == "Gene Expression").all():
            sys.exit(f"{ds} {lib}: molecule_info holds features other than Gene Expression")
        gg = np.unique(f["gem_group"][:])
        if len(gg) != 1:
            sys.exit(f"{ds} {lib}: {len(gg)} GEM groups in molecule_info")
        bc = np.char.add(f["barcodes"][:].astype(str), f"-{gg[0]}")
        pf = np.unique(f["barcode_info/pass_filter"][:, 0])
        ut = f["umi_type"][:]
        exon = ut == 1
        bi, fi = f["barcode_idx"][:][exon], f["feature_idx"][:][exon]
    if not np.array_equal(np.sort(bc[pf]), np.sort(called)):
        sys.exit(f"{ds} {lib}: molecule_info pass_filter barcodes differ from Cell Ranger's filtered matrix")
    M = sp.csr_matrix((np.ones(len(bi), np.int32), (bi, fi)), shape=(len(bc), len(ref)))  # duplicates summed
    tot = np.asarray(M.sum(axis=1)).ravel()
    is_called = np.zeros(len(bc), bool)
    is_called[pf] = True
    raw = np.where((tot > 0) | is_called)[0]
    write_10x_h5(d / "filtered_feature_bc_matrix.h5", M[pf], bc[pf], ref.gene_id.values, ref.gene_name.values)
    write_10x_h5(d / "raw_feature_bc_matrix.h5", M[raw], bc[raw], ref.gene_id.values, ref.gene_name.values)
    n_called_exon = int(tot[pf].sum())
    log(ds, lib, f"exon-only: {int(exon.sum())} of {len(ut)} UMIs ({exon.mean():.3f}); called cells keep "
        f"{n_called_exon} UMIs, median {np.median(tot[pf]):.0f}; raw barcodes {len(raw)}")
    return int(exon.sum()), len(ut)


def cellranger_library(ds, lib, cells, ref):
    """CHANGED: Cell Ranger route. Counts come from our Cell Ranger outputs, exon-only (CHANGED
    2026-11-01, see exon_only_matrices); cells.tsv holds the labelled cells that our Cell Ranger also called.
    cells: barcode (10x form), sample_id, donor_id."""
    cr = CELLRANGER[ds] / CR_DIR.get(ds, {}).get(lib, lib)
    with h5py.File(cr / "filtered_feature_bc_matrix.h5", "r") as f:
        called = f["matrix/barcodes"][:].astype(str)
        ids = f["matrix/features/id"][:].astype(str)
    if not np.array_equal(ids, ref.gene_id.values):
        sys.exit(f"{ds} {lib}: Cell Ranger features differ from the reference gene list")
    est = int(pd.read_csv(cr / "metrics_summary.csv", thousands=",")["Estimated Number of Cells"].iloc[0])
    if est != len(called):
        sys.exit(f"{ds} {lib}: {len(called)} filtered barcodes vs {est} estimated cells in metrics_summary.csv")
    if cells.barcode.duplicated().any():
        sys.exit(f"{ds} {lib}: duplicated labelled barcodes")
    inc = cells.barcode.isin(set(called)).values
    d = INPUT_ROOT / ds / lib
    d.mkdir(parents=True, exist_ok=True)
    n_exon, n_umis = exon_only_matrices(ds, lib, cr, d, called, ref)  # CHANGED 2026-11-01: exon-only
    cells[inc].to_csv(d / "cells.tsv", sep="\t", index=False)
    log(ds, lib, f"{len(called)} called by Cell Ranger; {len(cells)} labelled, {int((~inc).sum())} of them not called (dropped)")
    return {"library": lib, "filtered_h5": str(d / "filtered_feature_bc_matrix.h5"),
            "raw_h5": str(d / "raw_feature_bc_matrix.h5"), "cells_tsv": str(d / "cells.tsv"),
            "n_called": len(called), "n_labelled": len(cells), "n_labelled_not_called": int((~inc).sum()),
            "umis_exonic": n_exon, "umis_all": n_umis}


def perez(ref):
    f = h5py.File(RAW["Perez"] / "raw" / "GSE174188.h5ad", "r")
    obs = f["obs"]

    def col(k):
        g = obs[k]
        if isinstance(g, h5py.Group):
            return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))
        return g[:].astype(str)

    meta = pd.DataFrame({k: np.asarray(col(k)) for k in
                         ["library_uuid", "sample_uuid", "donor_id", "disease", "disease_state", "sex",
                          "development_stage", "Processing_Cohort", "ind_cov"]})
    # CELLxGENE appends merge suffixes to the 10x barcode (AAAC...-1-1-0-0-...); keep the 10x barcode
    meta["barcode"] = pd.Series(obs[obs.attrs["_index"]][:].astype(str)).str.extract(r"^([ACGT]+-\d+)")[0].values
    f.close()
    # CHANGED: counts from our Cell Ranger outputs; the h5ad gives only the labels, so its
    # raw/X and the gene mapping are no longer read (the deposited-matrix test route is removed)
    libs = []
    for lib, idx in meta.groupby("library_uuid", observed=True).indices.items():
        cells = pd.DataFrame({"barcode": meta.barcode.values[idx], "sample_id": meta.sample_uuid.values[idx],
                              "donor_id": meta.donor_id.values[idx]})
        libs.append(cellranger_library("Perez", lib, cells, ref))
    s = meta.drop_duplicates("sample_uuid")
    donors = pd.DataFrame({"donor_id": s.donor_id, "sample_id": s.sample_uuid,
                           "condition": s.disease.map(CONDITION), "sex_reported": s.sex,
                           "age": s.development_stage, "disease_state": s.disease_state,
                           "processing_cohort": s.Processing_Cohort, "author_id": s.ind_cov})
    return libs, donors, None, "Cell Ranger 10.0.0 molecule_info, exon-only"  # CHANGED: Cell Ranger route; exon-only 2026-11-01


# CHANGED: GSE157278's deposited-matrix test route removed (2026-10-28); excluded, and as a
# public-FASTQ dataset it would go through Cell Ranger if it were included again


def gse253568(ref):
    s = pd.read_csv(Path(__file__).resolve().parent / "inputs/GSE253568_samples.tsv", sep="\t", dtype=str)
    gl_file = DEV_ROOT / "labels" / "GSE253568" / "donor_labels.tsv"  # CHANGED 2026-11-02: per run root
    if not gl_file.exists():
        sys.exit(f"GSE253568: {gl_file} missing; run 00_genotype.py labels first")
    gl = pd.read_csv(gl_file, sep="\t", dtype=str)
    # CHANGED: counts from our Cell Ranger outputs, one library per capture (capture1 = GEM well 1);
    # labels are the deposited hashtag labels (barcode <10x barcode>_<capture>), and the
    # deposited matrices and the symbol mapping are no longer read (the test route is removed)
    gl["cr_barcode"] = gl.barcode.str.rsplit("_", n=1).str[0]
    libs = []
    for lib, L in sorted(gl.groupby("capture")):
        cells = pd.DataFrame({"barcode": L.cr_barcode.values, "sample_id": L.sample_id.values,
                              "donor_id": L.sample_id.values})
        # labelled cells not called by our Cell Ranger are dropped first (counted in cellranger_library);
        # then the hashtag-genotype rule: keep a cell only when vireo agrees (intersection rule, 2026-11-02)
        row = cellranger_library("GSE253568", lib, cells, ref)
        d = pd.read_csv(row["cells_tsv"], sep="\t", dtype=str)
        kept = set(L.cr_barcode[L.kept == "True"])
        rule = d.barcode.isin(kept)
        d[rule].to_csv(row["cells_tsv"], sep="\t", index=False)
        row["n_removed_genotype_rule"] = int((~rule).sum())
        log("GSE253568", lib, f"{int(rule.sum())} kept by the hashtag-genotype rule, {int((~rule).sum())} removed")
        libs.append(row)
    donors = pd.DataFrame({"donor_id": s.donor_id, "sample_id": s.sample_id,
                           "condition": s.condition.map(CONDITION), "sex_reported": "", "age": ""})
    return libs, donors, None, "Cell Ranger 10.0.0 molecule_info, exon-only"  # CHANGED: Cell Ranger route; exon-only 2026-11-01


def gse135779_adult(ref):
    r = RAW["GSE135779_adult"] / "raw/extracted"
    s = pd.read_csv(Path(__file__).resolve().parent / "inputs/GSE135779_adult_samples.tsv", sep="\t", dtype=str)
    libs, gm = [], None
    for g in s.gsm:
        stem = list(r.glob(f"{g}_*_matrix.mtx.gz"))[0].name.replace("_matrix.mtx.gz", "")
        X, bc, ft = read_mtx_dir(r / f"{stem}_matrix.mtx.gz", r / f"{stem}_barcodes.tsv.gz", r / f"{stem}_genes.tsv.gz")
        if gm is None:
            A, genes, gm = map_to_reference(ft[0].values, ft[1].values, ref)
            ids0 = ft[0].values
        elif not np.array_equal(ft[0].values, ids0):
            sys.exit(f"GSE135779: {g} has a different feature list")
        donor = s.title[s.gsm == g].iloc[0].split(" ")[0]
        cells = pd.DataFrame({"barcode": bc, "sample_id": donor, "donor_id": donor})
        libs.append(write_library("GSE135779_adult", g, (A @ X).tocsc(), bc, cells, genes))
    donors = pd.DataFrame({"donor_id": s.title.str.split(" ").str[0], "sample_id": s.title.str.split(" ").str[0],
                           "condition": s.geo_group.map(CONDITION), "sex_reported": "", "age": s.age, "gsm": s.gsm})
    return libs, donors, gm, "deposited Cell Ranger 2.1.0 filtered matrix"


def gse316482(ref):  # CHANGED 2026-10-29: deposited Cell Ranger filtered matrices, one library per sample
    r = RAW["GSE316482"] / "raw/extracted"
    s = pd.read_csv(Path(__file__).resolve().parent / "inputs/GSE316482_samples.tsv", sep="\t", dtype=str)
    libs, gm = [], None
    for _, S in s.iterrows():
        stem = f"{S.gsm}_{S.donor}_scRNAseq_all_cells"  # the RNA matrix of all cells (TCR and BCR files are not read)
        X, bc, ft = read_mtx_dir(r / f"{stem}_matrix.mtx.gz", r / f"{stem}_barcodes.tsv.gz", r / f"{stem}_features.tsv.gz")
        check_counts("GSE316482", S.gsm, X)
        if gm is None:
            A, genes, gm = map_to_reference(ft[0].values, ft[1].values, ref)
            ids0 = ft[0].values
        elif not np.array_equal(ft[0].values, ids0):
            sys.exit(f"GSE316482: {S.gsm} has a different feature list")
        cells = pd.DataFrame({"barcode": bc, "sample_id": S.donor, "donor_id": S.donor})
        libs.append(write_library("GSE316482", S.gsm, (A @ X).tocsc(), bc, cells, genes))
    n = sum(L["n_called"] for L in libs)
    if n != 275747:  # cells deposited (GEO), before the authors' QC
        sys.exit(f"GSE316482: {n} cells, expected 275747")
    donors = pd.DataFrame({"donor_id": s.donor, "sample_id": s.donor, "condition": s.geo_group.map(CONDITION),
                           "sex_reported": "", "age": "", "gsm": s.gsm})
    return libs, donors, gm, "deposited Cell Ranger filtered matrix"


def hra003613(ref):  # CHANGED 2026-10-29: Zenodo 10884425 h5ad, PBMC cells, raw counts in .raw
    f = h5py.File(RAW["HRA003613"] / "raw/SG_PBMC_immune_inte_seurat_rpca_pc20_neighbor50.h5ad", "r")

    def col(g, k):
        o = g[k]
        if isinstance(o, h5py.Group):
            return o["categories"][:].astype(str)[o["codes"][:]]
        return o[:].astype(str) if o.dtype.kind in "SO" else o[:]

    obs = f["obs"]
    names = obs[obs.attrs["_index"]][:].astype(str)
    tissue, sample = col(obs, "tissue"), col(obs, "sample")
    n_count, n_feat = col(obs, "nCount_RNA"), col(obs, "nFeature_RNA")
    rv = f["raw/var"][f["raw/var"].attrs["_index"]][:].astype(str)
    if not np.array_equal(rv, f["var"][f["var"].attrs["_index"]][:].astype(str)):
        sys.exit("HRA003613: raw/var differs from var")
    # gene names are the GRCh38-2020-A symbols made unique in Cell Ranger's order (R make.unique): recover the Ensembl IDs
    f20 = pd.read_csv(Path(__file__).resolve().parent / "reference/GRCh38-2020-A_features.tsv.gz", sep="\t", header=None, dtype=str)
    seen, mu = {}, []
    for x in f20[1]:
        mu.append(x if x not in seen else f"{x}.{seen[x]}")
        seen[x] = seen.get(x, 0) + 1
    if not np.array_equal(rv, np.array(mu)):
        sys.exit("HRA003613: gene names are not the GRCh38-2020-A list")
    A, genes, gm = map_to_reference(f20[0].values, f20[1].values, ref)
    g = f["raw/X"]
    R = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
    f.close()
    pb = tissue == "PBMC"
    R, names, sample, n_count, n_feat = R[pb], names[pb], sample[pb], n_count[pb], n_feat[pb]
    check_counts("HRA003613", "all", R)
    if not (np.array_equal(np.asarray(R.sum(1)).ravel(), n_count) and np.array_equal(np.diff(R.indptr), n_feat)):
        sys.exit("HRA003613: raw counts do not reproduce nCount_RNA and nFeature_RNA")
    if pb.sum() != 55203:
        sys.exit(f"HRA003613: {pb.sum()} PBMC cells, expected 55203")
    bc = pd.Series(names).str.replace(r"_\d+_\d+$", "", regex=True).values  # integration suffix removed: 10x barcode
    libs = []
    for smp in sorted(set(sample)):
        i = np.where(sample == smp)[0]
        cells = pd.DataFrame({"barcode": bc[i], "sample_id": smp, "donor_id": smp})
        libs.append(write_library("HRA003613", smp, (A @ R[i].T).tocsc(), bc[i], cells, genes))
    s = pd.Series(sorted(set(sample)))
    donors = pd.DataFrame({"donor_id": s, "sample_id": s, "condition": s.str.split("_").str[0].map(CONDITION),
                           "sex_reported": "female", "age": ""})  # Supplemental Table 2: all 8 PBMC donors female
    return libs, donors, gm, "deposited cells after the authors' QC and doublet removal"


if __name__ == "__main__":
    ds = sys.argv[1]
    ref = reference_genes()
    routes = {"Perez": perez, "GSE253568": gse253568, "GSE135779_adult": gse135779_adult,  # CHANGED: no GSE157278
              "GSE316482": gse316482, "HRA003613": hra003613}
    if ds not in routes:
        sys.exit(f"{ds}: no input route in this script (GSE157278 is excluded; a public-FASTQ "
                 "dataset needs its Cell Ranger route)")
    libs, donors, gm, source = routes[ds](ref)
    out = INPUT_ROOT / ds
    lt = pd.DataFrame(libs)
    lt["n_called_source"] = source
    lt.to_csv(out / "libraries.tsv", sep="\t", index=False)
    if donors.condition.isna().any():
        sys.exit(f"{ds}: unmapped condition label(s)")
    donors.to_csv(out / "donors.tsv", sep="\t", index=False)
    if gm is None:  # CHANGED: Cell Ranger route, counts already on the reference; no gene map
        log(f"{ds}: {len(lt)} libraries, {lt.n_called.sum()} cells called by Cell Ranger, "
            f"{lt.n_labelled.sum()} labelled, {lt.n_labelled_not_called.sum()} labelled but not called (dropped), "
            f"{donors.donor_id.nunique()} donors, {donors.sample_id.nunique()} samples")
        sys.exit(0)
    gm.to_csv(out / "gene_map.tsv", sep="\t", index=False)
    k = gm.mapped_by.value_counts()
    log(f"{ds}: {len(lt)} libraries, {lt.n_called.sum()} cells, {donors.donor_id.nunique()} donors, "
        f"{donors.sample_id.nunique()} samples; deposited features {int((gm.mapped_by != 'absent_from_deposit').sum())}: "
        f"mapped by Ensembl ID {k.get('ensembl', 0)}, by symbol {k.get('symbol', 0)}, unmapped {k.get('unmapped', 0)}, "
        f"summed {int(gm.summed.sum())}; reference genes kept {gm.gene_id[gm.mapped_by.isin(['ensembl', 'symbol'])].nunique()}, "
        f"absent from the deposit {k.get('absent_from_deposit', 0)}")
