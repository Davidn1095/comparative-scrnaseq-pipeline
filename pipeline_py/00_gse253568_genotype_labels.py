"""GSE253568 donor labels (intersection rule, 2026-11-02): among the cells the authors deposited with a
hashtag label, a cell keeps that label only when vireo assigns it to the same donor; it is removed when vireo leaves
it unassigned, assigns another donor or calls a doublet, and when our Cell Ranger run does not call it (not
genotyped). Cells vireo assigns that the authors did not deposit are not added.

Per capture: vireo is run with N = the number of hashtag donors in the capture (pipeline_py/genotype.sh;
cellsnp-lite on the Cell Ranger cells, 1000 Genomes SNPs). Capture1 is two GEM wells; only well 1, the
deposited one, is used, with its own Cell Ranger run and genotypes (capture1_well1). Each vireo donor is mapped to the
hashtag label held by most of its singlets; the mapping must be one-to-one, otherwise the script
stops. Status per cell: agree | unassigned (kept) or disagree | doublet | not_genotyped (removed); column 'kept'.

Output: $PIPELINE_OUT/labels/GSE253568/donor_labels.tsv and donor_labels_summary.tsv (per run root, so that a run
  with another rule keeps its own labels)
  capture, barcode (deposited form, <10x barcode>_<capture>), sample_id (hashtag), vireo_donor,
  vireo_label (mapped), status (agree | disagree | doublet | unassigned | not_genotyped)
Usage: python 00_gse253568_genotype_labels.py
"""
import gzip
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(os.environ.get("ATLAS_ROOT", Path(__file__).resolve().parent.parent))
OUT = Path(os.environ.get("PIPELINE_OUT", REPO / "out"))
RAW = Path(os.environ.get("DATA_ROOT", REPO / "data")) / "SjD/GSE253568"
GENO_DIR = {"capture1": "capture1_well1"}
LABELS = OUT / "labels/GSE253568"
LABELS.mkdir(parents=True, exist_ok=True)

smp = pd.read_csv(Path(__file__).resolve().parent / "inputs" / "GSE253568_samples.tsv", sep="\t", dtype=str)
lab = []
for s in smp.sample_id:
    with gzip.open(RAW / f"raw/extracted/{s}_barcodes.tsv.gz", "rt") as f:
        for line in f:
            b = line.strip()
            lab.append((f"capture{b.rsplit('_', 1)[1]}", b, b.rsplit("_", 1)[0], s))
lab = pd.DataFrame(lab, columns=["capture", "barcode", "cr_barcode", "sample_id"])
if lab.barcode.duplicated().any():
    sys.exit("a deposited barcode carries two hashtag labels")

out, summary = [], []
for cap, L in lab.groupby("capture", sort=True):
    n = L.sample_id.nunique()
    # capture1 is two GEM wells; the deposit holds well 1 only, genotyped on its own
    gdir = GENO_DIR.get(cap, cap)
    f = OUT / f"genotype/GSE253568/{gdir}/vireo_N{n}/donor_ids.tsv"  # genotype.sh
    if not f.exists():
        sys.exit(f"{cap}: {f} missing (vireo with N = {n} hashtag donors)")
    v = pd.read_csv(f, sep="\t", usecols=["cell", "donor_id"]).rename(columns={"cell": "cr_barcode", "donor_id": "vireo_donor"})
    L = L.merge(v, on="cr_barcode", how="left")
    sing = L[L.vireo_donor.notna() & ~L.vireo_donor.isin(["doublet", "unassigned"])]
    major = sing.groupby("vireo_donor").sample_id.agg(lambda x: x.value_counts().index[0])
    if major.nunique() != len(major) or len(major) != n:
        sys.exit(f"{cap}: vireo donors do not map one-to-one to the {n} hashtag donors: {major.to_dict()}")
    L["vireo_label"] = L.vireo_donor.map(major)
    L["status"] = "disagree"
    L.loc[L.vireo_label == L.sample_id, "status"] = "agree"
    L.loc[L.vireo_donor == "doublet", "status"] = "doublet"
    L.loc[L.vireo_donor == "unassigned", "status"] = "unassigned"
    L.loc[L.vireo_donor.isna(), "status"] = "not_genotyped"
    out.append(L)
    t = L.groupby(["sample_id", "status"]).size().unstack(fill_value=0)
    t.insert(0, "capture", cap)
    summary.append(t)

res = pd.concat(out)
KEEP = ["agree"]  # intersection rule: a cell is kept only when vireo and its hashtag agree
res["kept"] = res.status.isin(KEEP)
res[["capture", "barcode", "sample_id", "vireo_donor", "vireo_label", "status", "kept"]].to_csv(
    LABELS / "donor_labels.tsv", sep="\t", index=False)
s = pd.concat(summary).fillna(0)
cols = [c for c in ["agree", "disagree", "doublet", "unassigned", "not_genotyped"] if c in s]
s[cols] = s[cols].astype(int)
s["labelled"] = s[cols].sum(axis=1)
s["kept"] = s[[c for c in KEEP if c in s]].sum(axis=1)
s["kept_fraction"] = (s["kept"] / s["labelled"]).round(3)
s = s.reset_index()[["capture", "sample_id"] + cols + ["labelled", "kept", "kept_fraction"]]
pd.set_option("display.width", 250)
print(s.to_string(index=False))
print(f"\ntotal: {int(s.kept.sum())} of {int(s.labelled.sum())} labelled cells kept ({s.kept.sum() / s.labelled.sum():.3f})")
s.to_csv(LABELS / "donor_labels_summary.tsv", sep="\t", index=False)
