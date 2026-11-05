"""GSE253568 donor labels checked by genotype (intersection rule, 2026-11-02). Two parts, both in the atlas
container (cellsnp-lite 1.2.3 and vireoSNP 0.5.9 are in the image):

1. Genotypes of one Cell Ranger library (one job per library; 16 CPUs, about 1-3.5 h, mostly cellsnp-lite):
   cellsnp-lite on the library's BAM at the 1000 Genomes SNPs (allele frequency >= 5%, GRCh38), over the cells Cell Ranger
   called, then vireo once per donor number N. Input: $DATA_ROOT/SjD/GSE253568/cellranger/<library>/ (cellranger.sh with
   CREATE_BAM=true). Output: $PIPELINE_OUT/genotype/GSE253568/<library>/ (barcodes.tsv.gz, cellsnp/, vireo_N<N>/,
   vireo_N<N>.log). Libraries: capture1_well1 and capture2-5.
2. Labels (seconds): among the cells the authors deposited with a hashtag label, a cell keeps that label only when vireo
   assigns it to the same donor; it is removed when vireo leaves it unassigned, assigns another donor or calls a doublet,
   and when our Cell Ranger run does not call it (not genotyped). Cells vireo assigns that the authors did not deposit are
   not added. Per capture the vireo run with N = the number of hashtag donors in the capture is used; capture1 is two GEM
   wells, and only well 1, the deposited one, is used (capture1_well1). Each vireo donor is mapped to the hashtag label held
   by most of its singlets; the mapping must be one-to-one, otherwise the script stops. Status per cell: agree (kept) or
   unassigned | disagree | doublet | not_genotyped (removed); column 'kept'.
   Output, read by 00_inputs.py: $PIPELINE_OUT/labels/GSE253568/donor_labels.tsv and donor_labels_summary.tsv (per run
   root, so that a run with another rule keeps its own labels)
     capture, barcode (deposited form, <10x barcode>_<capture>), sample_id (hashtag), vireo_donor,
     vireo_label (mapped), status (agree | disagree | doublet | unassigned | not_genotyped), kept

Usage: python 00_genotype.py <library> [N ...]   part 1 for one library; N defaults to 1 2 3 4 5
       python 00_genotype.py labels             part 2, after part 1 for every library
"""
import gzip
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import pandas as pd

from config import DATA_ROOT, DEV_ROOT, log

CR = DATA_ROOT / "SjD" / "GSE253568" / "cellranger"
GT = DEV_ROOT / "genotype" / "GSE253568"
# 1000 Genomes phase 3 SNPs with allele frequency >= 5% on GRCh38, as cellsnp-lite distributes them
PANEL = Path(os.environ.get("SNP_PANEL", REPO / "refdata" / "genome1K.phase3.SNP_AF5e2.chr1toX.GRCh38.vcf.gz"))
THREADS = 16
REPO = Path(os.environ.get("ATLAS_ROOT", Path(__file__).resolve().parent.parent))
RAW = DATA_ROOT / "SjD/GSE253568"
GENO_DIR = {"capture1": "capture1_well1"}  # capture1: the deposited GEM well only
KEEP = ["agree"]  # intersection rule (2026-11-02)


def genotype(lib, n_list):
    cr, gt = CR / lib, GT / lib
    gt.mkdir(parents=True, exist_ok=True)
    log(subprocess.run(["cellsnp-lite", "--version"], capture_output=True, text=True).stdout.strip())
    import vireoSNP
    log("vireoSNP", vireoSNP.__version__)
    # barcode list: the cells Cell Ranger called, from the filtered .h5
    with h5py.File(cr / "filtered_feature_bc_matrix.h5", "r") as f:
        bc = f["matrix/barcodes"][:].astype(str)
    with gzip.open(gt / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(bc) + "\n")
    t0 = time.time()
    subprocess.run(["cellsnp-lite", "-s", str(cr / "possorted_genome_bam.bam"), "-b", str(gt / "barcodes.tsv.gz"),
                    "-R", str(PANEL), "-O", str(gt / "cellsnp"), "-p", str(THREADS), "--minMAF", "0.1",
                    "--minCOUNT", "20", "--cellTAG", "CB", "--UMItag", "UB", "--gzip"], check=True)
    with gzip.open(gt / "cellsnp" / "cellSNP.base.vcf.gz", "rt") as f:
        n_snp = sum(1 for line in f if not line.startswith("#"))
    log(f"cellsnp-lite {(time.time() - t0) / 60:.0f} min; SNPs kept: {n_snp}")
    for n in n_list:
        t = time.time()
        with open(gt / f"vireo_N{n}.log", "w") as out:
            r = subprocess.run(["vireo", "-c", str(gt / "cellsnp"), "-N", str(n), "-o", str(gt / f"vireo_N{n}"),
                                "--randSeed", "1", "-p", str(THREADS)], stdout=out, stderr=subprocess.STDOUT)
        if r.returncode:
            log(f"vireo N={n} exit {r.returncode}")
        tail = [s for s in (gt / f"vireo_N{n}.log").read_text().splitlines()
                if re.search("elbo|lower bound|loglik|final", s, re.I)][-3:]
        log(f"vireo N={n} {(time.time() - t) / 60:.0f} min", *tail)
    log(f"total {(time.time() - t0) / 60:.0f} min")


def labels():
    out_dir = DEV_ROOT / "labels" / "GSE253568"
    out_dir.mkdir(parents=True, exist_ok=True)
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
        f = GT / GENO_DIR.get(cap, cap) / f"vireo_N{n}" / "donor_ids.tsv"
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
    res["kept"] = res.status.isin(KEEP)
    res[["capture", "barcode", "sample_id", "vireo_donor", "vireo_label", "status", "kept"]].to_csv(
        out_dir / "donor_labels.tsv", sep="\t", index=False)
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
    s.to_csv(out_dir / "donor_labels_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "labels":
        labels()
    else:
        genotype(sys.argv[1], [int(n) for n in sys.argv[2:]] or [1, 2, 3, 4, 5])
