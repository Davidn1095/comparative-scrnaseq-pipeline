"""atlas-2026.10b: CellTypist imports, fetches Immune_All_Low and annotates simulated cells; the
packages of atlas-2026.10 are unchanged (checked against the freeze taken before the install)."""
import os, sys, tempfile, subprocess
import importlib.metadata as md

fail = []
before = set(open("/opt/atlas-2026.10.freeze.txt").read().split())
after = set(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout.split())
changed = sorted(before - after)
print("PKG celltypist", md.version("celltypist"))
print("ADDED", ", ".join(sorted(after - before)))
if changed:
    print("FAIL packages of atlas-2026.10 changed:", changed); fail.append("unchanged base")
try:
    os.environ.setdefault("CELLTYPIST_FOLDER", tempfile.mkdtemp())
    import numpy as np, anndata as ad, scanpy as sc, celltypist
    from celltypist import models
    models.download_models(model="Immune_All_Low.pkl")
    m = models.Model.load("Immune_All_Low.pkl")
    a = ad.AnnData(np.random.default_rng(1).poisson(1.0, (50, len(m.features))).astype("float32"))
    a.var_names = list(m.features)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    p = celltypist.annotate(a, model=m, majority_voting=False)
    ok = len(p.predicted_labels) == 50
    print("CHECK OK  " if ok else "CHECK FAIL", "celltypist annotate Immune_All_Low", m.description.get("version", ""))
    if not ok: fail.append("annotate")
except Exception as e:  # noqa: BLE001
    print("CHECK FAIL celltypist", repr(e)); fail.append("celltypist")
if fail:
    print("FAILED:", ", ".join(fail)); sys.exit(1)
print("ALL ATLAS-B LOAD TESTS PASSED")
