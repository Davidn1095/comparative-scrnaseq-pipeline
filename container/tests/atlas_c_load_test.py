"""atlas-2026.10c: rpy2 and anndata2ri import; R runs from Python and DecontX (seed 12345) and
scDblFinder (set.seed(123)) work on a simulated matrix passed through rpy2 as the sc-best-practices
notebooks pass it; the packages of atlas-2026.10b are unchanged."""
import subprocess, sys
import importlib.metadata as md

fail = []
before = set(open("/opt/atlas-2026.10b.freeze.txt").read().split())
after = set(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout.split())
print("ADDED", ", ".join(sorted(after - before)))
if before - after:
    print("FAIL packages of atlas-2026.10b changed:", sorted(before - after)); fail.append("unchanged base")
try:
    import numpy as np
    import anndata2ri  # noqa: F401
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.conversion import localconverter
    print("PKG rpy2", md.version("rpy2"), "| anndata2ri", md.version("anndata2ri"), "|", ro.r("R.version.string")[0])
    rng = np.random.default_rng(1)
    m = rng.negative_binomial(2, 0.3, size=(300, 400)).astype(np.float64)   # genes x cells
    with localconverter(ro.default_converter + numpy2ri.converter):
        ro.globalenv["m"] = m
    ro.r("""
      suppressMessages({library(SingleCellExperiment); library(decontX); library(scDblFinder)})
      rownames(m) <- paste0("g", seq_len(nrow(m))); colnames(m) <- paste0("c", seq_len(ncol(m)))
      sce <- SingleCellExperiment(list(counts = m))
      d <- decontX(sce, seed = 12345, verbose = FALSE)
      set.seed(123)
      s <- scDblFinder(SingleCellExperiment(list(counts = m)), dbr = 0.01)
      ok <- length(d$decontX_contamination) == 400 && length(s$scDblFinder.class) == 400
    """)
    ok = bool(ro.r("ok")[0])
    print("CHECK OK  " if ok else "CHECK FAIL", "DecontX and scDblFinder through rpy2")
    if not ok: fail.append("rpy2 R calls")
except Exception as e:  # noqa: BLE001
    print("CHECK FAIL rpy2", repr(e)); fail.append("rpy2")
if fail:
    print("FAILED:", ", ".join(fail)); sys.exit(1)
print("ALL ATLAS-C LOAD TESTS PASSED")
