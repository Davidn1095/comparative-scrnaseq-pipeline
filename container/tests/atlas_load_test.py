"""Import every Python package the atlas uses; print one 'PKG <name> <version>' line each.
Exits non-zero on any failure. With --gpu, also requires CUDA to be usable by torch and JAX-free scVI."""
import importlib, importlib.metadata as md, sys

mods = {"scvi": "scvi-tools", "scanpy": "scanpy", "anndata": "anndata", "mudata": "mudata",
        "pertpy": "pertpy", "patsy": "patsy", "scib_metrics": "scib-metrics", "numpyro": "numpyro",
        "jax": "jax", "torch": "torch", "decoupler": "decoupler", "pydeseq2": "pydeseq2",
        "numpy": "numpy", "scipy": "scipy", "pandas": "pandas", "h5py": "h5py",
        "openpyxl": "openpyxl", "sklearn": "scikit-learn", "igraph": "igraph", "leidenalg": "leidenalg"}
fail = []
for mod, dist in mods.items():
    try:
        importlib.import_module(mod)
        print("PKG", dist, md.version(dist))
    except Exception as e:  # noqa: BLE001
        print("FAIL", dist, repr(e)); fail.append(dist)
print("PKG python", sys.version.split()[0])
checks = {}
try:
    import pertpy as pt
    sc = pt.tl.Sccoda()
    import inspect
    sig = inspect.signature(sc.prepare)
    checks["pertpy Sccoda.prepare(formula, reference_cell_type='automatic')"] = (
        "formula" in sig.parameters and sig.parameters["reference_cell_type"].default == "automatic")
except Exception as e:  # noqa: BLE001
    checks[f"pertpy Sccoda ({e!r})"] = False
try:
    from scib_metrics.benchmark import Benchmarker  # noqa: F401
    checks["scib_metrics Benchmarker"] = True
except Exception as e:  # noqa: BLE001
    checks[f"scib_metrics Benchmarker ({e!r})"] = False
try:
    import scvi
    checks["scvi.model.SCVI"] = hasattr(scvi.model, "SCVI")
except Exception as e:  # noqa: BLE001
    checks[f"scvi.model.SCVI ({e!r})"] = False
if "--gpu" in sys.argv:
    import torch
    ok = torch.cuda.is_available()
    print("GPU torch", torch.__version__, "cuda", torch.version.cuda, "available", ok,
          torch.cuda.get_device_name(0) if ok else "")
    if ok:
        x = torch.randn(1024, 1024, device="cuda"); _ = (x @ x).sum().item()
    checks["torch CUDA usable"] = ok
for k, v in checks.items():
    print("CHECK OK  " if v else "CHECK FAIL", k)
    if not v: fail.append(k)
if fail:
    print("FAILED:", ", ".join(fail)); sys.exit(1)
print("ALL PYTHON LOAD TESTS PASSED")
