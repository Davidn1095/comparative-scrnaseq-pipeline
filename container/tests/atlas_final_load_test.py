"""atlas-final-2026.10 load test. Checks that the image has the versions of atlas-2026.10c plus cellsnp-lite 1.2.3
(htslib 1.24) and vireoSNP 0.5.9, against the locks exported from atlas-2026.10c (/opt/atlas-locks), then runs the
earlier load tests and a genotype smoke test. Versions are compared, not bytes: packages compiled at install time
(R packages from source, Python sdists such as ete4 and rpy2-rinterface) are not byte-identical between builds.
Exits 1 on any failure.
  1. pip freeze: equal to the lock plus vireoSNP==0.5.9, nothing else added, removed or changed
  2. R installed.packages (Package, Version, LibPath, Priority, Repository, RemoteSha): equal
  3. dpkg packages and versions: equal
  4. /opt entries: equal plus atlas-locks and cellsnp
  5. Python, R and Ubuntu versions: equal
  6. cellsnp-lite 1.2.3 (htslib 1.24), every library it loads resolved, non-system ones inside /opt/cellsnp
  7. vireo separates two simulated donors (simulated cellsnp-lite output, N = 2, randSeed 1)
  8. R ignores ~/.Renviron and ~/.Rprofile (R_ENVIRON_USER, R_PROFILE_USER), Python ignores the user site
  9. the atlas-2026.10, 10b and 10c load tests pass"""
import gzip
import os
import subprocess
import sys
import tempfile

L = "/opt/atlas-locks"
fail = []


def check(name, ok, detail=""):
    print(("CHECK OK   " if ok else "CHECK FAIL ") + name + (f" | {detail}" if detail else ""))
    if not ok:
        fail.append(name)


def run(cmd, **kw):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, **kw)
    except FileNotFoundError as e:  # a missing tool fails its check instead of stopping the test
        return subprocess.CompletedProcess(cmd, 127, "", str(e))


# 1. pip freeze
lock = set(open(f"{L}/atlas-2026.10c.pip.txt").read().split())
now = set(run([sys.executable, "-m", "pip", "list", "--format=freeze", "--disable-pip-version-check"]).stdout.split())
check("pip freeze = 10c + vireoSNP", now - lock == {"vireoSNP==0.5.9"} and not lock - now,
      f"added {sorted(now - lock)}, missing {sorted(lock - now)}")

# 2. R installed.packages
r = run(["Rscript", "-e", """
ip <- installed.packages(fields = c("Repository", "RemoteType", "RemoteSha"))
df <- data.frame(Package = ip[, "Package"], Version = ip[, "Version"], LibPath = ip[, "LibPath"],
                 Priority = ifelse(is.na(ip[, "Priority"]), "", ip[, "Priority"]),
                 Repository = ifelse(is.na(ip[, "Repository"]), "", ip[, "Repository"]),
                 RemoteSha = ifelse(is.na(ip[, "RemoteSha"]), "", ip[, "RemoteSha"]))
df <- df[order(df$LibPath, df$Package), ]
write.table(df, stdout(), sep = "\\t", quote = FALSE, row.names = FALSE)"""])
r_now = set(r.stdout.splitlines())
r_lock = set(open(f"{L}/atlas-2026.10c.R.tsv").read().splitlines())
check("R installed.packages = 10c", r.returncode == 0 and r_now == r_lock,
      f"only now {sorted(r_now - r_lock)[:8]}, only 10c {sorted(r_lock - r_now)[:8]}")

# 3. dpkg
d_now = set(run(["dpkg-query", "-W", "-f=${Package}\t${Version}\t${Architecture}\n"]).stdout.splitlines())
d_lock = set(open(f"{L}/atlas-2026.10c.dpkg.tsv").read().splitlines())
check("dpkg packages = 10c", d_now == d_lock, f"only now {sorted(d_now - d_lock)[:8]}, only 10c {sorted(d_lock - d_now)[:8]}")

# 4. /opt
opt_now = {"./" + e for e in os.listdir("/opt")}
opt_lock = set(open(f"{L}/atlas-2026.10c.opt.txt").read().split())
check("/opt = 10c + atlas-locks, cellsnp", opt_now - opt_lock == {"./atlas-locks", "./cellsnp"} and not opt_lock - opt_now,
      f"added {sorted(opt_now - opt_lock)}, missing {sorted(opt_lock - opt_now)}")

# 5. system versions
sys_now = [run(["python3", "--version"]).stdout.strip(), run(["R", "--version"]).stdout.splitlines()[0],
           [x for x in open("/etc/os-release").read().splitlines() if x.startswith("PRETTY_NAME")][0]]
sys_lock = open(f"{L}/atlas-2026.10c.system.txt").read().splitlines()
check("Python, R, Ubuntu versions = 10c", [s.strip() for s in sys_now] == [s.strip() for s in sys_lock], f"{sys_now}")

# 6. cellsnp-lite
v = run(["cellsnp-lite", "--version"])
check("cellsnp-lite 1.2.3 (htslib 1.24)", v.stdout.strip() == "cellsnp-lite 1.2.3 (htslib 1.24)", v.stdout.strip() + v.stderr.strip())
lp = run(["ldd", "/opt/cellsnp/bin/cellsnp-lite"])
ldd = lp.stdout.splitlines()
own = [x for x in ldd if "=>" in x and not any(s in x for s in ("/lib/x86_64-linux-gnu/", "/lib64/", "linux-vdso"))]
check("cellsnp-lite libraries resolved, own ones in /opt/cellsnp",
      lp.returncode == 0 and bool(ldd) and not any("not found" in x for x in ldd) and all("/opt/cellsnp/" in x for x in own),
      "; ".join(x.strip() for x in own) or lp.stderr.strip())

# 7. vireo on simulated cellsnp-lite output: 2 donors, 300 SNPs, 200 cells
try:
    import numpy as np
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix
    rng = np.random.default_rng(1)
    n_snp, n_cell = 300, 200
    geno = rng.integers(0, 3, size=(n_snp, 2))
    donor = np.repeat([0, 1], n_cell // 2)
    af = np.array([0.02, 0.5, 0.98])[geno[:, donor]]
    dp = rng.poisson(2.0, size=(n_snp, n_cell))
    adm = rng.binomial(dp, af)
    with tempfile.TemporaryDirectory() as t:
        c = os.path.join(t, "cellsnp")
        os.mkdir(c)
        with gzip.open(os.path.join(c, "cellSNP.base.vcf.gz"), "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for k in range(n_snp):
                f.write(f"1\t{1000 + 10 * k}\t.\tA\tG\t.\tPASS\t.\n")
        open(os.path.join(c, "cellSNP.samples.tsv"), "w").write("\n".join(f"cell{k}" for k in range(n_cell)) + "\n")
        mmwrite(os.path.join(c, "cellSNP.tag.AD.mtx"), csr_matrix(adm))
        mmwrite(os.path.join(c, "cellSNP.tag.DP.mtx"), csr_matrix(dp))
        o = os.path.join(t, "vireo")
        vr = run(["vireo", "-c", c, "-N", "2", "-o", o, "--randSeed", "1"])
        ids = [line.split("\t") for line in open(os.path.join(o, "donor_ids.tsv")).read().splitlines()[1:]]
        call = {x[0]: x[1] for x in ids}
        purity = 0
        for g in (0, 1):
            calls = [call[f"cell{k}"] for k in range(n_cell) if donor[k] == g]
            purity += max(calls.count(x) for x in set(calls))
        check("vireo separates 2 simulated donors", vr.returncode == 0 and purity / n_cell >= 0.95,
              f"purity {purity / n_cell:.3f}; {vr.stderr.strip()[-200:]}")
except Exception as e:  # noqa: BLE001
    check("vireo separates 2 simulated donors", False, repr(e))

# 8. home leaks
import site  # noqa: E402
r = run(["Rscript", "-e", "cat(Sys.getenv('R_ENVIRON_USER'), Sys.getenv('R_PROFILE_USER'), Sys.getenv('R_LIBS_USER'), .libPaths(), sep = '|')"])
check("R ignores ~/.Renviron, ~/.Rprofile", r.stdout.startswith("/dev/null|/dev/null|/opt/r-user-lib-unused|")
      and all(p.startswith("/usr/local/lib/R/") for p in r.stdout.split("|")[3:]), r.stdout)
check("Python user site disabled", not site.ENABLE_USER_SITE and os.environ.get("PYTHONNOUSERSITE") == "1")

# 9. earlier load tests
for cmd in (["Rscript", "/opt/atlas-tests/atlas_load_test.R"], [sys.executable, "/opt/atlas-tests/atlas_load_test.py"],
            [sys.executable, "/opt/atlas-tests/atlas_b_load_test.py"], [sys.executable, "/opt/atlas-tests/atlas_c_load_test.py"]):
    p = run(cmd)
    last = (p.stdout.strip().splitlines() or [""])[-1]
    check(os.path.basename(cmd[-1]), p.returncode == 0, last if p.returncode == 0 else last + " | " + p.stderr.strip()[-300:])

if fail:
    print("FAILED:", ", ".join(fail))
    sys.exit(1)
print("ALL ATLAS-FINAL LOAD TESTS PASSED")
