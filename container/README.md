# Atlas container (atlas-final-2026.10)

The software environment of the Python pipeline of the Genes & Immunity revision: R 4.6.1, Bioconductor 3.23,
Python 3.12, cellsnp-lite 1.2.3 and vireoSNP 0.5.9, with every version pinned.

| | |
|---|---|
| Image | `ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10` |
| Image digest | `sha256:09bbfeca81f191614ebc0a644847aea0a27f9b3f59ec9992bbe18a538359d213` |
| SIF (Apptainer) | `oras://ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10-sif` |
| SIF sha256 | `c7f26ae98de80940c7e4d4565c0a4e60d54344172aa0b877e631d9e7ca4cf89a` |

The image and the SIF were built on 2026-10-28 in github.com/Davidn1095/research-env and copied here unchanged
(`.github/workflows/copy-atlas-final.yml`), so the digest and the sha256 are those of the original build.

Contents:

- `Dockerfile.atlas-final`, `Dockerfile.atlas-final.dockerignore`: the recipe. Base image by digest, Ubuntu snapshot
  2026-10-23, Posit CRAN snapshot 2026-10-22, Bioconductor 3.23, a full pip freeze, cellsnp-lite's conda build by URL
  and MD5. Only comments and the `image.source` label differ from the recipe the published image was built from.
- `build/atlas_final_install.R`: R installs; every package must end at its locked version.
- `locks/`: package lists (pip, R, dpkg, `/opt`, system versions) and the cellsnp-lite conda lock.
- `tests/`: `atlas_final_load_test.py` (versions equal to the locks, vireo on simulated donors, home isolation) and
  the earlier tests it runs.

Use:

```
apptainer pull atlas-final-2026.10.sif oras://ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10-sif
sha256sum atlas-final-2026.10.sif    # c7f26ae9...
```

To rebuild from the recipe, run `build-atlas-final.yml` and then `build-atlas-final-sif.yml` (by hand,
workflow_dispatch). A rebuild gets a new digest; the published tag stays the record of the analysis.
