# Autoimmune Atlas

Analysis code for a single-cell RNA-seq comparison of peripheral blood
mononuclear cells (PBMCs) from systemic lupus erythematosus (SLE), Sjögren's
syndrome (SjS) and healthy donors.

**This repository holds the analysis code only.** It contains no data, no
results and no intermediate objects. All four datasets are public and are
named below, but two sets of hand-prepared inputs are **not included**: the
per-dataset sample sheets and the extracts from the Toro-Domínguez et al.
(2014) supplementary files. **The pipeline will therefore not run as-is from a
clone.** See [What you need to supply](#what-you-need-to-supply).

## Datasets

All data are public. Download them from the listed source.

| Source | Accession | Disease | Donors | Publication |
|---|---|---|---|---|
| CELLxGENE | collection `436154da-bcf1-4130-9c8b-120ff9a888f2` | SLE | 261 | Perez et al., *Science* 2022, [doi:10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970) |
| GEO | [GSE162577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162577) | SLE | 3 | Deng et al., *EBioMedicine* 2021, [doi:10.1016/j.ebiom.2021.103477](https://doi.org/10.1016/j.ebiom.2021.103477) |
| GEO | [GSE157278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157278) | SjS | 10 | Hong et al., *Front. Immunol.* 2021, [doi:10.3389/fimmu.2020.594658](https://doi.org/10.3389/fimmu.2020.594658) |
| GEO | [GSE253568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253568) | SjS | 17 | McDermott et al., *J. Autoimmun.* 2025, [doi:10.1016/j.jaut.2025.103419](https://doi.org/10.1016/j.jaut.2025.103419) |

The scripts expect each dataset under `data/<disease>/<accession>/`, with the
downloaded files in `raw/` (10x MTX, HDF5 or H5AD) and a sample sheet in `meta/`.

## What you need to supply

1. **Sample sheets (not included).** `scripts/01_qc_preproc.R` reads
   `data/<disease>/<accession>/meta/samples.tsv`, one row per sample, with the
   columns `dataset_id`, `batch_id`, `sample_id`, `gsm_id`, `donor_id`,
   `condition`, `age_group`, `sex` and `sex_source`. These map each library to
   its donor and condition. They were prepared by hand from the repository
   metadata of each dataset; where a dataset records no sex, `01_qc_preproc.R`
   infers it from XIST and Y-chromosome gene expression and fills it in.
2. **Toro-Domínguez et al. (2014) supplementary extracts (not included).**
   `scripts/06b_toro_comparison.R` compares the results with the shared bulk
   PBMC signature of Toro-Domínguez et al., *Arthritis Res. Ther.* 16 (2014) 489,
   [doi:10.1186/s13075-014-0489-x](https://doi.org/10.1186/s13075-014-0489-x).
   It reads three CSV files extracted from that article's additional files, all
   in `results/manuscript1_figures/toro_comparison/`:
   `toro2014_sheet1_shared_signature.csv`, `toro2014_sheet4_go_enrichment.csv`
   and `toro2014_addfile3_gain_loss.csv`.
3. **MSigDB Hallmark gene sets.** `scripts/07_pathways.R` and
   `scripts/08_donor_classifiers.R` read the MSigDB 2025.1 human database in
   msigdbr's format, `msigdb.2025.1.Hs.rds`. Set `MSIGDB_RDS` to its path; the
   default is msigdbr's user data directory. The scripts use only the 50 Hallmark
   sets.
4. **SingleR reference (fetched automatically).** On first run
   `scripts/02_annotate.R` downloads the Monaco immune reference with
   `celldex::MonacoImmuneData()` and caches it under `data/references/`. This
   step needs internet access.

## Software

The analyses were run in a Singularity container providing R 4.3.3,
Bioconductor 3.18 and Python 3.10. The image definition is at
<https://github.com/Davidn1095/research-env>.

Main packages and versions: Seurat 5.4, SingleR 2.4, celldex 1.12,
scDblFinder 1.16, DESeq2 1.42, limma 3.58, speckle 1.2 (propeller), fgsea 1.28,
glmnet 5.0, ranger 0.18, xgboost 3.1 and scvi-tools 1.3. Training the scVI model
(`04b_train_scvi.py`) needs a GPU.

## Configuration

Settings are read from environment variables. `scripts/env_setup.sh` is a
template that sets the paths.

| Variable | Used by | Meaning |
|---|---|---|
| `ATLAS_ROOT` | all | Repository root; `data/` and `results/` sit beneath it. Defaults to the working directory. |
| `ACC_DIR` | 01-03 | Dataset directory, `data/<disease>/<accession>`. Defaults to the working directory; accession and disease are taken from its path. |
| `DENOISE_METHOD` | 08, 09 | Required. `limma_modulescore` is the configuration used for the classifier. |
| `SCVI_EXTRACT` | 04a | Set to `1` to write the inputs for `04b_train_scvi.py`. |
| `MSIGDB_RDS` | 07, 08 | Path of the MSigDB 2025.1 Hs cache. |
| `REF_RDS` | 02 | Optional path of an existing SingleR reference. |
| `ATLAS_EXTRA_R_LIB` | 05, 08, 09 | Optional extra R library searched first. |

Analysis thresholds (donor and cell floors, significance cut-offs) are defined
in `scripts/00_config.R`, and most can be overridden by environment variables.

## Pipeline

Run the per-dataset scripts once for each of the four datasets, then the
atlas-level scripts in order.

```
Per dataset (ACC_DIR = data/<disease>/<accession>)
  01_qc_preproc.R              QC (adaptive IQR thresholds, scDblFinder), normalisation
  02_annotate.R                SingleR cell-type annotation (Monaco reference)
  03_pseudobulk.R              Donor-level pseudobulk counts per cell type

Atlas level
  04a_integrate.R              Merge the four datasets (SCVI_EXTRACT=1 writes scVI inputs)
  04b_train_scvi.py            scVI latent space, used for the UMAP visualisation only
  05_composition.R             Cell-type composition (propeller)
  06_differential_analysis.R   Pseudobulk differential expression (DESeq2, ~ dataset + condition)
  07_pathways.R                Hallmark gene set enrichment (fgsea on the DESeq2 Wald statistic)
  06b_toro_comparison.R        Comparison with the Toro-Domínguez 2014 signature (after 07)

Donor-level classification
  08_donor_classifiers.R       Pathway-by-cell-type features; elastic net, random forest
                               and XGBoost under repeated stratified cross-validation
  09_shap_importance.R         Exact linear SHAP attribution for the elastic net

Figures and tables
  10_manuscript_figures.R      Main figures
  11_supplementary_tables.R    Supplementary tables
```

Shared settings and helper functions are in `scripts/00_config.R` and
`scripts/00_utils.R`.
