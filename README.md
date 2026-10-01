# Autoimmune Atlas

Analysis code for a single-cell RNA-seq comparison of peripheral blood
mononuclear cells (PBMCs) from systemic lupus erythematosus (SLE), Sjögren's
disease (SjD) and healthy donors.

This repository holds the analysis code only. It contains no data, no
results and no intermediate objects. All four datasets are public and are
named below.

## Datasets

All data are public. Download them from the listed source.

| Source | Accession | Disease | Donors | Publication |
|---|---|---|---|---|
| CELLxGENE | collection `436154da-bcf1-4130-9c8b-120ff9a888f2` | SLE | 261 | Perez et al., *Science* 2022, [doi:10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970) |
| GEO | [GSE162577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162577) | SLE | 3 | Deng et al., *EBioMedicine* 2021, [doi:10.1016/j.ebiom.2021.103477](https://doi.org/10.1016/j.ebiom.2021.103477) |
| GEO | [GSE157278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157278) | SjD | 10 | Hong et al., *Front. Immunol.* 2021, [doi:10.3389/fimmu.2020.594658](https://doi.org/10.3389/fimmu.2020.594658) |
| GEO | [GSE253568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253568) | SjD | 17 | McDermott et al., *J. Autoimmun.* 2025, [doi:10.1016/j.jaut.2025.103419](https://doi.org/10.1016/j.jaut.2025.103419) |

## Pipeline

The scripts in `scripts/` run in numerical order; `00_config.R` and `00_utils.R`
hold the shared settings and helpers.

| Script | Step |
|---|---|
| `01_qc_preproc.R` | Per-dataset quality control, doublet removal, normalisation |
| `02_annotate.R` | Cell-type annotation with SingleR and the Monaco immune reference |
| `03_pseudobulk.R` | Pseudobulk counts per donor, condition and cell type |
| `04a_integrate.R` | Merge of the four datasets into the atlas |
| `04b_train_scvi.py` | scVI latent space and UMAP, used for visualisation only |
| `05_composition.R` | Cell-type composition (propeller) |
| `06_differential_analysis.R` | Pseudobulk differential expression (DESeq2) per cell type |
| `06b_toro_comparison.R` | Recovery of the shared bulk signature of Toro-Domínguez et al. 2014 |
| `07_pathways.R` | Hallmark pathway enrichment (fgsea) |
| `07b_sjd_single_dataset.R` | SjD against healthy controls within GSE253568 alone, a robustness check |
| `08_donor_classifiers.R` | Donor-level classifiers on Hallmark pathway scores per cell type |
| `09_shap_importance.R` | SHAP attribution of the elastic-net classifier |
| `10_manuscript_figures.R` | Manuscript figures |
| `11_supplementary_tables.R` | Supplementary tables |
