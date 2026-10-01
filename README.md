# Autoimmune Atlas

Analysis code for a single-cell RNA-seq comparison of peripheral blood
mononuclear cells (PBMCs) from systemic lupus erythematosus (SLE), Sjögren's
disease (SjD) and healthy donors. The scripts keep the earlier abbreviation SjS
(`SjS`, `sjs`) in variable, column, folder and file names.

**This repository holds the analysis code only.** It contains no data, no
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

## Gene symbols and gene exclusion

The datasets name some genes differently (for example RIGI and DDX58), so before
any step that combines datasets every dataset's genes are mapped onto the
approved symbols of the HGNC complete set of 28 September 2026: by Ensembl gene
ID where the dataset deposits one, otherwise by approved, previous or alias
symbol when that symbol names exactly one approved gene (`harmonise_genes()` in
`00_utils.R`, applied in 03 and 04a). The Hallmark gene sets, the exclusion
lists and the Toro-Domínguez signature use the same mapping.

Mitochondrial genes, HGNC ribosomal protein genes, haemoglobin genes,
Y-chromosome genes, XIST, TSIX, the erythroid markers SLC4A1, ALAS2, AHSP and
CA1, and the dissociation-induced genes of van den Brink et al. 2017 with HSPA6
are excluded from differential expression, pathway enrichment and the classifier
(`is_excluded_gene()` in `00_config.R`).

## Reference files

`scripts/reference/` holds the fixed inputs, each with its source:

- `hgnc_complete_set_2026-09-28.txt.gz`: the HGNC complete set used for gene
  symbol harmonisation (provenance in the accompanying `.README`).
- `dissociation_genes_vandenbrink2017.tsv`: the dissociation-induced genes of
  van den Brink et al., *Nat. Methods* 2017
  ([doi:10.1038/nmeth.4437](https://doi.org/10.1038/nmeth.4437)), mapped to
  human orthologs.
- `toro2014/`: the shared SLE, rheumatoid arthritis and SjD signature, its GO
  terms and its gain genes from Additional file 1 of Toro-Domínguez et al.,
  *Arthritis Res. Ther.* 2014
  ([doi:10.1186/s13075-014-0489-x](https://doi.org/10.1186/s13075-014-0489-x)),
  read by `06b_toro_comparison.R`.

## Example: SjD robustness within GSE253568 (07b)

Run from the repository root once 03 and 04a have produced the pseudobulk counts
and the atlas. The settings are environment variables:

```sh
ATLAS_ROOT="$PWD" ROBUSTNESS_DATASET=GSE253568 MIN_DONORS=8 MIN_CELLS=100 \
MIN_CELLS_CLASS=10 MIN_CELLS_PER_DONOR=0 N_CORES=8 \
Rscript --vanilla scripts/07b_sjd_single_dataset.R
```
