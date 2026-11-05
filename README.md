# comparative-scrnaseq-pipeline

Code for a cell-type-resolved comparison of systemic lupus erythematosus (SLE) and Sjögren's disease (SjD). Public
single-cell RNA sequencing datasets of peripheral blood mononuclear cells, each with its own healthy controls, are
reprocessed from raw reads where available (Cell Ranger, exonic counts) or from the deposited count matrices, quality
controlled per sample, integrated with scVI and annotated from marker genes into one atlas. Cell-type composition
(scCODA), pseudobulk gene expression (PyDESeq2) and pathway enrichment (fgsea, MSigDB Hallmark) are compared between
each disease and its healthy controls, and SLE and SjD are compared indirectly through their controls, with an
equivalence test. The pipeline follows the single-cell best-practices book (Heumos et al.); each script names the
notebook it ports, and every departure is marked `# CHANGED:`. The repository holds code only, no data or results.

## Datasets

| Accession | Disease | Source | Publication |
|---|---|---|---|
| GSE174188 | SLE | CELLxGENE collection 436154da-bcf1-4130-9c8b-120ff9a888f2 (cell labels); FASTQ files from the Human Cell Atlas Data Portal | Perez et al., Science 2022 |
| GSE135779 (adult arm) | SLE | GEO | Nehar-Belaid et al., Nat Immunol 2020 |
| GSE253568 | SjD | GEO (hashtag labels); FASTQ files from SRA | McDermott et al., J Autoimmun 2025 |
| GSE316482 | SjD | GEO | no associated publication |
| HRA003613 (PBMC arm) | SjD | count matrix on Zenodo (10884425); reads on GSA-Human (controlled access) | Xiang et al., iScience 2023 |

Downloads go under `data/<disease>/<accession>/` (`DATA_ROOT`, default `<repo>/data`), as `00_deposited_input.py` reads them:

| Folder | Contents |
|---|---|
| `data/SLE/GSE174188/raw/GSE174188.h5ad` | the CELLxGENE h5ad of the collection above (barcodes, samples, donors) |
| `data/SLE/GSE174188/cellranger/<library>/` | Cell Ranger outputs of each library (`cellranger.sh`, FASTQ_DIR mode) |
| `data/SLE/GSE135779/raw/extracted/` | the GEO per-sample files of the adult arm (`<GSM>_<sample>_matrix.mtx.gz`, `_barcodes.tsv.gz`, `_genes.tsv.gz`) |
| `data/SjD/GSE253568/raw/extracted/` | the GEO per-sample hashtag barcodes (`<GSM>_<sample>_barcodes.tsv.gz`) |
| `data/SjD/GSE253568/cellranger/<capture>/` | Cell Ranger outputs per capture (`cellranger.sh`, SRA mode; capture1 is GEM well 1, `capture1_well1`, run with `CREATE_BAM=true` for `genotype.sh`) |
| `data/SjD/GSE316482/raw/extracted/` | the GEO per-sample files (`<GSM>_<sample>_scRNAseq_all_cells_matrix.mtx.gz`, `_barcodes.tsv.gz`, `_features.tsv.gz`) |
| `data/SjD/HRA003613/raw/SG_PBMC_immune_inte_seurat_rpca_pc20_neighbor50.h5ad` | the Zenodo h5ad (integer raw counts in `.raw`) |

`pipeline_py/inputs/` holds the sample tables built from GEO (`GSE135779_adult_samples.tsv`, `GSE253568_samples.tsv`,
`GSE316482_samples.tsv`), the manual cluster annotation of the published atlas (`manual_annotation.tsv`) and the
numbering of the supplementary items (`supplementary_numbering.tsv`).

## Software

All steps after read alignment run in one container:

- Image: `ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10`
  (digest `sha256:09bbfeca81f191614ebc0a644847aea0a27f9b3f59ec9992bbe18a538359d213`)
- Apptainer/Singularity: `apptainer pull oras://ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10-sif`
  (SIF sha256 `c7f26ae98de80940c7e4d4565c0a4e60d54344172aa0b877e631d9e7ca4cf89a`)

The recipe, package locks and load tests are in `container/`; the build workflows are in `.github/workflows/`.
Cell Ranger 10.0.0 with the GRCh38-2024-A reference, the SRA Toolkit and pigz run outside the container (`cellranger.sh`).
Fig 1 needs pdflatex (TeX Live 2024) and the PNG previews Ghostscript, both on the host.

Reference files in `pipeline_py/reference/`: the MSigDB Hallmark 2025.1 gene sets are included (CC BY 4.0). Two files are
not redistributed and must be downloaded first, as their READMEs describe: the CellTypist curated-marker table
(`celltypist_v2_Basic_celltype_information.xlsx`) and the GRCh38-2020-A feature list (`GRCh38-2020-A_features.tsv.gz`).
Genotype demultiplexing (`genotype.sh`) needs the 1000 Genomes phase 3 SNPs with allele frequency of at least 5% on
GRCh38, as used by cellsnp-lite (`SNP_PANEL`).

## Running

Environment (defaults in brackets; `<repo>` is this repository):

| Variable | Meaning |
|---|---|
| `PIPELINE_OUT` | output folder [`<repo>/out`] |
| `DATA_ROOT` | raw downloads [`<repo>/data`] |
| `REF_DIR` | Cell Ranger GRCh38-2024-A reference [`<repo>/refdata/refdata-gex-GRCh38-2024-A`] |
| `ATLAS_SIF` | the container image [`<repo>/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif`] |
| `ATLAS_RUN` | folder of the integration run under `$PIPELINE_OUT/py/atlas/` [`02_integration`] |
| `PIPELINE_BIND` | extra folders to bind into the container, comma-separated |
| `APPTAINER_MODULE`, `TEX_MODULE`, `CELLRANGER_MODULE`, `SRA_MODULE` | environment modules to load first, if your system uses them |
| `PIPELINE_LOCAL_TMP` | folder for per-job temporary files [`$TMPDIR`, else `/tmp`] |

Each step runs with `bash pipeline_py/run.sh <script> [args]`, or on Slurm with
`bash pipeline_py/run.sh --sbatch "<sbatch options>" <script> [args]`, which submits from a read-only snapshot of
`pipeline_py/`. `cellranger.sh` and `genotype.sh` run on the host (`bash` or `sbatch`).

| Order | Script | Step |
|---|---|---|
| 0 | `cellranger.sh <dataset> <disease> <library> [SRA runs]` | Cell Ranger per library (GSE174188, GSE253568) |
| 0 | `genotype.sh capture1_well1 GSE253568 SjD`, then `00_gse253568_genotype_labels.py` | GSE253568 hashtag labels checked by genotype (cellsnp-lite, vireo) |
| 0 | `00_deposited_input.py <dataset>` | per-library inputs of each dataset |
| 1 | `01_quality_control.py <dataset> [--stage qc / library / merge]` | quality control, ambient RNA (DecontX), doublets (scDblFinder) |
| 2 | `02_integration.py` | merge, gene filter, normalisation, highly variable genes, scVI, Leiden |
| 3 | `03_annotation.py`, then `03_annotation.py --apply` | marker and CellTypist tables; the manual labels of `inputs/manual_annotation.tsv` |
| 4 | `04_differential_expression.py`; `04b_percohort_de.py` | pseudobulk differential expression; per-cohort analysis and meta-analysis |
| 5 | `05_composition.py`; `05b_reference_sensitivity.py <base> <reference>`, then `05b_reference_sensitivity.py --summary` | scCODA; refits with five further reference cell types |
| 6 | `06_pathways.py` | fgsea, MSigDB Hallmark 2025.1 |
| figures | `make_figures.sh`; `make_supplementary.sh` | Figs 1-4 and Table 1; supplementary figures and tables |

Dataset names in the code follow `pipeline_py/datasets.tsv`; `Perez` is GSE174188 and `GSE135779_adult` the adult arm of
GSE135779.

## Citation

Núñez-Nepomuceno D, Liñares-Blanco J, Alarcón-Riquelme ME, Carmona-Sáez P. A computational framework for
cell-type-resolved comparative transcriptomics of systemic lupus erythematosus and Sjögren's disease. In preparation.
