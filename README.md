# comparative-scrnaseq-pipeline

Code for the analyses in "A computational framework for cell-type-resolved comparative transcriptomics of systemic lupus erythematosus and Sjögren's disease" (Núñez-Nepomuceno et al., in preparation).

Public single-cell RNA-seq datasets of peripheral blood mononuclear cells from systemic lupus erythematosus (SLE) and Sjögren's disease (SjD), each with its own healthy controls, are integrated into one atlas. Cell-type composition, gene expression and pathway enrichment are compared between each disease and its controls and between the two diseases.

## Data

| Dataset | Disease | Source |
|---|---|---|
| GSE174188 | SLE | GEO, CELLxGENE |
| GSE135779 | SLE | GEO |
| GSE253568 | SjD | GEO |
| GSE316482 | SjD | GEO |
| HRA003613 | SjD | Zenodo (10884425) |

## Software

Analyses were run in the container image `ghcr.io/davidn1095/comparative-scrnaseq-pipeline:atlas-final-2026.10`, whose recipe and package versions are in `container/`. Main packages: scanpy 1.12.4, anndata 0.13.4, scvi-tools 1.5.1, leidenalg 0.12.0, CellTypist 1.7.1, decoupler 2.2.0, pertpy 1.4.0, PyDESeq2 0.5.4, celda 1.28.0 (DecontX), scDblFinder 1.26.9, fgsea 1.38.0, cellsnp-lite 1.2.3 and vireoSNP 0.5.9. Cell Ranger 10.0.0 was run outside the container.

The scripts in `pipeline_py/` are numbered in the order they run, from data input (00: `00_genotype.py`, `00_inputs.py`) to pathway enrichment (06), followed by the figure scripts (`07_figure_data.py`, `08_figures.R`), which write every figure and table. `cellranger.sh` runs Cell Ranger for each 10x library of the datasets reprocessed from reads.

## Citation

Núñez-Nepomuceno D, Liñares-Blanco J, Alarcón-Riquelme ME, Carmona-Sáez P. A computational framework for cell-type-resolved comparative transcriptomics of systemic lupus erythematosus and Sjögren's disease. In preparation.
