#!/bin/bash
# Genotype demultiplexing of one Cell Ranger library: cellsnp-lite on its BAM, vireo for each N in N_LIST, both in the
# atlas container. Input: $DATA_ROOT/<disease>/<DATASET>/cellranger/<LIB>/ (cellranger.sh with CREATE_BAM=true).
# Output: $PIPELINE_OUT/genotype/<DATASET>/<LIB>/ (read by 00_gse253568_genotype_labels.py).
# Usage: bash pipeline_py/genotype.sh <LIB> <DATASET> <DISEASE>, or sbatch it with your own options
#        (e.g. sbatch -c 16 --mem=64G -t 24:00:00 --signal=B:USR1@300 pipeline_py/genotype.sh capture1_well1 GSE253568 SjD)
#        N_LIST (env) sets the vireo donor numbers, default "1 2 3 4 5".
# Environment (defaults in brackets; <repo> is the folder above pipeline_py/): ATLAS_ROOT [<repo>], DATA_ROOT [<repo>/data],
#   PIPELINE_OUT [<repo>/out], ATLAS_SIF [<repo>/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif],
#   SNP_PANEL (required: 1000 Genomes phase 3 SNPs, allele frequency >= 5%, GRCh38, VCF; see README),
#   APPTAINER_MODULE (environment module loaded first, if set)
set -euo pipefail
LIB=$1; DS=$2; DIS=$3; N_LIST=${N_LIST:-1 2 3 4 5}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${ATLAS_ROOT:-$(dirname "$HERE")}
CR=${DATA_ROOT:-$ROOT/data}/$DIS/$DS/cellranger/$LIB; GT=${PIPELINE_OUT:-$ROOT/out}/genotype/$DS/$LIB; mkdir -p $GT
SIF=${ATLAS_SIF:-$ROOT/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif}
PANEL=${SNP_PANEL:?set SNP_PANEL to the 1000 Genomes SNP VCF (README)}
if [ -n "${APPTAINER_MODULE:-}" ] && command -v module >/dev/null 2>&1; then module load "$APPTAINER_MODULE"; fi
# own temporary folder, cleanup on exit/USR1/TERM/INT (job_localscratch.sh)
source "$HERE/job_localscratch.sh"
export APPTAINER_BIND=$CR,$GT,$(dirname "$PANEL"),$JOBTMP
export APPTAINERENV_PYTHONNOUSERSITE=1 APPTAINERENV_R_ENVIRON_USER=/dev/null APPTAINERENV_R_PROFILE_USER=/dev/null
img() { apptainer exec "$SIF" "$@"; }
echo "node $(hostname) | cpu $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs) | start $(date)"; img cellsnp-lite --version; img python -c "import vireoSNP; print('vireoSNP', vireoSNP.__version__)"
# barcode list from the filtered .h5
img python -c "import h5py, gzip; bc = h5py.File('$CR/filtered_feature_bc_matrix.h5')['matrix/barcodes'][:].astype(str); gzip.open('$GT/barcodes.tsv.gz', 'wt').write(chr(10).join(bc) + chr(10))"
t0=$(date +%s)
NP=${SLURM_CPUS_PER_TASK:-16}
ls_run apptainer exec "$SIF" cellsnp-lite -s $CR/possorted_genome_bam.bam -b $GT/barcodes.tsv.gz \
  -R $PANEL -O $GT/cellsnp -p $NP --minMAF 0.1 --minCOUNT 20 --cellTAG CB --UMItag UB --gzip
echo "cellsnp-lite $(( ($(date +%s)-t0)/60 )) min; SNPs kept: $(zcat $GT/cellsnp/cellSNP.base.vcf.gz | grep -vc '^#')"
for N in $N_LIST; do
  t=$(date +%s)
  ls_run apptainer exec "$SIF" vireo -c $GT/cellsnp -N $N -o $GT/vireo_N$N --randSeed 1 -p $NP > $GT/vireo_N$N.log 2>&1 || echo "vireo N=$N exit $?"
  echo "vireo N=$N $(( ($(date +%s)-t)/60 )) min"; grep -iE 'elbo|lower bound|logLik|final' $GT/vireo_N$N.log | tail -3
done
echo "total $(( ($(date +%s)-t0)/60 )) min; end $(date)"; echo GENOTYPE_DONE
