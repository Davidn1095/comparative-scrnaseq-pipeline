#!/bin/bash
# Cell Ranger 10.0.0 count for one 10x library of a public-FASTQ dataset (GSE174188, GSE253568), against the
# GRCh38-2024-A reference with default settings. Runs on the host (Cell Ranger is not in the container).
# Usage:
#   FASTQ_DIR=<dir> SAMPLES=<name>[,<name>] bash pipeline_py/cellranger.sh <DATASET> <DISEASE> <LIB>
#       10x-named FASTQs already in one folder (GSE174188: the FASTQs of each library, named by its sample names)
#   bash pipeline_py/cellranger.sh <DATASET> <DISEASE> <LIB> <SRR> [<SRR> ...]
#       SRA runs of one library (GSE253568); each run is converted with fasterq-dump --include-technical, and the reads
#       are assigned by length: 26-28 bp = R1 (barcode + UMI), the longest = R2 (cDNA); index reads are dropped
#   Either can be submitted with sbatch and your own options, e.g. sbatch -c 32 --mem=128G -t 3-00:00:00 --tmp=700G
#   --signal=B:USR1@300 pipeline_py/cellranger.sh ...
# Environment (defaults in brackets; <repo> is the folder above pipeline_py/):
#   ATLAS_ROOT [<repo>], DATA_ROOT [<repo>/data], REF_DIR [<repo>/refdata/refdata-gex-GRCh38-2024-A]
#   SRA_DIR [$DATA_ROOT/<DISEASE>/<DATASET>/sra]: folder of the prefetched runs, <SRR>/<SRR>.sra
#   CREATE_BAM [false]: true keeps the BAM (needed by genotype.sh, e.g. GSE253568 capture1_well1)
#   CELLRANGER_MODULE, SRA_MODULE: environment modules to load first, if set; otherwise cellranger, fasterq-dump and
#   pigz must be on PATH
#   PIPELINE_LOCAL_TMP [$TMPDIR, else /tmp]: work folder for the FASTQ conversion and Cell Ranger (needs ~250 GB free)
# Output: $DATA_ROOT/<DISEASE>/<DATASET>/cellranger/<LIB>/ with the filtered and raw .h5 matrices, molecule_info.h5,
#   metrics_summary.csv and web_summary.html (and the BAM when CREATE_BAM=true), read by 00_deposited_input.py.
set -euo pipefail
DS=$1; DIS=$2; LIB=$3; shift 3; RUNS=("$@")
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${ATLAS_ROOT:-$(dirname "$HERE")}
DATA=${DATA_ROOT:-$ROOT/data}
REF=${REF_DIR:-$ROOT/refdata/refdata-gex-GRCh38-2024-A}
SRA=${SRA_DIR:-$DATA/$DIS/$DS/sra}
OUT=$DATA/$DIS/$DS/cellranger/$LIB
CREATE_BAM=${CREATE_BAM:-false}; FASTQ_DIR=${FASTQ_DIR:-}; SAMPLES=${SAMPLES:-$LIB}
source "$HERE/job_localscratch.sh"   # own work folder, deleted on exit, USR1, TERM or INT
FQ=$JOBTMP/fastq; W=$JOBTMP/work; mkdir -p "$FQ" "$W"
if command -v module >/dev/null 2>&1; then
  [ -n "${CELLRANGER_MODULE:-}" ] && module load "$CELLRANGER_MODULE"
  [ -n "${SRA_MODULE:-}" ] && module load "$SRA_MODULE"
fi
cellranger --version
NCORE=${SLURM_CPUS_PER_TASK:-$(nproc)}; NMEM=$(( ${SLURM_MEM_PER_NODE:-131072} / 1024 - 8 ))   # 8 GB left for the system
echo "$DS $LIB on $(hostname) | cellranger --localcores $NCORE --localmem $NMEM | start $(date)"

if [ -n "$FASTQ_DIR" ]; then
  ls_run cp -p "$FASTQ_DIR"/* "$FQ"/
else
  [ ${#RUNS[@]} -gt 0 ] || { echo "give FASTQ_DIR or at least one SRA run" >&2; exit 1; }
  i=0
  for r in "${RUNS[@]}"; do
    i=$((i+1)); L=$(printf 'L%03d' $i)
    ls_run fasterq-dump --split-files --include-technical -e "$NCORE" -t "$JOBTMP/tmp_$r" -O "$FQ" "$SRA/$r/$r.sra"
    R1=""; R2=""; R2LEN=0
    for f in "$FQ"/${r}_*.fastq; do len=$(sed -n 2p "$f" | tr -d '\n' | wc -c); echo "  $(basename "$f"): read length $len"
      if [ "$len" -ge 26 ] && [ "$len" -le 28 ] && [ -z "$R1" ]; then R1=$f; elif [ "$len" -gt $R2LEN ]; then R2=$f; R2LEN=$len; fi; done
    [ -n "$R1" ] && [ -n "$R2" ] && [ "$R1" != "$R2" ] || { echo "cannot identify R1/R2 for $r" >&2; exit 1; }
    for f in "$FQ"/${r}_*.fastq; do [ "$f" = "$R1" ] || [ "$f" = "$R2" ] || rm -f "$f"; done
    ls_run pigz -p "$NCORE" "$R1" "$R2"
    mv "$R1.gz" "$FQ/${LIB}_S1_${L}_R1_001.fastq.gz"; mv "$R2.gz" "$FQ/${LIB}_S1_${L}_R2_001.fastq.gz"
    rm -rf "$JOBTMP/tmp_$r"
  done
fi

cd "$W"
ls_run cellranger count --id="$LIB" --transcriptome="$REF" --fastqs="$FQ" --sample="$SAMPLES" \
  --create-bam "$CREATE_BAM" --localcores "$NCORE" --localmem "$NMEM" --disable-ui
# only the outputs the pipeline reads are kept
KEEP="raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 molecule_info.h5 metrics_summary.csv web_summary.html"
[ "$CREATE_BAM" = true ] && KEEP="$KEEP possorted_genome_bam.bam possorted_genome_bam.bam.bai"
mkdir -p "$OUT"; for o in $KEEP; do ls_run cp -a "$W/$LIB/outs/$o" "$OUT/"; done
for o in $KEEP; do [ -s "$OUT/$o" ] || { echo "MISSING OUTPUT $o" >&2; exit 1; }; done
echo "outputs in $OUT; end $(date)"
