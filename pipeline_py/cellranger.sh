#!/bin/bash
# Cell Ranger 10.0.0 count for one 10x library (GSE174188 = Perez, GSE253568), the script that was run, with its paths
# taken from the environment variables of config.py: ATLAS_ROOT (default: the folder above pipeline_py/), REF_DIR
# (default <repo>/refdata/refdata-gex-GRCh38-2024-A) and PIPELINE_OUT (default <repo>/out), whose parent folder holds
# downloads/<DATASET>/ (sra/<SRR>/<SRR>.sra or fastq_ready/), cellranger/, logs/ and queue/. Options as run:
# GRCh38-2024-A, Cell Ranger defaults (no --include-introns, --expect-cells or --chemistry). Cell Ranger, the SRA
# Toolkit and pigz must be on PATH (or set CELLRANGER_MODULE and SRA_MODULE to environment modules to load).
# Warning: as run, the script deletes its inputs after success (FASTQ_DIR mode: once copied to the work folder).
# Usage: sbatch cellranger.sh <DATASET> <LIB> <RUN1> [<RUN2> ...]   (runs of ONE 10x library -> one Cell Ranger sample)
#        FASTQ_DIR=<dir> SAMPLES=<name>[,<name>] sbatch cellranger.sh <DATASET> <LIB>   (10x-named FASTQs, e.g. Perez/HCA)
# As run:
#   Perez, each h5ad library_uuid once its FASTQs were downloaded: FASTQ_DIR=<downloads>/Perez/<library>
#     SAMPLES=<its HCA FASTQ name(s), two for 272dbb99 and e16d2c39> CREATE_BAM=false sbatch --cpus-per-task=32 --mem=64G
#     --qos=<short or medium, from the run-time estimate> --time=<estimate> cellranger.sh Perez <library>
#   GSE253568, SRA runs prefetched into <downloads>/GSE253568/sra/: CREATE_BAM=true for the captures
#     00_genotype.py genotypes (false for capture1_well2) sbatch --cpus-per-task=32 --mem=64G cellranger.sh GSE253568 <capture> <SRR ...>
#     (capture1_well1: SRR27636185 SRR27636186 SRR27636167 SRR27636168; capture2-5 ran on 2026-10-25 with an earlier
#     revision of this script, the same cellranger count options, their BAMs kept for genotyping)
#   00_inputs.py reads the outputs from data/<disease>/<dataset>/cellranger/<library>/ (copy them there).
# FASTQ_DIR mode (2026-10-27): at job start the FASTQs are copied to the node-local disk, every file's byte size is
# checked against the source, and only then is the shared storage copy deleted (logged in logs/cleanup.log), so a running
# job no longer holds space on shared storage. A copy mismatch or too little local disk aborts the job with nothing
# deleted. If the job fails after the deletion, queue/requeue_library.sh puts the library's files back in the download
# queue and records the failure in queue/FAILED.txt (Cell Ranger is then resubmitted by hand). DRY_RUN=1 stops after
# the verified copy, before any deletion or Cell Ranger run (start-up test).
# CREATE_BAM=true to keep the BAM (only for genotyping, e.g. GSE266852); default false (README).
# Temporary data (FASTQ conversion, Cell Ranger work dir) live on the node-local disk; only outs/ go to shared storage.
#SBATCH --time=3-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --tmp=700G
#SBATCH --signal=B:USR1@300
set -euo pipefail
DS=$1; LIB=$2; shift 2; RUNS=("$@")
ROOT=${ATLAS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
B=$(dirname "${PIPELINE_OUT:-$ROOT/out}")   # downloads/, cellranger/, logs/, queue/
G=$B/downloads/$DS; OUT=$B/cellranger/$DS/$LIB
# own work folder (FASTQ conversion, Cell Ranger work dir), deleted on exit; ls_run runs a long command
JOBTMP=$(mktemp -d "${TMPDIR:-/tmp}/cellranger_XXXXXX"); LS_KILL=()
trap 'for p in "${LS_KILL[@]}"; do kill $p 2>/dev/null; done; rm -rf "$JOBTMP"' EXIT; trap 'exit 1' TERM INT
ls_run() { "$@"; }
LOC=$JOBTMP; FQ=$LOC/fastq; W=$LOC/work; mkdir -p $FQ $W
REF=${REF_DIR:-$ROOT/refdata/refdata-gex-GRCh38-2024-A}
MON=$B/logs/cellranger/$DS/peak_${LIB}_${SLURM_JOB_ID}.tsv
echo "node $(hostname) start $(date); $DS library $LIB runs ${RUNS[*]}; local $LOC ($(df -h $LOC | tail -1 | awk '{print $4}') free)"
( printf 'time\tlocal_GB\tlocal_files\n' > $MON; while true; do printf '%s\t%s\t%s\n' "$(date +%H:%M)" "$(du -sBG $LOC 2>/dev/null | cut -f1 | tr -d G)" "$(find $LOC 2>/dev/null | wc -l)" >> $MON; sleep 120; done ) & LS_KILL+=($!)
T0=$(date +%s)
CREATE_BAM=${CREATE_BAM:-false}; FASTQ_DIR=${FASTQ_DIR:-}; SAMPLES=${SAMPLES:-$LIB}
# Cell Ranger cores and memory follow the allocation (8 GB left for the system); defaults are the #SBATCH header
NCORE=${SLURM_CPUS_PER_TASK:-32}; NMEM=$(( ${SLURM_MEM_PER_NODE:-131072} / 1024 - 8 ))
echo "cellranger --localcores $NCORE --localmem $NMEM"
[ -z "${SRA_MODULE:-}" ] || module load "$SRA_MODULE"
# Cell Ranger must load and run on this node before any input is copied or deleted
{ [ -z "${CELLRANGER_MODULE:-}" ] || module load "$CELLRANGER_MODULE"; } && cellranger --version || { echo "ABORT: cellranger module not usable on $(hostname -s); nothing copied or deleted"; exit 1; }
i=0
if [ -n "$FASTQ_DIR" ]; then
  RUNS=(); SRC=$FASTQ_DIR
  need=$(du -sb $SRC | cut -f1); avail=$(df -B1 --output=avail $LOC | tail -1 | tr -d ' ')
  echo "FASTQ dir $SRC, samples $SAMPLES: $(ls $SRC | wc -l) files, $(numfmt --to=iec $need); local disk free $(numfmt --to=iec $avail)"
  # FASTQs + Cell Ranger work dir (peak 156-226 GB for 60-91 GB Perez libraries) must fit: abort before touching anything
  [ $avail -ge $(( need + 250 * 1024 * 1024 * 1024 )) ] || { echo "ABORT: local disk too small for the FASTQs + 250 GB work space; nothing deleted"; exit 1; }
  t=$(date +%s); ls_run /usr/bin/cp -p $SRC/* $FQ/
  bad=0; for f in $SRC/*; do b=$(basename $f); s1=$(stat -c %s $f); s2=$(stat -c %s $FQ/$b 2>/dev/null || echo missing)
    [ "$s1" = "$s2" ] || { echo "COPY MISMATCH $b: shared storage $s1 bytes, local $s2"; bad=1; }; done
  [ "$(ls $SRC | wc -l)" = "$(ls $FQ | wc -l)" ] || { echo "COPY MISMATCH: file counts differ"; bad=1; }
  [ $bad -eq 0 ] || { echo "ABORT: copy to local disk failed verification; shared storage FASTQs kept"; exit 1; }
  echo "FASTQs copied to $FQ in $(( ($(date +%s)-t)/60 )) min, sizes verified ($(numfmt --to=iec $need) on shared storage to be freed)"
  if [ "${DRY_RUN:-0}" = 1 ]; then echo "DRY RUN: copy verified; shared storage FASTQs kept; stopping before deletion and Cell Ranger"; echo DRY_RUN_DONE; exit 0; fi
  if [ -e $B/queue/keep/${DS}__${LIB} ]; then
    echo "inputs KEPT on shared storage: $B/queue/keep/${DS}__${LIB} exists (earlier failure, to be reviewed); Cell Ranger reads the local copy"
  else
  for f in $SRC/*; do echo "$(date '+%F %T') deleted $f ($(stat -c %s $f) bytes) | $DS $LIB FASTQ moved to node-local disk of job $SLURM_JOB_ID ($(hostname -s)) after size-verified copy | cr_library.sh" >> $B/logs/cleanup.log; done
  rm -rf "$SRC"; echo "FASTQ dir $SRC deleted from shared storage after verified copy"
  requeue_on_fail() { bash $B/queue/requeue_library.sh "$DS" "$LIB" "$SLURM_JOB_ID" "Cell Ranger job failed after its FASTQs were moved off shared storage ($1)" "$B/logs/cellranger/$DS/cr_${LIB}_${SLURM_JOB_ID}.log"; }
  LS_ON_FAIL=requeue_on_fail
  fi
fi
for r in "${RUNS[@]}"; do
  i=$((i+1)); L=$(printf 'L%03d' $i); t=$(date +%s)
  if [ -f $G/fastq_ready/${r}_1.fastq.gz ]; then
    ls_run cp $G/fastq_ready/${r}_1.fastq.gz $FQ/${LIB}_S1_${L}_R1_001.fastq.gz; ls_run cp $G/fastq_ready/${r}_2.fastq.gz $FQ/${LIB}_S1_${L}_R2_001.fastq.gz
  else
    # --include-technical: some SRA runs (e.g. GSE253568) store the barcode read as technical; reads are then
    # assigned by length: 26-28 bp = R1 (barcode + UMI), the longest = R2 (cDNA); index reads are dropped
    ls_run fasterq-dump --split-files --include-technical -e $NCORE -t $LOC/tmp_$r -O $FQ $G/sra/$r/$r.sra
    R1=""; R2=""; R2LEN=0
    for f in $FQ/${r}_*.fastq; do len=$(sed -n 2p $f | tr -d '\n' | wc -c); echo "  $(basename $f): read length $len"
      if [ $len -ge 26 ] && [ $len -le 28 ] && [ -z "$R1" ]; then R1=$f; elif [ $len -gt $R2LEN ]; then R2=$f; R2LEN=$len; fi; done
    [ -n "$R1" ] && [ -n "$R2" ] && [ "$R1" != "$R2" ] || { echo "cannot identify R1/R2 for $r"; exit 1; }
    for f in $FQ/${r}_*.fastq; do [ "$f" = "$R1" ] || [ "$f" = "$R2" ] || rm -f $f; done
    ls_run pigz -p $NCORE $R1 $R2
    mv $R1.gz $FQ/${LIB}_S1_${L}_R1_001.fastq.gz; mv $R2.gz $FQ/${LIB}_S1_${L}_R2_001.fastq.gz
    rm -rf $LOC/tmp_$r
  fi
  echo "$r -> $L done in $(( ($(date +%s)-t)/60 )) min; $(du -ch $FQ/${LIB}_S1_${L}_R*_001.fastq.gz | tail -1 | cut -f1)"
done
T1=$(date +%s); echo "FASTQ prep $(( (T1-T0)/60 )) min; library FASTQ total $(du -sh $FQ | cut -f1)"
cellranger --version   # loaded and checked at job start
cd $W
ls_run cellranger count --id=$LIB --transcriptome=$REF --fastqs=$FQ --sample=$SAMPLES \
  --create-bam $CREATE_BAM --localcores $NCORE --localmem $NMEM --disable-ui
T2=$(date +%s); echo "cellranger $(( (T2-T1)/60 )) min; work dir $(du -sh $W/$LIB | cut -f1), $(find $W/$LIB | wc -l) files"
# Only the outputs we keep go to shared storage (2026-10-26): both .h5 matrices, molecule_info.h5, the metrics and the web
# summary, plus the BAM when CREATE_BAM=true (genotyping). cloupe, analysis/ and the MEX folders (same data as the .h5)
# are not read by any pipeline and stay on the node-local disk, deleted with it.
KEEP="raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 molecule_info.h5 metrics_summary.csv web_summary.html"
[ "$CREATE_BAM" = true ] && KEEP="$KEEP possorted_genome_bam.bam possorted_genome_bam.bam.bai"
mkdir -p $OUT; for o in $KEEP; do ls_run cp -a $W/$LIB/outs/$o $OUT/; done
# Inputs are deleted only if every kept output is present and non-empty in $OUT
for o in $KEEP; do
  [ -s $OUT/$o ] || { echo "MISSING OUTPUT $o: inputs kept"; exit 1; }
done
echo "all required outputs present"
ls -la $OUT; du -sh $OUT; echo "kept files: $(find $OUT | wc -l)"
# a library marked by the failure watcher (queue/keep/<DATASET>__<LIB>) keeps its inputs until the mark is removed by hand
if [ -e $B/queue/keep/${DS}__${LIB} ]; then echo "inputs KEPT: $B/queue/keep/${DS}__${LIB} exists (earlier failure, to be reviewed)"; else
for r in "${RUNS[@]}"; do rm -rf $G/sra/$r $G/fastq_ready/${r}_*.fastq.gz; done; echo ".sra / ready FASTQ for these runs deleted from shared storage"
fi
echo "peak local use: $(sort -t$'\t' -k2,2n $MON | tail -1)"; echo "total $(( ($(date +%s)-T0)/60 )) min; end $(date)"
echo CELLRANGER_DONE
