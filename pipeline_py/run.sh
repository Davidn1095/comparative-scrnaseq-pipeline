#!/bin/bash
# Run one step of the Python pipeline (pipeline_py/) in the atlas container (atlas-final-2026.10).
#   bash pipeline_py/run.sh <script> [args...]            run here (workstation, login node or inside a job)
#   bash pipeline_py/run.sh --sbatch "<sbatch options>" <script> [args...]
# --sbatch first copies pipeline_py/ to a read-only snapshot under $PIPELINE_OUT/snapshots/, so a queued or running job
# is never affected by later edits, then submits the step from it.
# Environment (defaults in brackets; <repo> is the folder above pipeline_py/):
#   ATLAS_ROOT          repository root [<repo>]
#   PIPELINE_OUT        output folder [<repo>/out]
#   REF_DIR             Cell Ranger GRCh38-2024-A reference [<repo>/refdata/refdata-gex-GRCh38-2024-A]
#   ATLAS_SIF           the container image [<repo>/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif]
#   PIPELINE_BIND       extra folders to bind into the container, comma-separated [none]; the repository, PIPELINE_OUT,
#                       REF_DIR and the job's temporary folder are always bound
#   APPTAINER_MODULE    environment module to load before running apptainer, if any [none]
#   PIPELINE_LOCAL_TMP  folder for per-job temporary files [$TMPDIR, else /tmp]
#   PIPELINE_CONSTRAINT, PIPELINE_GPU_CONSTRAINT  Slurm --constraint for CPU and GPU steps [none]
set -euo pipefail
# inside a job, Slurm runs a spool copy of this file: the snapshot is passed as PIPELINE_SNAPSHOT
HERE=${PIPELINE_SNAPSHOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
ROOT=${ATLAS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
OUT=${PIPELINE_OUT:-$ROOT/out}
REF=${REF_DIR:-$ROOT/refdata/refdata-gex-GRCh38-2024-A}
SIF=${ATLAS_SIF:-$ROOT/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif}

if [[ "${1:-}" == "--sbatch" ]]; then
  opts=$2; shift 2
  snap=$OUT/snapshots/$(date +%Y%m%d_%H%M%S)_$$
  mkdir -p "$snap" "$OUT/logs"
  (cd "$HERE" && tar --exclude=__pycache__ -cf - .) | (cd "$snap" && tar -xf -)
  chmod -R a-w "$snap"
  name=$(basename "$1" | sed 's/\.[^.]*$//')_${2:-all}
  # node type: results are identical on one node type but can differ between processor (and GPU) types. To pin one job,
  # put --constraint=<feature> in the options; to pin a whole run, export PIPELINE_CONSTRAINT (CPU steps) and
  # PIPELINE_GPU_CONSTRAINT (GPU steps). Options that already set a constraint are left as given.
  node=()
  if [[ "$opts" != *constraint* && "$opts" != *-C* ]]; then
    if [[ "$opts" == *gpu* ]]; then c=${PIPELINE_GPU_CONSTRAINT:-}; else c=${PIPELINE_CONSTRAINT:-}; fi
    [ -n "$c" ] && node=(--constraint="$c")
  fi
  # stdout carries only the job ID; anything else sbatch prints (some sites wrap it) goes to stderr
  # shellcheck disable=SC2086
  if ! sub=$(PIPELINE_SNAPSHOT="$snap" ATLAS_ROOT="$ROOT" PIPELINE_OUT="$OUT" REF_DIR="$REF" ATLAS_SIF="$SIF" \
      sbatch --parsable --signal=B:USR1@300 "${node[@]}" $opts \
      --job-name="$name" --output="$OUT/logs/${name}_%j.log" "$snap/run.sh" "$@"); then
    printf 'run.sh: sbatch failed: %s\n' "$(printf '%s' "$sub" | sed 's/\x1b\[[0-9;]*[A-Za-z]//g')" >&2; exit 1
  fi
  sub=$(printf '%s\n' "$sub" | sed 's/\x1b\[[0-9;]*[A-Za-z]//g')
  j=$(printf '%s\n' "$sub" | grep -E '^[0-9]+(;.*)?$' | tail -1 | cut -d';' -f1)   # --parsable: <id>[;<cluster>]
  printf '%s\n' "$sub" | grep -vE '^[0-9]+(;.*)?$|^$' >&2 || true
  [[ $j =~ ^[0-9]+$ ]] || { echo "run.sh: no job ID in the sbatch output" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\n' "$name" "$j" "$OUT/logs/${name}_$j.log" "$(date '+%F %T')" >> "$OUT/logs/submitted.tsv"
  echo "$j"
  exit
fi

script=$1; shift
if [ -n "${APPTAINER_MODULE:-}" ] && command -v module >/dev/null 2>&1; then module load "$APPTAINER_MODULE"; fi
command -v apptainer >/dev/null 2>&1 || { echo "run.sh: apptainer not found (set APPTAINER_MODULE or add it to PATH)" >&2; exit 1; }
[ -f "$SIF" ] || { echo "run.sh: container image not found: $SIF (set ATLAS_SIF)" >&2; exit 1; }
gpu=""; [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]] && gpu="--nv"
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
# home is mounted in the container: keep ~/.local (Python user site) and ~/.Renviron / ~/.Rprofile out of every step
export APPTAINERENV_PYTHONNOUSERSITE=1 APPTAINERENV_R_ENVIRON_USER=/dev/null APPTAINERENV_R_PROFILE_USER=/dev/null
# math-library threads capped at the job's CPUs: uncapped, R and NumPy threads oversubscribe the node
n=${SLURM_CPUS_PER_TASK:-2}
export APPTAINERENV_OMP_NUM_THREADS=$n APPTAINERENV_OPENBLAS_NUM_THREADS=$n APPTAINERENV_MKL_NUM_THREADS=$n APPTAINERENV_RCPP_PARALLEL_NUM_THREADS=$n
# own temporary folder, deleted on exit or on USR1 (--signal=B:USR1@300), TERM or INT; ls_run lets the signal act
source "$HERE/job_localscratch.sh"   # in pipeline_py/, so frozen with the snapshot
# numba (imported by decoupler in 04-06) needs a writable cache folder; the container's site-packages is read-only
mkdir -p "$TMPDIR/numba"; export APPTAINERENV_NUMBA_CACHE_DIR=$TMPDIR/numba
mkdir -p "$OUT"
bind=$ROOT,$OUT,$JOBTMP; [ -d "$REF" ] && bind=$bind,$REF; [ "$HERE" != "$ROOT/pipeline_py" ] && bind=$bind,$HERE
[ -n "${PIPELINE_BIND:-}" ] && bind=$bind,$PIPELINE_BIND
export APPTAINER_BIND=$bind
if [[ "$script" == *.sh ]]; then ATLAS_ROOT=$ROOT PIPELINE_OUT=$OUT REF_DIR=$REF ATLAS_SIF=$SIF ls_run bash "$HERE/$script" "$@"; exit; fi
case "$script" in
  *.R)  cmd=(Rscript --vanilla "$HERE/$script") ;;
  *.py) cmd=(python -u "$HERE/$script") ;;
esac
echo "host $(hostname) | snapshot $HERE | $script $* | $(date)"
# node type of this step: processor, and GPU when one is used
echo "cpu $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)${gpu:+ | gpu $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)}"
ls_run apptainer exec $gpu \
  --env PIPELINE_DIR="$HERE",PIPELINE_OUT="$OUT",ATLAS_ROOT="$ROOT",REF_DIR="$REF" \
  "$SIF" "${cmd[@]}" "$@"
