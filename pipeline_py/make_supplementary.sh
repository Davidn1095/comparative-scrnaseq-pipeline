#!/bin/bash
# Supplementary tables and figures of one atlas run. Not a pipeline step: it reads the outputs of steps 03-06 and 04b. Submit through run.sh, which freezes pipeline_py/
# in a read-only snapshot:
#   ATLAS_RUN=<run> bash pipeline_py/run.sh --sbatch "-c 4 --mem=64G -t 01:00:00" make_supplementary.sh
#   ATLAS_RUN = folder of 02_integration.py (default 02_integration); FIG_DRAFT=1 adds a red DRAFT mark
# Steps:
#   09_supplementary_data.py  per-donor cell-type counts -> $PIPELINE_OUT/py/atlas/<run>/09_supplementary/data/  (atlas container)
#   10_supplementary.R        tables (xlsx) and figures (PDF) -> manuscript/supplementary/               (atlas container)
#   PNG previews (150 dpi) of the figures, host Ghostscript if available -> $PIPELINE_OUT/py/atlas/<run>/09_supplementary/
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)   # the snapshot, when run by run.sh --sbatch
ROOT=${ATLAS_ROOT:-$(dirname "$here")}
OUT=${PIPELINE_OUT:-$ROOT/out}
SIF=${ATLAS_SIF:-$ROOT/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif}
export ATLAS_RUN=${ATLAS_RUN:-02_integration} PIPELINE_DIR=$here PIPELINE_OUT=$OUT ATLAS_ROOT=$ROOT
export SUPP_OUT=$ROOT/manuscript/supplementary
PREVIEW=$OUT/py/atlas/$ATLAS_RUN/09_supplementary
mkdir -p "$SUPP_OUT" "$PREVIEW"
echo "supplementary: run $ATLAS_RUN | scripts $here | output $SUPP_OUT | previews $PREVIEW"

# run.sh has set TMPDIR, APPTAINER_BIND and the APPTAINERENV_ variables; the container inherits this environment
for s in 09_supplementary_data.py 10_supplementary.R; do
  case $s in *.py) cmd=(python -u) ;; *.R) cmd=(Rscript --vanilla) ;; esac
  apptainer exec "$SIF" "${cmd[@]}" "$here/$s"
done

command -v gs >/dev/null 2>&1 || { echo "Ghostscript not found: no PNG previews"; echo "done $(date)"; exit 0; }
for f in supp_fig_composition_cohorts supp_fig_de_cohorts supp_fig_atlas; do
  gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -sOutputFile="$PREVIEW/$f.png" "$SUPP_OUT/$f.pdf"
  echo "saved $SUPP_OUT/$f.pdf, preview $PREVIEW/$f.png"
done
echo "done $(date)"
