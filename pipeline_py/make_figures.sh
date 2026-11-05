#!/bin/bash
# Manuscript figures and Table 1 of one atlas run (style in figure_style.R). Not a pipeline step: it reads the outputs of steps 02-06. Submit through run.sh, which freezes
# pipeline_py/ (these scripts included) in a read-only snapshot:
#   ATLAS_RUN=<run> bash pipeline_py/run.sh --sbatch "-c 2 --mem=32G -t 00:45:00" make_figures.sh
#   ATLAS_RUN = folder of 02_integration.py (default 02_integration); FIG_DRAFT=1 adds a red DRAFT mark to Figs 2-4 and
#   Table 1 (Fig 1 has none)
# Steps:
#   07_figure_data.py  figure inputs -> $PIPELINE_OUT/py/atlas/<run>/07_figures/data/   (atlas container)
#   08_figures.R       fig2_atlas, fig3_de, fig4_pathways (PDF) and table1.xlsx          (atlas container)
#   fig1_pipeline.tex  Fig 1 (PDF), pdflatex from TeX Live 2024 on the host (not in the container; TEX_MODULE, if set, is
#                      loaded as an environment module first)
#   PNG previews (150 dpi) of every figure, host Ghostscript if available (not in the container)
# Outputs, overwritten on each run: figure PDFs in manuscript/figures/; Table 1 in manuscript/tables/table1.xlsx (the
# journal takes each table as a separate Excel workbook); PNG previews in $PIPELINE_OUT/py/atlas/<run>/07_figures/.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)   # the snapshot, when run by run.sh --sbatch
ROOT=${ATLAS_ROOT:-$(dirname "$here")}
OUT=${PIPELINE_OUT:-$ROOT/out}
SIF=${ATLAS_SIF:-$ROOT/comparative-scrnaseq-pipeline_atlas-final-2026.10.sif}
export ATLAS_RUN=${ATLAS_RUN:-02_integration} PIPELINE_DIR=$here PIPELINE_OUT=$OUT ATLAS_ROOT=$ROOT
export FIG_OUT=$ROOT/manuscript/figures TABLE_OUT=$ROOT/manuscript/tables
PREVIEW=$OUT/py/atlas/$ATLAS_RUN/07_figures
mkdir -p "$FIG_OUT" "$TABLE_OUT" "$PREVIEW"
echo "figures: run $ATLAS_RUN | scripts $here | figures $FIG_OUT | table $TABLE_OUT | previews $PREVIEW"

# run.sh has set TMPDIR, APPTAINER_BIND and the APPTAINERENV_ variables; the container inherits this environment
for s in 07_figure_data.py 08_figures.R; do
  case $s in *.py) cmd=(python -u) ;; *.R) cmd=(Rscript --vanilla) ;; esac
  apptainer exec "$SIF" "${cmd[@]}" "$here/$s"
done

if [ -n "${TEX_MODULE:-}" ]; then type module >/dev/null 2>&1 || source /etc/profile >/dev/null 2>&1; module load "$TEX_MODULE"; fi
tex=${TMPDIR:-/tmp}/fig1_$$; mkdir -p "$tex"; cp "$here/fig1_pipeline.tex" "$tex/"
(cd "$tex" && pdflatex -interaction=nonstopmode -halt-on-error fig1_pipeline.tex >/dev/null) \
  || { tail -30 "$tex/fig1_pipeline.log"; exit 1; }
cp "$tex/fig1_pipeline.pdf" "$FIG_OUT/"

command -v gs >/dev/null 2>&1 || { echo "Ghostscript not found: no PNG previews"; echo "done $(date)"; exit 0; }
for f in fig1_pipeline fig2_atlas fig3_de fig4_pathways; do
  gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -sOutputFile="$PREVIEW/$f.png" "$FIG_OUT/$f.pdf"
  echo "saved $FIG_OUT/$f.pdf, preview $PREVIEW/$f.png"
done
echo "done $(date)"
