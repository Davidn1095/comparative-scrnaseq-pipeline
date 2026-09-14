#!/bin/bash
# Environment template for running the pipeline. Edit the values for your
# system, then source this file from the repository root before running the
# scripts:
#
#   source scripts/env_setup.sh
#
# ATLAS_ROOT         Repository root. data/ and results/ are read and written
#                    beneath it. Defaults to the current directory.
# SIF                Singularity/Apptainer image providing R 4.3.3,
#                    Bioconductor 3.18 and Python 3.10. Image definition:
#                    https://github.com/Davidn1095/research-env
# ATLAS_EXTRA_R_LIB  Optional extra R library searched first, for packages
#                    installed outside the image. Leave empty if unused.
# MSIGDB_RDS         Optional path of the MSigDB 2025.1 Hs cache in msigdbr's
#                    format. Defaults to msigdbr's user data directory.

export ATLAS_ROOT="${ATLAS_ROOT:-$(pwd)}"
export SIF="${SIF:-}"
export ATLAS_EXTRA_R_LIB="${ATLAS_EXTRA_R_LIB:-}"
export MSIGDB_RDS="${MSIGDB_RDS:-}"
[ -z "$MSIGDB_RDS" ] && unset MSIGDB_RDS

echo "Environment configured:"
echo "  ATLAS_ROOT:        $ATLAS_ROOT"
echo "  SIF:               ${SIF:-<not set>}"
echo "  ATLAS_EXTRA_R_LIB: ${ATLAS_EXTRA_R_LIB:-<not set>}"
echo "  MSIGDB_RDS:        ${MSIGDB_RDS:-<msigdbr default>}"
