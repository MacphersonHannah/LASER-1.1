#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Auto activates the LASER runner env, no manual `conda activate` needed.
if ! command -v conda &> /dev/null; then
  echo "ERROR: conda not found on PATH." >&2
  exit 1
fi
eval "$(conda shell.bash hook)"
conda activate "${SCRIPT_DIR}/envs_prebuilt/core_utils_env"

CACHE_BASE="$(pwd)/.sing"
mkdir -p "$CACHE_BASE/cache" "$CACHE_BASE/tmp"


export SINGULARITY_CACHEDIR="$CACHE_BASE/cache"
export SINGULARITY_TMPDIR="$CACHE_BASE/tmp"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
export TMPDIR="$SINGULARITY_TMPDIR"

CONDA_CACHE_BASE="${CACHE_BASE}/conda_cache"
export CONDA_PKGS_DIRS="${CONDA_CACHE_BASE}/conda-pkgs"
export NXF_CONDA_CACHEDIR="${CONDA_CACHE_BASE}/nxf-conda"
mkdir -p "$CONDA_PKGS_DIRS" "$NXF_CONDA_CACHEDIR"

RUN_NAME="${1:-default}"
CONFIG_FILE="${2:-}"
if [[ -n "$CONFIG_FILE" ]]; then
  CONFIG_FILE="$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${SCRIPT_DIR}/run_${RUN_NAME}"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

export NXF_HOME="$(pwd)/.nextflow_${RUN_NAME}"

WORKDIR="work_${RUN_NAME}"
OUTDIR="${RUN_NAME}"
TS="$(date +%Y%m%d_%H%M%S)"

echo ">>> Running Nextflow with run name: $RUN_NAME"
echo ">>> Workdir: $WORKDIR"
echo ">>> Outdir:  $OUTDIR"
echo ">>> Stamp:   $TS"

# Copy config into output directory for reproducibility
if [[ -n "$CONFIG_FILE" ]]; then
  mkdir -p "$OUTDIR"
  cp "$CONFIG_FILE" "${OUTDIR}/config_used_${RUN_NAME}_${TS}.config"
fi

# --outdir below is always derived from RUN_NAME (this script's first argument)
NF_CMD=(
  nextflow run "${SCRIPT_DIR}/main.nf"
  -profile conda
  -work-dir "$WORKDIR"
  --outdir "$OUTDIR"
  -with-report   "${OUTDIR}/report_${RUN_NAME}_${TS}.html"
  -with-timeline "${OUTDIR}/timeline_${RUN_NAME}_${TS}.html"
  -with-trace    "${OUTDIR}/trace_${RUN_NAME}_${TS}.txt"
  -resume
)

if [[ -n "$CONFIG_FILE" ]]; then
  NF_CMD+=(-c "$CONFIG_FILE")
fi

"${NF_CMD[@]}"

### ./run_LASER.sh run1 run1.config   # goes into work_run1 and run1 folders, uses run1.config
