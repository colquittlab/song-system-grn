#!/bin/bash
#SBATCH --job-name=rctd
#SBATCH --array=1-80%40
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/rctd_%A_%a.out
#SBATCH --error=logs/rctd_%A_%a.err
# ---------------------------------------------------------------------------
# Confirmed on phoenix (prism/Slurm) via srun + sacct: this cluster needs no
# --account, and partition "medium" accepts and dispatches this job's
# resources (4 CPUs, 12G, 1h) — medium allows up to 12h. Adjust --partition
# (and re-add --account) if deploying to a different cluster.
#SBATCH --partition=medium
# ---------------------------------------------------------------------------
#
# Sizing rationale (measured locally, not guessed):
#   fitPixels on 20k cells = ~17 min compute + ~4 min gather = ~21 min wall
#   spacexr uses ~2.4 cores effectively regardless of how many it is given, so
#   4 cpus-per-task is the efficient point; extra cores are wasted.
#   Memory: a pre-sliced 20k-cell chunk peaks near 6 GB; 12G leaves headroom.
#   --time=01:00:00 is ~3x the measured 21 min.
#
#   %40 throttles to 40 concurrent tasks. Raise or remove for a faster run;
#   80 concurrent finishes the whole dataset in roughly one task's wall time.
#
# Submit:
#   mkdir -p logs
#   sbatch rctd_array.sbatch /path/to/hpc_rctd
#
# Resume after partial failure (only missing indices):
#   sbatch --array=$(Rscript 04_missing.R /path/to/hpc_rctd) rctd_array.sbatch /path/to/hpc_rctd

set -euo pipefail

BASE="${1:?usage: sbatch rctd_array.sbatch <hpc_rctd_dir>}"
BASE="$(cd "$BASE" && pwd)"
CHUNK_DIR="$BASE/chunks"
OUT_DIR="$BASE/results"

# --- environment -----------------------------------------------------------
# Adjust to however conda/micromamba is exposed on your cluster.
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate xenium-rctd
else
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate xenium-rctd
fi

# spacexr spawns PSOCK workers; keep BLAS single-threaded so the two do not
# oversubscribe the allocated cpus.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# NOTE: do NOT derive this from ${BASH_SOURCE[0]}. Slurm copies the batch
# script to a node-local spool dir (/var/lib/slurm/slurmd/job<id>/), so
# BASH_SOURCE resolves there and the sibling R scripts are not present.
# Caused 160/160 task failures across two runs on phoenix (all exited in 5-7s
# with "cannot open file '.../job<id>/02_fit_chunk.R'") before this fix.
SCRIPT_DIR="$BASE/scripts"
if [ ! -f "$SCRIPT_DIR/02_fit_chunk.R" ]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$PWD}/scripts"
fi
[ -f "$SCRIPT_DIR/02_fit_chunk.R" ] || { echo "cannot locate 02_fit_chunk.R (tried $BASE/scripts and ${SLURM_SUBMIT_DIR:-$PWD}/scripts)" >&2; exit 1; }

Rscript "$SCRIPT_DIR/02_fit_chunk.R" \
    "$CHUNK_DIR" "$OUT_DIR" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
