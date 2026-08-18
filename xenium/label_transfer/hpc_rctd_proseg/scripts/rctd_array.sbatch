#!/bin/bash
#SBATCH --job-name=rctd-proseg
#SBATCH --array=1-79%40
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/rctd_%A_%a.out
#SBATCH --error=logs/rctd_%A_%a.err
# NOTE on the two lines above: SLURM opens these files before this script's own
# body runs, so `logs/` (relative to wherever `sbatch` is invoked from) must
# already exist at submit time. It ships pre-created (with a `.keep`
# placeholder, since rclone/rsync can drop empty directories) by
# 01_prepare_chunks.R. If submission fails with "couldn't open output file",
# run `mkdir -p logs` yourself first — most likely the sync dropped it.
# ---------------------------------------------------------------------------
# Confirmed on phoenix (prism/Slurm) via srun + sacct: this cluster needs no
# --account, and partition "medium" accepts and dispatches this job's
# resources (4 CPUs, 12G, 1h) — medium allows up to 12h. Adjust --partition
# (and re-add --account) if deploying to a different cluster.
#SBATCH --partition=medium
# ---------------------------------------------------------------------------
#
# RCTD over the PROSEG-resegmented Xenium data (all ten sections). Same
# sizing rationale as the delivered-segmentation run in
# xenium/label_transfer/hpc/ (measured locally, not guessed):
#   fitPixels on 20k cells = ~17 min compute + ~4 min gather = ~21 min wall
#   spacexr uses ~2.4 cores effectively regardless of how many it is given, so
#   4 cpus-per-task is the efficient point; extra cores are wasted.
#   Memory: a pre-sliced 20k-cell chunk peaks near 6 GB; 12G leaves headroom.
#   --array count depends on the exact chunk count 01_prepare_chunks.R reports
#   (~1.59M proseg cells / 20k ~= 80, close to but not guaranteed identical to
#   the delivered-segmentation run's count) — check its "DONE" message and
#   update --array=1-N above before submitting if it differs from 80.
#
#   %40 throttles to 40 concurrent tasks. Raise or remove for a faster run.
#
# Submit (self-contained directory, built by 01_prepare_chunks.R, moved with
# rclone — see README.md):
#   cd /scratch/$USER/hpc_rctd_proseg
#   mkdir -p logs                                   # defensive, see note above
#   sbatch scripts/rctd_array.sbatch .
#
# Resume after partial failure (only missing indices):
#   sbatch --array=$(Rscript scripts/04_missing.R .) scripts/rctd_array.sbatch .

set -euo pipefail

BASE="${1:?usage: sbatch rctd_array.sbatch <hpc_rctd_proseg_dir>}"
BASE="$(cd "$BASE" && pwd)"
mkdir -p "$BASE/logs" "$BASE/results"   # results/: nothing to lose by re-asserting; logs/: helps a resubmit
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
