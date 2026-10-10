#!/bin/bash
#SBATCH --job-name=ct_hybrid_score
#SBATCH --partition=medium
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ct_score_%A_%a.out
#SBATCH --error=logs/ct_score_%A_%a.err
# One array task per part of the motif list; submit_cistarget_build.sh sets --array=1-NPARTS and exports CT_ROOT and NPARTS.
# Locally the whole scoring took ~19 h on 40 cores (~8.6 motifs/min), so 6 parts is ~3.2 h each on 40 cores; --time leaves margin.
# A part that already finished is skipped (see score_part.sh), so a requeued task is cheap.
set -euo pipefail
set +u   # conda's activate scripts read unset variables
source /private/home/${USER}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
bash "${CT_ROOT}/prism/score_part.sh" "${SLURM_ARRAY_TASK_ID}" "${NPARTS}" "${SLURM_CPUS_PER_TASK}"
