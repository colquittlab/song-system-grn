#!/bin/bash
#SBATCH --job-name=ct_hybrid_combine
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ct_combine_%j.out
#SBATCH --error=logs/ct_combine_%j.err
# Runs after every scoring task has succeeded (submit_cistarget_build.sh sets the dependency). Each full matrix is
# 4 bytes x 499,348 regions x 10,249 motifs = 20 GB; combining and ranking hold a few of them at once.
set -euo pipefail
set +u
source /private/home/${USER}/miniforge3/etc/profile.d/conda.sh
conda activate scenicplus
set -u
bash "${CT_ROOT}/prism/combine_and_rank.sh" "${NPARTS}"
