#!/bin/bash
#SBATCH --job-name=pc_samap_adult_gg
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/logs/pcsamap_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/logs/pcsamap_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
cd /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult
$HOME/miniforge3/envs/SAMap/bin/python per_cell_samap.py
