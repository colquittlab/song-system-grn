#!/bin/bash
#SBATCH --job-name=gg_adult_export
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=900G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/export_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/export_%j.err
set -euo pipefail
echo "host=$(hostname) start=$(date)"
$HOME/miniforge3/envs/r44/bin/Rscript \
  /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/export_seurat.R
echo "end=$(date)"
