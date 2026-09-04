#!/bin/bash
#SBATCH --job-name=va_gg_adult
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/validate_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/validate_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
$HOME/miniforge3/envs/saturn/bin/python /private/groups/colquittlab/saturn/snrna-bf-dev_snrna-lamanno-dev/validate_transfer.py \
  --saturn_h5ad "$ADATA" --out_dir "$OUT" --ref_species chicken --query_species finch \
  --label_key "${LABELKEY:-ref_labels}" \
  --query_h5ad "$QH5" --ref_h5ad "$REFH5"
