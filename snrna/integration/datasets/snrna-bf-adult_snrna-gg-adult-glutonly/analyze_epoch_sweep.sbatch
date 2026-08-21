#!/bin/bash
#SBATCH --job-name=an_ep
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/logs/analyze_ep%x_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/logs/analyze_ep%x_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16

ADATA=${ADATA:?set ADATA}
OUT=${OUT:?set OUT}
GGA=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly

echo "host=$(hostname) start=$(date) ADATA=$ADATA"
$HOME/miniforge3/envs/saturn/bin/python \
  /private/groups/colquittlab/saturn/snrna-bf-dev_snrna-lamanno-dev/analyze_integration.py \
  --adata "$ADATA" --out_dir "$OUT" --label_col labels --strip_species_prefix \
  --ref_species chicken --query_species finch \
  --annotation "$GGA/data/gg_glutonly_label_annotation.csv"
echo "end=$(date)"
