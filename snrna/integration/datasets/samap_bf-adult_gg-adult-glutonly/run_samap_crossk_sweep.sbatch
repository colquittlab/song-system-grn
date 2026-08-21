#!/bin/bash
#SBATCH --job-name=samap_ck
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly/logs/sweep_ck%x_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly/logs/sweep_ck%x_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
cd /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly
echo "host=$(hostname) start=$(date) CROSSK=${CROSSK}"
$HOME/miniforge3/envs/SAMap/bin/python run_samap.py --ncpus 32 --numiters 3 --crossk "${CROSSK}" \
  --out_dir "results_sweep/ck${CROSSK}"
echo "end=$(date)"
