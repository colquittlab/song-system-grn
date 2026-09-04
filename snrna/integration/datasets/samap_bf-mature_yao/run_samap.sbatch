#!/bin/bash
#SBATCH --job-name=samap_yao
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-mature_yao/logs/samap_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-mature_yao/logs/samap_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
cd /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-mature_yao
echo "host=$(hostname) start=$(date)"
$HOME/miniforge3/envs/SAMap/bin/python run_samap.py --ncpus 32
echo "end=$(date)"
