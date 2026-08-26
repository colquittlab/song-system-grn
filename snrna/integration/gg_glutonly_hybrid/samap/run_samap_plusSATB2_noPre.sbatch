#!/bin/bash
#SBATCH --job-name=samap_plusSATB2_noPre
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/logs/samap_plusSATB2_noPre_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/logs/samap_plusSATB2_noPre_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
cd /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly
echo "host=$(hostname) start=$(date)"
$HOME/miniforge3/envs/SAMap/bin/python run_samap.py --ncpus 32 \
    --finch_h5ad /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data/bf_adult_glut_hybrid_plusSATB2.h5ad \
    --chicken_h5ad /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult-glutonly/data/gg_adult_ex_noPre.h5ad \
    --out_dir /private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/samap/results/plusSATB2_noPre
echo "end=$(date)"
