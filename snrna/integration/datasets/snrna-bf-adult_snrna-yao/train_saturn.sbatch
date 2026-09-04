#!/bin/bash
#SBATCH --job-name=saturn_bf_yao_adult
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-yao/logs/saturn_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-yao/logs/saturn_%j.err

# SATURN: adult finch song system (hybrid labels, lonStrDom2) x Yao et al. 2023 mouse
# telencephalon (ABC Atlas). No --score_adatas / --ct_map_path -- label transfer done
# post hoc from the joint embedding.

set -euo pipefail

BASE=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-yao
PY=$HOME/miniforge3/envs/saturn/bin/python
SATURN_REPO=$HOME/group/repos/SATURN
SATURN=$SATURN_REPO/train-saturn.py

N_MACRO=${N_MACRO:-2000}
N_HV=${N_HV:-8000}
SEED=${SEED:-0}
EPOCHS=${EPOCHS:-50}
POLL=${POLL:-25}
IN_LABEL=${IN_LABEL:-cluster}
REF_LABEL=${REF_LABEL:-cluster}
RUN=${RUN:-macro${N_MACRO}_hv${N_HV}_seed${SEED}_ep${EPOCHS}}
WORK=$BASE/train/saturn_results/$RUN
mkdir -p "$WORK"

echo "host=$(hostname) run=$RUN start=$(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv

DATA_CSV=${DATA_CSV:-$BASE/data/bf_yao_adult_hybrid_run.csv}

$PY "$SATURN" \
    --in_data="$DATA_CSV" \
    --in_label_col="$IN_LABEL" \
    --ref_label_col="$REF_LABEL" \
    --num_macrogenes="$N_MACRO" \
    --hv_genes="$N_HV" \
    --centroids_init_path="$WORK/bf_yao_adult_hybrid_centroids.pkl" \
    --work_dir="$WORK" \
    --seed="$SEED" \
    --epochs="$EPOCHS" \
    --polling_freq="$POLL" \
    --equalize_triplets_species \
    --balance_species_cells

echo "end=$(date)"
ls -la "$WORK/saturn_results/" || true
