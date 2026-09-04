#!/bin/bash
#SBATCH --job-name=saturn_bf_gg_adult
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/saturn_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/logs/saturn_%j.err

# SATURN: adult finch song system (Bengalese finch, lonStrDom2) x adult chicken
# (Zaremba et al. 2025). No --score_adatas / --ct_map_path -- the cross-species
# correspondence is the unknown; label transfer is done post hoc from the joint embedding.

set -euo pipefail

BASE=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult
PY=$HOME/miniforge3/envs/saturn/bin/python
SATURN_REPO=$HOME/group/repos/SATURN
if [ -n "${BP:-}" ]; then
    SATURN=$BASE/train-saturn-bpfix.py
    export PYTHONPATH="$SATURN_REPO${PYTHONPATH:+:$PYTHONPATH}"
else
    SATURN=$SATURN_REPO/train-saturn.py
fi

N_MACRO=${N_MACRO:-2000}
N_HV=${N_HV:-8000}
SEED=${SEED:-0}
EPOCHS=${EPOCHS:-50}
POLL=${POLL:-25}
IN_LABEL=${IN_LABEL:-cluster}
REF_LABEL=${REF_LABEL:-cluster}
BP=${BP:-}
RUN=${RUN:-macro${N_MACRO}_hv${N_HV}_seed${SEED}_ep${EPOCHS}${BP:+_bp}$([ "$IN_LABEL" != cluster ] && echo _$IN_LABEL)}
WORK=$BASE/train/saturn_results/$RUN
mkdir -p "$WORK"

echo "host=$(hostname) run=$RUN start=$(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv

DATA_CSV=${DATA_CSV:-$BASE/data/bf_gg_adult_run.csv}

$PY "$SATURN" \
    --in_data="$DATA_CSV" \
    --in_label_col="$IN_LABEL" \
    --ref_label_col="$REF_LABEL" \
    --num_macrogenes="$N_MACRO" \
    --hv_genes="$N_HV" \
    --centroids_init_path="$WORK/bf_gg_adult_centroids.pkl" \
    --work_dir="$WORK" \
    --seed="$SEED" \
    --epochs="$EPOCHS" \
    --polling_freq="$POLL" \
    --equalize_triplets_species \
    --balance_species_cells \
    ${BP:+--balance_pretrain}

echo "end=$(date)"
ls -la "$WORK/saturn_results/" || true
