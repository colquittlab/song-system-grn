#!/bin/bash
#SBATCH --job-name=metaneighbor_pilot
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/logs/mn_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/logs/mn_%j.err
set -euo pipefail
mkdir -p /private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/logs
echo "host=$(hostname) start=$(date)"
$HOME/miniforge3/envs/metaneighbor/bin/Rscript /private/groups/colquittlab/finch-integration-toolkit/R/run_metaneighbor.R \
    /private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/data \
    /private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/metaneighbor/results
echo "end=$(date)"
