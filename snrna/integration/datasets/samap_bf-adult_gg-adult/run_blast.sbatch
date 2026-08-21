#!/bin/bash
#SBATCH --job-name=blast_bfgg_adult
#SBATCH --partition=medium
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/logs/blast_%j.out
#SBATCH --error=/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult/logs/blast_%j.err
set -euo pipefail
cd /private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult
BB=$HOME/miniforge3/envs/SAMap/bin   # blast is not on compute nodes; use the env copy
$BB/makeblastdb -in data/bf_adult.faa -dbtype prot -out data/bfdb >/dev/null
$BB/makeblastdb -in data/gg_grcg7b.faa -dbtype prot -out data/ggdb >/dev/null
echo "bf->gg $(date)"
$BB/blastp -query data/bf_adult.faa -db data/ggdb -outfmt 6 -max_hsps 1 -evalue 1e-6 \
       -num_threads 48 -out maps/bfgg/bf_to_gg.txt
echo "gg->bf $(date)"
$BB/blastp -query data/gg_grcg7b.faa -db data/bfdb -outfmt 6 -max_hsps 1 -evalue 1e-6 \
       -num_threads 48 -out maps/bfgg/gg_to_bf.txt
echo "done $(date)"; wc -l maps/bfgg/*.txt
