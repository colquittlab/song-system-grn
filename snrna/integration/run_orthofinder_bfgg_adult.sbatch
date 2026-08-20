#!/bin/bash
#SBATCH --job-name=of_bfgg_adult
#SBATCH --partition=long
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/colquittlab/saturn/zaremba_composite/logs/ofggadult_%j.out
#SBATCH --error=/private/groups/colquittlab/saturn/zaremba_composite/logs/ofggadult_%j.err
set -euo pipefail
cd /private/groups/colquittlab/saturn/zaremba_composite
# of254b: python 3.10 + numpy 1.24 + scipy 1.10 -- the bioconda default env ships python
# 3.14, whose numpy 2.x breaks OrthoFinder 2.5.4 at the gene-distance step.
export PATH="$HOME/miniforge3/envs/of254b/bin:$PATH"
rm -rf orthofinder_bfgg_adult/proteomes/OrthoFinder
$HOME/miniforge3/envs/of254b/bin/orthofinder -f orthofinder_bfgg_adult/proteomes -t 32 -a 1
echo "done $(date)"
find orthofinder_bfgg_adult/proteomes/OrthoFinder -name "*__v__*.tsv"
