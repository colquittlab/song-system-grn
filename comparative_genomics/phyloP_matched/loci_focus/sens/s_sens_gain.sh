#!/bin/bash
#SBATCH -p medium
#SBATCH -t 8:00:00
#SBATCH --mem 64G
#SBATCH -c 2
#SBATCH -J sensgain
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/sens/gain-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
cd /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus
python sens_gain.py
echo SENSGAINDONE
