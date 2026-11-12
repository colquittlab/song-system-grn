#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 32G
#SBATCH -c 2
#SBATCH -J phylop_density
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/phylop_density-%j.out
set -euo pipefail
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; set +u; conda activate hic; set -u
cd /private/groups/colquittlab/song-system-grn-worktrees/loci-focus/comparative_genomics/phyloP_matched/loci_focus
python phylop_density.py ${TRACKS:+--tracks $TRACKS} ${CENTER:+--center $CENTER} ${SMOOTH:+--smooth $SMOOTH} ${SPANS:-domain 1p5mb gw}
