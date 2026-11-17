#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 16G
#SBATCH -c 1
#SBATCH -J domdepth
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/domain_depth-%j.out
set -euo pipefail
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; set +u; conda activate hic; set -u
# DHX35-TOP1 domain +/-10 kb
python /private/groups/colquittlab/song-system-grn-worktrees/loci-focus/comparative_genomics/phyloP_matched/loci_focus/domain_depth.py chr20 12541036 13557533 /private/warm-archive/colquittlab/assembly/cactus/avian363/maf/maf_chunked/chr20/chunk_chr20_chr20.13.maf /private/warm-archive/colquittlab/assembly/cactus/avian363/maf/maf_chunked/chr20/chunk_chr20_chr20.14.maf
