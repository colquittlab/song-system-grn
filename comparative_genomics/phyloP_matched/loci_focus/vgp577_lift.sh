#!/bin/bash
# lonStrDom2 elements -> bLonStr1 (liftOver chain) -> zebra finch GCA_048771995.1 (halLiftover, VGP 577-way)
#SBATCH -p medium
#SBATCH -t 12:00:00
#SBATCH --mem 32G
#SBATCH -c 1
#SBATCH -J vgplift
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/vgp577/logs/lift-%j.out
set -eu
W=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/vgp577
CH=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/chain/lonStrDom2ToBLonStr1.over.chain
HAL=/private/warm-archive/colquittlab/assembly/cactus/vgp577/vgp-577way-v1.hal
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:/private/groups/colquittlab/opt/ucsc_utilities:$PATH
liftOver -minMatch=0.5 $W/elements_lsd2.bed $CH $W/elements_bLonStr1.bed $W/elements_bLonStr1.unmapped
echo "bLonStr1: $(wc -l < $W/elements_bLonStr1.bed) lifted, $(grep -vc '^#' $W/elements_bLonStr1.unmapped) unmapped"
halLiftover $HAL GCA_046129695.1 $W/elements_bLonStr1.bed GCA_048771995.1 $W/elements_zf_frag.bed
echo "zebra finch fragments: $(wc -l < $W/elements_zf_frag.bed)"
