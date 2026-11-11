#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 32G
#SBATCH -J paired
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/paired-%j.out
set -eu
BW=/private/groups/colquittlab/opt/ucsc_utilities/bigWigAverageOverBed
H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; cd /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus
for t in wholetree subtree222; do
  $BW $H/${t}_perbase_lsd2.bw gwrest.bed ${t}_gwrest.tab
  echo "  $t gwrest done"
done
echo DONE
