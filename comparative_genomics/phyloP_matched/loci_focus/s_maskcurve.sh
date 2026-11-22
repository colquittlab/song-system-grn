#!/bin/bash
#SBATCH -p medium
#SBATCH -t 10:00:00
#SBATCH --mem 48G
#SBATCH -c 8
#SBATCH -J maskcurve
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/maskcurve-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; FA=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.fa
cd $F
# sequences for the 416 locus DARs, named by element so hits map back
cut -f1-4 loci_dars.bed > mc_dars.bed
bedtools getfasta -fi $FA -bed mc_dars.bed -nameOnly -fo mc_dars.fa
echo "  sequences: $(grep -c '^>' mc_dars.fa)"
for t in ar1 sw37 meta192 all7446; do
  if [ ! -s mc_$t.txt ]; then
    fimo --verbosity 1 --thresh 1e-4 --text mask_$t.meme mc_dars.fa > mc_$t.txt
  fi
  echo "  $t hits: $(( $(wc -l < mc_$t.txt) - 1 ))"
done
echo SCANDONE
