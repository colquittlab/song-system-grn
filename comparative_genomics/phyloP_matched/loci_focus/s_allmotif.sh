#!/bin/bash
#SBATCH -p medium
#SBATCH -t 12:00:00
#SBATCH --mem 96G
#SBATCH -c 16
#SBATCH -J allmotif
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/allmotif-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:/private/groups/colquittlab/opt/ucsc_utilities:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; FA=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.fa
cd $F
# all RA and GABA-4 DARs
cut -f1-4 allelem_cls.bed | sort -k1,1 -k2,2n > am_dars.bed
echo "  DARs: $(wc -l < am_dars.bed)"
[ -s am_dars.fa ] || bedtools getfasta -fi $FA -bed am_dars.bed -nameOnly -fo am_dars.fa
if [ ! -s am_hits.txt ]; then
  fimo --verbosity 1 --thresh 1e-4 --text mask_meta192.meme am_dars.fa > am_hits.txt
fi
echo "  hits: $(( $(wc -l < am_hits.txt) - 1 ))"
# score every DAR once
for t in subtree222 wholetree; do
  [ -s ${t}_am_dars.tab ] || bigWigAverageOverBed $H/${t}_perbase_lsd2.bw am_dars.bed ${t}_am_dars.tab
done
echo "  DAR scoring done"
python3 build_hits.py
for t in subtree222 wholetree; do
  [ -s ${t}_am_hits.tab ] || bigWigAverageOverBed $H/${t}_perbase_lsd2.bw am_merged.bed ${t}_am_hits.tab
done
echo ALLDONE
