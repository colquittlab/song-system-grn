#!/bin/bash
#SBATCH -p medium
#SBATCH -t 8:00:00
#SBATCH --mem 32G
#SBATCH -c 8
#SBATCH -J arscan
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/ar-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; FA=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.fa
cd $F
awk '$5=="RA_over_Arco"'      allelem_cls.bed | cut -f1-4 > ar_ra.bed
awk '$5=="GABA4_over_GABA1"'  allelem_cls.bed | cut -f1-4 > ar_g4.bed
cat ar_ra.bed ar_g4.bed | sort -k1,1 -k2,2n > ar_all.bed
# matched background: same sizes, same chromosomes, not overlapping any DAR
grep -E "^chr([0-9]+|[0-9]+A|Z)\b" /private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.sizes | sort -k1,1 > ar_gen.sizes
bedtools shuffle -i ar_all.bed -g ar_gen.sizes -excl ar_all.bed -chrom -seed 1 \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"bg"NR}' | sort -k1,1 -k2,2n > ar_bg.bed
for s in ra g4 bg; do
  bedtools getfasta -fi $FA -bed ar_$s.bed -nameOnly -fo ar_$s.fa
  echo "  $s: $(grep -c '^>' ar_$s.fa) sequences"
done
for m in jaspar hoco; do
  for s in ra g4 bg; do
    fimo --verbosity 1 --thresh 1e-4 --text ar_$m.meme ar_$s.fa > fimo_${m}_${s}.txt 2>/dev/null
    echo "  $m $s: $(( $(wc -l < fimo_${m}_${s}.txt) - 1 )) hits"
  done
done
