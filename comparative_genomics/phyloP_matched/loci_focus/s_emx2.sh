#!/bin/bash
#SBATCH -p medium
#SBATCH -t 8:00:00
#SBATCH --mem 32G
#SBATCH -c 8
#SBATCH -J emx2tf
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/emx2-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; FA=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.fa
cd $F
cat e_ra.bed e_ast.bed | sort -k1,1 -k2,2n > e_all.bed
bedtools shuffle -i e_all.bed -g ar_gen.sizes -excl e_all.bed -chrom -seed 11 \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"BG."NR}' | sort -k1,1 -k2,2n > e_bg.bed
for s in ra ast bg; do
  bedtools getfasta -fi $FA -bed e_$s.bed -nameOnly -fo e_$s.fa
  echo "  $s: $(grep -c '^>' e_$s.fa)"
done
for tf in Fosl2 Mafb Nfia Nr2e1 Rfx4; do
  for s in ra ast bg; do
    fimo --verbosity 1 --thresh 1e-4 --text tf_$tf.meme e_$s.fa > efimo_${tf}_${s}.txt 2>/dev/null
  done
  echo "  $tf done"
done
