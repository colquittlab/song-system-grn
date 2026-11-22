#!/bin/bash
#SBATCH -p medium
#SBATCH -t 6:00:00
#SBATCH --mem 24G
#SBATCH -c 8
#SBATCH -J tfscan
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/tf-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
export PATH=/private/groups/colquittlab/opt/cactus-bin-v3.1.2/bin:$PATH
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; FA=/private/warm-archive/colquittlab/assembly/cactus/avian363/vgp577/fa/lonStrDom2.fa
cd $F
cat spec_GlutRA_gw.bed spec_Astro1_gw.bed | sort -k1,1 -k2,2n > spec_all.bed
bedtools shuffle -i spec_all.bed -g ar_gen.sizes -excl spec_all.bed -chrom -seed 7 \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"bg"NR}' | sort -k1,1 -k2,2n > spec_bg.bed
for s in GlutRA_gw Astro1_gw; do bedtools getfasta -fi $FA -bed spec_$s.bed -nameOnly -fo spec_$s.fa; done
bedtools getfasta -fi $FA -bed spec_bg.bed -nameOnly -fo spec_bg.fa
for tf in Fosl2 Mafb Nfia Nr2e1 Rfx4; do
  for s in GlutRA_gw Astro1_gw bg; do
    fimo --verbosity 1 --thresh 1e-4 --text tf_$tf.meme spec_$s.fa > fimo_${tf}_${s}.txt 2>/dev/null
    echo "  $tf $s: $(( $(wc -l < fimo_${tf}_${s}.txt) - 1 ))"
  done
done
