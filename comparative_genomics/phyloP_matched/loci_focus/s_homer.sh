#!/bin/bash
#SBATCH -p medium
#SBATCH -t 8:00:00
#SBATCH --mem 32G
#SBATCH -c 8
#SBATCH -J homerscan
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/homer-%j.out
set -eu
source /private/home/bcolquit/miniforge3/etc/profile.d/conda.sh; conda activate crested; set +u
F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus; D=/private/home/bcolquit/group/motifs/homer/ra_pos_vs_neg/homerResults_manual
cd $F
for m in motif2 motif4 motif8; do
  for s in ra g4 bg; do
    fimo --verbosity 1 --thresh 1e-4 --text $D/$m.meme ar_$s.fa > fimo_${m}_${s}.txt 2>/dev/null
    echo "  $m $s: $(( $(wc -l < fimo_${m}_${s}.txt) - 1 )) hits"
  done
done
