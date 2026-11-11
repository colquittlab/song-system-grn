#!/bin/bash
#SBATCH -p medium
#SBATCH -t 10:00:00
#SBATCH --mem 48G
#SBATCH -c 4
#SBATCH -J rescore
#SBATCH -o /private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/rescore-%j.out
set -eu
BW=/private/groups/colquittlab/opt/ucsc_utilities/bigWigAverageOverBed
H=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/hub_upload; F=/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus
cd $F
for t in wholetree subtree222; do
  for b in ar_hits motif2_hits motif4_hits arpos jaspos ardar_mafb ardar_rest gw_motifs gw_nonmotif; do
    [ -s $b.bed ] || continue
    cut -f1-4 $b.bed > .tmp_$b.bed
    $BW $H/${t}_perbase_lsd2.bw .tmp_$b.bed ${t}_${b}.tab
    echo "  $t $b: $(wc -l < ${t}_${b}.tab)"
  done
  for tf in Fosl2 Mafb Nfia Nr2e1 Rfx4; do
    [ -s emx2hits_$tf.bed ] && $BW $H/${t}_perbase_lsd2.bw emx2hits_$tf.bed ${t}_e_$tf.tab
  done
  cp ${t}_ar_hits.tab ${t}_ar.tab
  cp ${t}_motif2_hits.tab ${t}_motif2.tab
  cp ${t}_motif4_hits.tab ${t}_motif4.tab
  cp ${t}_ardar_mafb.tab ${t}_ardar.tab
  cp ${t}_ardar_rest.tab ${t}_ardarrest.tab
  cp ${t}_gw_motifs.tab ${t}_gwm.tab
  cp ${t}_gw_nonmotif.tab ${t}_gwnm.tab
done
rm -f .tmp_*.bed
echo DONE
