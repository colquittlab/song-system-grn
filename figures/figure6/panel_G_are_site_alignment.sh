#!/usr/bin/env bash
# Figure 6, panel G -- the ARE-full site of element 3 (chr20:12,946,574-12,946,585) across the 363-bird alignment: pruned cladogram with
# reconstructed presence, one representative genome per clade (Climacteris rufus shown for its pair: it carries the site), motif score and genomes
# with the site; 143 mm wide, compact layout: maf_switch/plot_site_alignment.py. Input: the lobe alignment extract (make_inputs.sh).
#   bash figures/figure6/panel_G_are_site_alignment.sh
source "$(dirname "$0")/common.sh"
d=$(mktemp -d)
(cd "$MAF_DATA" && "$PY_ALN" "$LOCI/maf_switch/plot_site_alignment.py" "chr20:12,946,375-12,946,665" full --aln 363 \
   --coords chr20:12946574-12946585 --cls shared --maf data/lobes_gap_pvalb1.maf \
   --paper 143 --row-mm 2.18 --flank 8 --compact --show Ptilonorhynchus_violaceus=Climacteris_rufus \
   --outdir "$d" --out site)
for e in pdf svg png; do put "$d/site.$e" "panel_G_are_site_alignment.$e"; done
rm -rf "$d"
