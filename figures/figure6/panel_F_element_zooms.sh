#!/usr/bin/env bash
# Figure 6, panel F -- four elements side by side (400-bp flanks): accessibility (RA, C1H-1, PVALB-1, LGE-1; 0-30), the regulator motif hits of each
# element in a grid (AR motifs red, oscine gains bold, the panel-G site boxed), oscine and avian phyloP. Tracks read from the lab track hub
# (network access); motif hits from results/mafb_domain/element_motif_gain.tsv and element_ar_sites.tsv (make_inputs.sh): zoom_panels.py.
#   bash figures/figure6/panel_F_element_zooms.sh
source "$(dirname "$0")/common.sh"
(cd "$LOCI" && "$PY" zoom_panels.py)
for e in pdf svg png; do put "$RESULTS/ucsc/zoom_panels_pvalb1.$e" "panel_F_element_zooms.$e"; done
