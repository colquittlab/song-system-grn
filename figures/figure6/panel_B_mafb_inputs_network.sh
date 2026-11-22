#!/usr/bin/env bash
# Figure 6, panel B -- predicted SCENIC+ inputs to MAFB in C1H-RA (vs C1H-1) and GABAergic neurons, on one circle.
# grn/analysis/mafb_inputs_network_circle.py (config2 eRegulon tables of the glut-ra_glut-c1h-1_allruns_hybrid and ra-arco-hvc-nc_gaba_hybrid runs,
# under /hdd/jupyter on lark: run this panel there). KLF15, weakly expressed in C1H-RA, is left out.
#   bash figures/figure6/panel_B_mafb_inputs_network.sh
source "$(dirname "$0")/common.sh"
d=$(mktemp -d)
MAFB_NET_EXCLUDE=${MAFB_NET_EXCLUDE:-RA_all:KLF15} "${FIG6_PY_SCENIC:-python3}" "$REPO/grn/analysis/mafb_inputs_network_circle.py" "$d"
for e in svg pdf png; do put "$d/mafb_inputs_network_circle.$e" "panel_B_mafb_inputs_network.$e"; done
rm -rf "$d"
