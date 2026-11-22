#!/usr/bin/env bash
# Figure 6, panel C -- the MAFB domain (DHX35-TOP1): RA + C1H-1 peak-to-gene links to MAFB (1-Mb run, r >= 0.45), accessibility in RA, C1H-1,
# PVALB-1 and LGE-1 (each scaled from its median to its maximum; TOP1 promoter clipped), gene models; the zoom-ins of panel F as numbered grey bars
# and a 100-kb scale bar. Rendered by the UCSC Genome Browser from the lab track hub (needs network access, an API key in ~/.ucsc_api_key and
# Inkscape on the PATH of a login shell), then tidied: ucsc_dar_batch.py overview zooms.
#   bash figures/figure6/panel_C_mafb_domain_overview.sh
source "$(dirname "$0")/common.sh"
(cd "$LOCI" && "$PY" ucsc_dar_batch.py overview zooms)
for e in pdf svg; do put "$RESULTS/ucsc/mafb_domain_overview_zooms_pvalb1.$e" "panel_C_mafb_domain_overview.$e"; done
