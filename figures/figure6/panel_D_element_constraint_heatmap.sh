#!/usr/bin/env bash
# Figure 6, panel D -- accessibility-weighted mean avian and oscine phyloP of the MAFB-domain DARs (RA, PVALB-1, shared) aligned in >= 80% of
# oscine and non-oscine genomes (n = 66), grouped by phyloP class and Ward-clustered within each: plot_element_cluster.R fig6 depth=0.8.
# Input: results/mafb_domain/element_rule_classes.tsv and element_threshold_classes.tsv (make_inputs.sh).
#   bash figures/figure6/panel_D_element_constraint_heatmap.sh
source "$(dirname "$0")/common.sh"
(cd "$LOCI" && "$RSCRIPT" plot_element_cluster.R fig6 depth=0.8)
for e in pdf png; do put "$RESULTS/mafb_domain/fig_element_cluster_thr_noacc_depth0.8_nolabels_fig6.$e" "panel_D_element_constraint_heatmap.$e"; done
