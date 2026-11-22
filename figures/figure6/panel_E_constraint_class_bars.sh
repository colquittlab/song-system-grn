#!/usr/bin/env bash
# Figure 6, panel E -- share of each phyloP class among RA, PVALB-1 and shared DARs (elements of panel D), with Holm-adjusted pairwise Fisher's
# exact p: plot_element_class_bars.R fig6. Input: results/mafb_domain/element_threshold_classes.tsv (make_inputs.sh).
#   bash figures/figure6/panel_E_constraint_class_bars.sh
source "$(dirname "$0")/common.sh"
(cd "$LOCI" && "$RSCRIPT" plot_element_class_bars.R fig6)
for e in pdf png; do put "$RESULTS/mafb_domain/fig_element_class_bars_stack_depth0.8_mini_fig6.$e" "panel_E_constraint_class_bars.$e"; done
