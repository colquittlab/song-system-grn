#!/usr/bin/env Rscript
# Figure 4, panel A -- GRN activity per transcription factor in each cell type (positive, +/+ eRegulons).
#
# Value: gene-based AUC of each eRegulon's target genes (AUCell, computed by SCENIC+), averaged over the cells of each hybrid-labeled cell type; the direct
# and extended eRegulons of a TF share one row (their AUC values are pooled), and each row is min-max scaled to 0-1 across cell types. Rows are clustered
# with Ward's method on Euclidean distances; columns are in the fixed cell-type order. Same computation and drawing as plot_auc_heatmaps() in
# grn/analysis/scenicplus_hybrid_all_template.qmd (rendered by grn/analysis/all.qmd), which wrote
#   <scenicplus output>/config37/positive_gene_auc_heatmap_names.pdf     (and _no-names.pdf)
#
# Input: the config37 SCENIC+ MuData (scplusmdata[_slim].h5mu; see common.R).
# Output: FIG4_OUT/panel_A_positive_gene_auc_heatmap_names.pdf and _no-names.pdf
#
#   Rscript figures/figure4/panel_A_grn_activity_heatmap.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(ComplexHeatmap); library(circlize)})

auc <- read_gene_auc()
auc$cluster <- factor(as.character(auc$cluster), levels = celltype_order(levels(auc$cluster)))
m <- positive_auc_long(auc) %>%
  group_by(cluster, tf) %>% summarize(value_mean = mean(value), .groups = "drop") %>%
  group_by(tf) %>% mutate(value_mean_scale = (value_mean - min(value_mean)) / (max(value_mean) - min(value_mean))) %>% ungroup() %>%
  select(-value_mean) %>% pivot_wider(names_from = tf, values_from = value_mean_scale) %>% arrange(cluster) %>%
  column_to_rownames("cluster") %>% as.matrix() %>% t()
cat(nrow(m), "TFs x", ncol(m), "cell types\n")

d <- as.matrix(dist(m, method = "euclidean")); d[is.na(d)] <- 1; diag(d) <- 0
hc <- hclust(as.dist(d), method = "ward.D")
cols <- colorRamp2(breaks = seq(min(m), max(m), length.out = 9), colors = scales::brewer_pal(palette = "Greys", direction = 1)(9))

# with TF names (the file used for the figure)
h <- 10; w <- h * ncol(m) / nrow(m)
hm <- Heatmap(m, cluster_columns = FALSE, cluster_rows = hc, show_row_names = TRUE, col = cols, name = "target gene\nenrichment",
              row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), height = unit(h, "in"), width = unit(w, "in"))
pdf(file.path(FIG4_OUT, "panel_A_positive_gene_auc_heatmap_names.pdf"), height = h + 10, width = w + 10); print(hm); invisible(dev.off())
# without names
h <- 4; w <- 3
hm <- Heatmap(m, cluster_columns = FALSE, cluster_rows = hc, show_row_names = FALSE, col = cols, name = "target gene\nenrichment",
              height = unit(h, "in"), width = unit(w, "in"))
pdf(file.path(FIG4_OUT, "panel_A_positive_gene_auc_heatmap_no-names.pdf"), height = h + 2, width = w + 2); print(hm); invisible(dev.off())
write_csv(as_tibble(m[hc$order, ], rownames = "tf") %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(FIG4_OUT, "panel_A_positive_gene_auc_scaled.csv"))
message("saved ", file.path(FIG4_OUT, "panel_A_positive_gene_auc_heatmap_names.pdf"))
