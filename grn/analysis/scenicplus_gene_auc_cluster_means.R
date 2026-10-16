#!/usr/bin/env Rscript
# Per-cell-type mean of the gene-based regulon AUC (+/+ eRegulons, config37), the table behind positive_gene_auc_heatmap*.pdf, so it can be analyzed without
# re-reading the 360 MB h5mu. Same steps as plot_auc_heatmaps in scenicplus_hybrid_all_template.qmd: direct and extended regulons of a TF share a row (their AUC
# columns are pooled: mean), cells are averaged within each hybrid cluster. The values are the raw means (the heatmap min-max scales each row across cell
# types, reproduced in the *_scaled file). Cells without a cluster are dropped.
#
#   Rscript scenicplus_gene_auc_cluster_means.R
suppressMessages({library(rhdf5); library(tidyverse); library(here)})
MAIN <- "config37"
F <- file.path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/scenicplus", MAIN, "outs")
ereg_file <- if (file.exists(file.path(F, "scplusmdata.h5mu"))) file.path(F, "scplusmdata.h5mu") else file.path(F, "scplusmdata_slim.h5mu")
HERE <- here::here("grn/analysis")
extract_matrix <- function(obs_name) {
  cur <- h5read(ereg_file, name = obs_name)
  mat <- data.frame(t(cur[["X"]])); colnames(mat) <- cur[["var"]]$`_index`; rownames(mat) <- cur[["obs"]]$Cell
  mat
}
gene_auc <- map(c("direct", "extended"), ~ extract_matrix(paste0("/mod/", .x, "_gene_based_AUC"))) %>% bind_cols() %>% rownames_to_column("_index")
obs <- h5read(ereg_file, name = "/mod/scRNA_counts/obs")
cl <- obs$cluster
cluster <- factor(cl$codes, levels = 0:(length(cl$categories) - 1), labels = cl$categories)
grp <- tibble(`_index` = obs$`_index`, cluster = cluster)
gene_auc <- gene_auc %>% left_join(grp, by = "_index") %>% filter(!is.na(cluster))
pos <- gene_auc %>% select(cluster, matches("\\+\\/\\+"))
colnames(pos) <- sub("_(direct|extended)_\\+\\/\\+_\\([0-9]+[a-z]\\)", "", colnames(pos))
long <- pos %>% pivot_longer(-cluster, names_to = "tf") %>% group_by(cluster, tf) %>% summarize(value_mean = mean(value), n_cells = n(), .groups = "drop")
M <- long %>% select(-n_cells) %>% pivot_wider(names_from = cluster, values_from = value_mean)
write_csv(M %>% mutate(across(where(is.double), ~ signif(.x, 5))), file.path(HERE, paste0("scenicplus_gene_auc_cluster_means_", MAIN, ".csv")))
S <- M %>% pivot_longer(-tf, names_to = "cluster") %>% group_by(tf) %>% mutate(value = (value - min(value)) / (max(value) - min(value))) %>% ungroup() %>% pivot_wider(names_from = cluster, values_from = value)
write_csv(S %>% mutate(across(where(is.double), ~ signif(.x, 5))), file.path(HERE, paste0("scenicplus_gene_auc_cluster_means_scaled_", MAIN, ".csv")))
cat(nrow(M), "TFs x", ncol(M) - 1, "cell types;", n_distinct(gene_auc$`_index`), "cells\n")
cat("cell types:", paste(setdiff(colnames(M), "tf"), collapse = ", "), "\n")
