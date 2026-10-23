#!/usr/bin/env Rscript
# Figure 3, panel G -- expression log2 fold-change (song vs non-song) of the TFs of the panel F families, compared between the three song neurons.
#
# avg_log2FC of each song neuron against its non-song partner (RA vs C1H-1, HVCra vs D2-1, HVCx vs D2-4) from the pseudobulk DESeq2 table of
# multiome/seurat/deg/song-surround_deg_glut.qmd, for the genes matching the TF-family patterns below (genes absent from a contrast's table count as 0).
# One scatter per pair of song neurons, with Pearson R in the top left. Same computation and drawing as the points_deg chunk of
# multiome/archr/chromvar.qmd, which wrote
#   multiome/archr/chromvar_hybrid/points_deg_<A>_vs_<B>.pdf
#
# Inputs: the DEG table (FIG3_DEG_TABLE) and the gene names of the hybrid-labeled multiome object the notebook reads (FIG3_HYBRID_LABEL_OBJ).
# Output: FIG3_OUT/panel_G_points_deg_<A>_vs_<B>.pdf (x3), panel_G_tf_log2fc.csv
#
#   Rscript figures/figure3/panel_G_tf_deg_scatter.R
source(here::here("figures/figure3/common.R"))
suppressMessages(library(cowplot))

# DEG_TABLE, HYBRID_LABEL_OBJ, TF_TERMS, read_song_deg(), hybrid_features(): common.R
deg <- read_song_deg()
features <- hybrid_features()
genes_to_use <- unlist(map(TF_TERMS, ~ sort(grep(paste(.x, collapse = "|"), features, value = TRUE))))

deg_filt_mat <- deg %>% filter(gene %in% genes_to_use) %>% select(gene, compare, avg_log2FC) %>%
  filter(grepl("RA|HVC", compare)) %>%
  pivot_wider(names_from = "compare", values_from = "avg_log2FC", values_fill = 0)
write_csv(deg_filt_mat %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(FIG3_OUT, "panel_G_tf_log2fc.csv"))
cat(nrow(deg_filt_mat), "TF genes\n")

annotate_pearson <- function(gg, d, x, y, size_pt = 6) {
  r <- cor(d[[x]], d[[y]], use = "complete.obs")
  message(sprintf("  %s vs %s: R = %.4f (n = %d)", x, y, r, nrow(d)))
  gg + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, label = sprintf("R = %.2f", r), size = size_pt / .pt,
                family = "")   # the pdf() device default (Helvetica), as in the notebook; figure_theme.R makes Arial the geom default
}
to_plot_list <- list(c("Glut-CACNA1H-RA", "Glut-DACH2-HVCra"), c("Glut-CACNA1H-RA", "Glut-DACH2-HVCx"), c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx"))
for (to_plot in to_plot_list) {
  gg <- ggplot(deg_filt_mat, aes(.data[[to_plot[1]]], .data[[to_plot[2]]])) +
    geom_hline(yintercept = 0, linetype = 1, linewidth = .25) +
    geom_vline(xintercept = 0, linetype = 1, linewidth = .25) +
    geom_point(size = 0.5) +
    theme_bw() +
    theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6))
  gg <- annotate_pearson(gg, deg_filt_mat, to_plot[1], to_plot[2])
  save_plot(file.path(FIG3_OUT, sprintf("panel_G_points_deg_%s_vs_%s.pdf", to_plot[1], to_plot[2])), gg, base_height = 1.5, base_asp = 1)
}
message("saved ", FIG3_OUT, "/panel_G_points_deg_*.pdf")
