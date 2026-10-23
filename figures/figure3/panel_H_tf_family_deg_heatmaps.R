#!/usr/bin/env Rscript
# Figure 3, panel H -- expression log2 fold-change, song vs non-song, of the members of four TF families (MAF, TCF, BHLH, AP-1).
#
# avg_log2FC of each song neuron against its non-song partner (columns RA vs C1H-1, HVCra vs D2-1, HVCx vs D2-4) from the pseudobulk DESeq2 table of
# multiome/seurat/deg/song-surround_deg_glut.qmd; genes absent from a contrast's table count as 0. One heatmap per family; rows clustered
# (ComplexHeatmap defaults: complete linkage, Euclidean distance), columns fixed. The color scale is shared: PRGn, green at the minimum and purple at the
# maximum log2FC over all TF-family genes (the panel G genes), white at 0. Same computation and drawing as the heatmap_deg loop of
# multiome/archr/chromvar.qmd, which wrote
#   multiome/archr/chromvar_hybrid/heatmap_deg_<family>.pdf     (the notebook writes all 22 families; the figure uses these four)
#
# Inputs: the DEG table (FIG3_DEG_TABLE) and the gene names of the hybrid-labeled multiome object the notebook reads (FIG3_HYBRID_LABEL_OBJ).
# Output: FIG3_OUT/panel_H_heatmap_deg_{MAF,TCF,BHLH,AP1}.pdf, panel_H_tf_family_log2fc.csv (rows in the drawn order)
#
#   Rscript figures/figure3/panel_H_tf_family_deg_heatmaps.R
source(here::here("figures/figure3/common.R"))
suppressMessages({library(ComplexHeatmap); library(circlize)})

FAMILIES <- c("MAF", "TCF", "BHLH", "AP1")   # the four in the figure, top to bottom
# legend titles: ComplexHeatmap's automatic names in the notebook session; the figure replaced the legends with one color bar
LEGEND_NAME <- c(MAF = "matrix_16", TCF = "matrix_25", BHLH = "matrix_20", AP1 = "matrix_13")

# DEG_TABLE, HYBRID_LABEL_OBJ, TF_TERMS, read_song_deg(), hybrid_features(): common.R
deg <- read_song_deg()
features <- hybrid_features()

lfc_matrix <- function(genes) deg %>% filter(gene %in% genes) %>% select(gene, compare, avg_log2FC) %>% filter(grepl("RA|HVC", compare)) %>%
  pivot_wider(names_from = "compare", values_from = "avg_log2FC", values_fill = 0) %>% column_to_rownames("gene") %>% as.matrix()

# shared color scale: range over every TF-family gene
all_mat <- lfc_matrix(unlist(map(TF_TERMS, ~ sort(grep(paste(.x, collapse = "|"), features, value = TRUE)))))
pal <- scales::brewer_pal(type = "div", palette = 3)(9)
cols <- colorRamp2(c(min(all_mat), 0, max(all_mat)), colors = pal[c(9, 5, 1)])
cat(sprintf("color scale: %.3f (green) .. 0 .. %.3f (purple), %d genes\n", min(all_mat), max(all_mat), nrow(all_mat)))

drawn <- map_dfr(FAMILIES, function(fam) {
  genes <- sort(grep(paste(TF_TERMS[[fam]], collapse = "|"), features, value = TRUE))
  m <- lfc_matrix(genes)[, c("Glut-CACNA1H-RA", "Glut-DACH2-HVCra", "Glut-DACH2-HVCx")]
  height <- nrow(m) * 0.1; width <- ncol(m) * 0.1
  hm <- Heatmap(m, col = cols, name = LEGEND_NAME[[fam]], cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
                height = unit(height, "in"), width = unit(width, "in"))
  pdf(file.path(FIG3_OUT, sprintf("panel_H_heatmap_deg_%s.pdf", fam)), height = height + 2, width = width + 2)
  hm <- draw(hm); invisible(dev.off())
  ord <- row_order(hm)
  as_tibble(m[ord, , drop = FALSE], rownames = "gene") %>% mutate(family = fam, .before = 1)
})
print(as.data.frame(drawn %>% mutate(across(where(is.double), ~ round(.x, 2)))))
write_csv(drawn %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(FIG3_OUT, "panel_H_tf_family_log2fc.csv"))
message("saved ", FIG3_OUT, "/panel_H_heatmap_deg_*.pdf")
