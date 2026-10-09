#!/usr/bin/env Rscript
# Heatmap of pseudobulk accessibility (z across cell types) for ALL regions behind the genes that are up in both RA and PVALB-1, regardless
# of DAR significance (matrix: scenicplus_mafb_up_both_regions_accessibility.py). Rows are clustered two ways, as for the gene_auc heatmaps:
# Ward on Euclidean distance, and average linkage on correlation distance (1 - Pearson r across cell types). Columns are in the project's
# hybrid order (not clustered). Left tracks: whether the region is a DAR in RA vs C1H-1 / PVALB-1 vs LAMP5, and a focal gene (KCNC1, PVALB,
# ERBB4; a region linked to several genes takes the first of these it matches, otherwise "other").
#
#   Rscript scenicplus_mafb_region_accessibility_heatmap.R
suppressMessages({library(ComplexHeatmap); library(circlize); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
d <- read_csv(here::here("grn/analysis/scenicplus_mafb_up_both_regions_accessibility_config37.csv"), show_col_types = FALSE)
mat <- as.matrix(d %>% select(starts_with("z_")))
colnames(mat) <- sub("^z_", "", colnames(mat)); rownames(mat) <- d$region
FOCAL <- c("KCNC1", "PVALB", "ERBB4")
focal <- sapply(strsplit(d$genes, ";"), function(g) { h <- FOCAL[FOCAL %in% g]; if (length(h)) h[1] else "other" })
focal <- factor(focal, levels = c(FOCAL, "other"))
cat("regions:", nrow(mat), "| per focal gene:", paste(names(table(focal)), table(focal), collapse = "; "), "\n")

ZC <- 2.5
col_fun <- colorRamp2(c(-ZC, 0, ZC), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH))
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)
dar_col <- c(`DAR` = FIG_INK_SECONDARY, `not DAR` = FIG_GRID)
focal_col <- setNames(c(FIG_PAL[1:3], FIG_NEUTRAL_FILL), levels(focal))
ra <- factor(ifelse(d$DAR_RA == 1, "DAR", "not DAR"), levels = names(dar_col))
pv <- factor(ifelse(d$DAR_PV1 == 1, "DAR", "not DAR"), levels = names(dar_col))

row_hc <- function(m, method) {
  if (method == "euclidean") hclust(dist(m, "euclidean"), "ward.D") else hclust(as.dist(1 - cor(t(m))), "average")
}
for (method in c("euclidean", "correlation")) {
  hc <- row_hc(mat, method)
  ha <- rowAnnotation(`RA vs C1H-1` = ra, `PVALB-1 vs LAMP5` = pv, `focal gene` = focal,
                      col = list(`RA vs C1H-1` = dar_col, `PVALB-1 vs LAMP5` = dar_col, `focal gene` = focal_col),
                      annotation_name_gp = gp5, simple_anno_size = unit(0.09, "in"),
                      annotation_legend_param = list(title_gp = gp6, labels_gp = gp5, grid_height = unit(0.09, "in"), grid_width = unit(0.09, "in")))
  hm <- Heatmap(mat, name = "z", col = col_fun, cluster_rows = hc, cluster_columns = FALSE, show_row_names = FALSE,
                column_names_gp = gp5, row_dend_width = unit(0.35, "in"), left_annotation = ha,
                width = unit(2.1, "in"), height = unit(4, "in"),
                heatmap_legend_param = list(title = "Accessibility\n(z across cell types)", title_gp = gp6, labels_gp = gp5,
                                            legend_height = unit(0.7, "in"), grid_width = unit(0.09, "in"), at = c(-ZC, 0, ZC), labels = c("≤−2.5", "0", "≥2.5")))
  stem <- file.path(OUT, MAIN, paste0("mafb_up_both_regions_accessibility_heatmap_", method))
  cairo_pdf(paste0(stem, ".pdf"), width = 5.2, height = 6.4, family = FIG_FONT)
  draw(hm, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
  invisible(dev.off())
  png(paste0(stem, ".png"), width = 5.2, height = 6.4, units = "in", res = 600, type = "cairo", family = FIG_FONT)
  draw(hm, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
  invisible(dev.off())
  cat("saved", stem, "\n")
}
