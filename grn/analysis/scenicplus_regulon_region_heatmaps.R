#!/usr/bin/env Rscript
# One heatmap per regulon of the regions behind its targets: how accessible each region is in RA, C1H-1, PVALB-1 and LAMP5, and the linear DIFFERENCE
# in accessibility for the two contrasts (RA - C1H-1 and PVALB-1 - LAMP5), not the log ratio. Accessibility is pseudobulk CPM (fragments per million
# fragments in the cell type) from the 40-topic cisTopic object (scenicplus_region_accessibility_all_regions.py); the difference is in CPM, so a region
# that goes from 1 to 2 CPM and one that goes from 10 to 20 CPM are not scored alike (a ratio would treat them as equal).
# Regulons are the config37 +/+ regulons of the cross-regulon comparisons (scenicplus_accessibility_concordance_across_regulons.R):
#   whole_regulon/   regions linked to all of the TF's +/+ targets (TFs with >= MIN_REGIONS_ALL regions)
#   up_in_both/      regions linked to the TF's targets with RNA log2FC > UP_LFC in both contrasts (>= MIN_GENES_UP genes and >= MIN_REGIONS_UP regions)
# Rows are clustered (Ward, Euclidean) on the two differences; the accessibility block follows the same order. In the whole-regulon heatmaps the left
# bar marks regions linked to an up-in-both gene. Output: <out>/regulon_region_heatmaps/<view>/<TF>.pdf and an index csv.
#
#   Rscript scenicplus_regulon_region_heatmaps.R
suppressMessages({library(ComplexHeatmap); library(circlize); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
UP_LFC <- 1; MIN_REGIONS_ALL <- 100; MIN_GENES_UP <- 5; MIN_REGIONS_UP <- 30; MAX_CLUSTER <- 12000
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
D <- file.path(OUT, MAIN, "regulon_region_heatmaps"); for (v in c("whole_regulon", "up_in_both")) dir.create(file.path(D, v), showWarnings = FALSE, recursive = TRUE)
HERE <- here::here("grn/analysis")

acc <- read_tsv(file.path(D, "region_accessibility_cpm.tsv.gz"), col_types = "cdddd", progress = FALSE)
A <- as.matrix(acc[, c("RA", "C1H1", "PV1", "LAMP5")]); rownames(A) <- acc$region
dif <- cbind(RA = A[, "RA"] - A[, "C1H1"], PV1 = A[, "PV1"] - A[, "LAMP5"])
DC <- as.numeric(quantile(abs(dif), 0.99)); AC <- as.numeric(quantile(log2(A + 1), 0.99))     # shared color caps (99th percentile over all regions)
cat(sprintf("color caps: difference +/-%.1f CPM, accessibility log2(CPM+1) %.1f\n", DC, AC))
col_diff <- colorRamp2(c(-DC, 0, DC), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH))
col_acc <- colorRamp2(c(0, AC), c("#e6eef9", "#0f3c78"))
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)

e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% distinct(TF, Gene, Region) %>% filter(Region %in% rownames(A))
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE)
tested <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low)]
up_both <- g$gene[!is.na(g$lfc_RA_vs_C1H1) & !is.na(g$lfc_PV1_vs_low) & g$lfc_RA_vs_C1H1 > UP_LFC & g$lfc_PV1_vs_low > UP_LFC]

draw_regulon <- function(tf, view, regions, flag = NULL) {
  n <- length(regions)
  if (n > MAX_CLUSTER) regions <- regions[order(dif[regions, 1])]     # too many rows to cluster: sort by the RA difference
  dd <- dif[regions, , drop = FALSE]; aa <- log2(A[regions, , drop = FALSE] + 1)
  ord <- if (n <= MAX_CLUSTER) hclust(dist(dd), "ward.D") else NULL
  rho <- suppressWarnings(cor(dd[, 1], dd[, 2], method = "spearman"))
  colnames(dd) <- c("RA − C1H-1", "PVALB-1 − LAMP5"); colnames(aa) <- c("RA", "C1H-1", "PVALB-1", "LAMP5")
  ra <- if (!is.null(flag)) rowAnnotation(`up-in-both gene` = factor(ifelse(flag, "yes", "no"), levels = c("yes", "no")), col = list(`up-in-both gene` = c(yes = FIG_INK_SECONDARY, no = FIG_GRID)),
                                         annotation_name_gp = gp5, simple_anno_size = unit(0.09, "in"), show_legend = FALSE) else NULL
  h_acc <- Heatmap(aa, name = "accessibility", col = col_acc, cluster_rows = if (!is.null(ord)) ord else FALSE, cluster_columns = FALSE, show_row_names = FALSE, column_names_gp = gp5,
                   show_row_dend = !is.null(ord) && n <= 2500, row_dend_width = unit(0.3, "in"), left_annotation = ra, use_raster = TRUE, raster_quality = 2,
                   width = unit(0.95, "in"), height = unit(4.2, "in"), column_title = paste0(tf, ": ", n, " regions"), column_title_gp = gp6,
                   heatmap_legend_param = list(title = "Accessibility\n(log2 CPM+1)", title_gp = gp6, labels_gp = gp5, legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in")))
  h_dif <- Heatmap(dd, name = "difference", col = col_diff, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, column_names_gp = gp5, use_raster = TRUE, raster_quality = 2,
                   width = unit(0.5, "in"), height = unit(4.2, "in"), column_title = sprintf("ρ = %.2f", rho), column_title_gp = gp5,
                   heatmap_legend_param = list(title = "Difference\n(CPM)", title_gp = gp6, labels_gp = gp5, legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in"),
                                               at = c(-DC, 0, DC), labels = c(sprintf("≤−%.0f", DC), "0", sprintf("≥%.0f", DC))))
  cairo_pdf(file.path(D, view, paste0(tf, ".pdf")), width = 4.6, height = 5.6, family = FIG_FONT)
  # the accessibility block's row order comes from clustering the difference block; the difference heatmap follows it (main_heatmap)
  draw(h_acc + h_dif, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"), row_dend_side = "left", main_heatmap = "accessibility")
  invisible(dev.off())
  tibble(TF = tf, view = view, n_regions = n, rho_differences = rho, median_diff_RA = median(dd[, 1]), median_diff_PV1 = median(dd[, 2]),
         frac_up_in_both_CPM = mean(dd[, 1] > 0 & dd[, 2] > 0))
}

idx <- list()
for (tf in sort(unique(e$TF))) {
  et <- e %>% filter(TF == tf, Gene %in% tested)
  all_regs <- unique(et$Region)
  if (length(all_regs) >= MIN_REGIONS_ALL) {
    flag <- all_regs %in% unique(et$Region[et$Gene %in% up_both])
    idx[[length(idx) + 1]] <- draw_regulon(tf, "whole_regulon", all_regs, flag)
  }
  G <- intersect(unique(et$Gene), up_both); rG <- unique(et$Region[et$Gene %in% G])
  if (length(G) >= MIN_GENES_UP && length(rG) >= MIN_REGIONS_UP) idx[[length(idx) + 1]] <- draw_regulon(tf, "up_in_both", rG)
}
IDX <- bind_rows(idx)
write_csv(IDX %>% mutate(across(where(is.double), ~ round(.x, 4))), file.path(D, "index.csv"))
cat("heatmaps:", table(IDX$view)[["whole_regulon"]], "whole regulon,", table(IDX$view)[["up_in_both"]], "up in both ->", D, "\n")
