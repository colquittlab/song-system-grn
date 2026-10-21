#!/usr/bin/env Rscript
# Figure 1, panel G -- marker expression, dissected-region composition and regional specificity of each cell type.
#
# Four parts, all in the row order CT_ORDER (common.R):
#   1. Dot plot of curated markers (scCustomize::DotPlot_scCustom on the SCT assay: colour = average scaled expression, size = % of cells
#      expressing; dots shown for >= 5 % of cells). Same as snrna/deg/dotplot_manual_markers_hybrid.qmd -> dotplot_manual-markers.pdf.
#   2. Relative fraction of each cell type in each dissected region ("rel. fraction in region"): the region x cell-type count table normalized to
#      region totals, then to cell-type totals.
#   3. snRNA-seq regional specificity: KL divergence (nats) of each cell type's normalized region profile (2) from uniform. 2 and 3 are the same as
#      snrna/trees/celltypes_hclust_all_hybrid.qmd -> correlation_average/cluster_proportions_vert_with-bar_ordered.pdf.
#   4. Spatial (Xenium) regional specificity ("spatial trans."): for each cell type and section, confident RCTD calls with >= 400 transcripts are
#      binned on a 400 um grid; the Jensen-Shannon divergence (bits) of 100 subsampled cells' bin distribution from that of all cells, minus the
#      same for 100 cells drawn at random from all cells (each the mean of 30 draws), averaged over the sections with >= 100 cells of the type.
#      Glut-DACH2-HVCra-Int has < 100 such cells in every section and is not scored (N.D.). Same as
#      xenium/label_transfer/spatial_variation_by_type.R (score) + spatial_variation_aligned.R -> spatial_variation/spatial_variation_aligned_bars.pdf,
#      which draws 2-4 side by side; the figure takes 2 and 3 from part 3's file and 4 from this one.
#
# Inputs: the snRNA-seq object with the hybrid labels; the RCTD calls, proseg centroids and transcript counts (common.R).
# Outputs: FIG1_OUT/panel_G_dotplot_manual_markers.pdf
#          FIG1_OUT/panel_G_region_fraction_snrna_specificity.pdf (+ panel_G_region_fraction.csv)
#          FIG1_OUT/panel_G_spatial_specificity_aligned_bars.{pdf,png} (+ panel_G_spatial_specificity_by_type[_section].csv)
#
#   Rscript figures/figure1/panel_G_markers_regional_specificity.R
source(here::here("figures/figure1/common.R"))
suppressMessages({library(Seurat); library(scCustomize); library(cowplot); library(egg); library(ComplexHeatmap); library(circlize)})
theme_set(theme_cowplot())

obj <- read_snrna()
stopifnot("doublet assignments present" = !any(grepl("/", obj$assignment)))
bar_colors <- unname(CLASS_COLORS[celltype_class(CT_ORDER)])

# ---- 1. marker dot plot -----------------------------------------------------------------------------------------------------------------------------
MARKERS <- c("SLC17A6", "SYT1", "DACH2", "SATB2", "LHX9", "CACNA1H",
             "NECTIN3", "DCX", "SOX4", "NR2E1",
             "GAD1", "FOXP2", "NPY", "CEMIP", "NXPH1", "TTLL5", "PENK", "CRH", "LAMP5", "CALB2", "LHX8",
             "SLC1A3", "SPEF2", "PLP1", "PDGFRA", "CSF1R", "FLI1")
stopifnot(all(MARKERS %in% rownames(obj)))
# reversed levels put the first type of CT_ORDER at the top
obj@meta.data[[CELLTYPE_COL]] <- factor(obj@meta.data[[CELLTYPE_COL]], levels = rev(CT_ORDER))
height <- length(CT_ORDER) * 0.08; width <- length(MARKERS) * 0.08
gg <- DotPlot_scCustom(obj, features = MARKERS, group.by = CELLTYPE_COL, flip_axes = FALSE, dot.scale = 2, dot.min = 0.05, col.min = 0) +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.text.y = element_text(size = 6), axis.line = element_blank(), axis.ticks = element_blank())
gg <- egg::set_panel_size(p = gg, width = unit(width, "in"), height = unit(height, "in"))
save_plot(file.path(FIG1_OUT, "panel_G_dotplot_manual_markers.pdf"), gg, base_height = height + 5, base_width = width + 5)
message("saved ", file.path(FIG1_OUT, "panel_G_dotplot_manual_markers.pdf"))

# ---- 2, 3. relative fraction per region and snRNA-seq regional specificity --------------------------------------------------------------------------
tab <- table(obj$position, as.character(obj@meta.data[[CELLTYPE_COL]]))[POSITIONS, CT_ORDER]
rm(obj); invisible(gc())
tab <- sweep(tab, 1, rowSums(tab), FUN = "/")
tab <- sweep(tab, 2, colSums(tab), FUN = "/")
tab_spec <- apply(tab, 2, calc_tissue_specificity)
as_tibble(t(unclass(tab)), rownames = "celltype") %>% mutate(snrna_specificity = tab_spec[celltype]) %>%
  write_csv(file.path(FIG1_OUT, "panel_G_region_fraction.csv"))

cols <- colorRamp2(breaks = seq(0, 1, length.out = 9), colors = scales::brewer_pal(palette = "Greys")(9))
# 0.08 in per row and column, as in celltypes_hclust_all_hybrid.qmd (the heatmap body is 0.48 in wide, 3.76 in tall)
height <- nrow(tab) * 0.08; width <- ncol(tab) * 0.08
hm <- Heatmap(t(tab), cluster_rows = FALSE, show_row_names = TRUE, row_names_side = "left", cluster_columns = FALSE, col = cols,
              border = TRUE, border_gp = gpar(lwd = 0.5), height = unit(width, "in"), width = unit(height, "in"),
              column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6))
ann <- HeatmapAnnotation(ab = anno_barplot(tab_spec, 0, which = "row", border = TRUE, gp = gpar(fill = bar_colors),
                                           axis_param = list(gp = gpar(fontsize = 5))), which = "row")
pdf(file.path(FIG1_OUT, "panel_G_region_fraction_snrna_specificity.pdf"), height = width + 4, width = height + 4)
draw(hm + ann)
invisible(dev.off())
message("saved ", file.path(FIG1_OUT, "panel_G_region_fraction_snrna_specificity.pdf"))

# ---- 4. spatial (Xenium) regional specificity -------------------------------------------------------------------------------------------------------
BIN_UM <- 400; N_FIX <- 100; B <- 30
set.seed(1)
xd <- read_xenium_calls() %>%
  filter(spot_class %in% CONFIDENT_CLASSES, nCount >= MIN_TRANSCRIPTS) %>%
  mutate(bin = paste(floor(x_centroid / BIN_UM), floor(y_centroid / BIN_UM)))
jsd <- function(p, q) {
  m <- (p + q) / 2
  kl <- function(a, b) { i <- a > 0; sum(a[i] * log2(a[i] / b[i])) }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}
by_section <- map_dfr(unique(xd$section_id), function(sid) {
  sec <- xd %>% filter(section_id == sid)
  bins <- sort(unique(sec$bin))
  qv <- table(factor(sec$bin, levels = bins)); q <- as.numeric(qv) / sum(qv)
  tys <- sec %>% count(first_type) %>% filter(n >= N_FIX) %>% pull(first_type)
  map_dfr(tys, function(ty) {
    tb <- sec$bin[sec$first_type == ty]
    o <- replicate(B, { s <- sample(tb, N_FIX); jsd(as.numeric(table(factor(s, levels = bins))) / N_FIX, q) })
    nl <- replicate(B, { s <- sample(sec$bin, N_FIX); jsd(as.numeric(table(factor(s, levels = bins))) / N_FIX, q) })
    tibble(section_id = sid, first_type = ty, n_cells = length(tb), obs = mean(o), null = mean(nl), delta = mean(o) - mean(nl))
  })
})
by_type <- by_section %>% group_by(first_type) %>%
  summarize(n_sections = n(), delta_sd = sd(delta), delta = mean(delta), median_n = median(n_cells), .groups = "drop") %>% arrange(desc(delta))
write_csv(by_section, file.path(FIG1_OUT, "panel_G_spatial_specificity_by_type_section.csv"))
write_csv(by_type, file.path(FIG1_OUT, "panel_G_spatial_specificity_by_type.csv"))
delta <- by_type$delta[match(CT_ORDER, by_type$first_type)]   # NA (empty slot) for types not scored
message("not scored: ", paste(CT_ORDER[is.na(delta)], collapse = ", "))

hm_h <- 3.765                                       # body height of the heatmap in part 3's file
hm <- Heatmap(t(tab), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, row_names_side = "left", col = cols,
              border = TRUE, border_gp = gpar(lwd = 0.5), height = unit(hm_h, "in"), width = unit(height, "in"),
              column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title = "prop", labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6)))
bar_ann <- function(name, v) {
  bars <- anno_barplot(v, 0, which = "row", border = FALSE, gp = gpar(fill = bar_colors), width = unit(1, "cm"), axis_param = list(gp = gpar(fontsize = 5)))
  do.call(HeatmapAnnotation, c(setNames(list(bars), name), list(which = "row", annotation_name_gp = gpar(fontsize = 6), gap = unit(2, "mm"))))
}
ht <- hm + bar_ann("snRNA\nspecificity", tab_spec) + bar_ann("Xenium\nlocalisation", delta)
out <- file.path(FIG1_OUT, "panel_G_spatial_specificity_aligned_bars")
pdf(paste0(out, ".pdf"), height = hm_h + 2, width = 5.2); draw(ht, merge_legend = TRUE); invisible(dev.off())
png(paste0(out, ".png"), height = (hm_h + 2) * 300, width = 5.2 * 300, res = 300); draw(ht, merge_legend = TRUE); invisible(dev.off())
message("saved ", out)
