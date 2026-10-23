#!/usr/bin/env Rscript
# Figure 3, panel D -- song vs. non-song differentially accessible regions (DARs) in the three song/non-song glutamatergic pairs.
#
# Regions: the peaks called a DAR (FDR < 0.01 and |log2FC| >= 0.6) in any of the three ArchR marker tests RA vs C1H-1, HVCra vs D2-1, HVCx vs D2-4
# (getMarkerFeatures on the PeakMatrix, bias-matched Wilcoxon), n = 97,420. Value: mean accessibility of each region in each of the six cell types
# (getGroupSE on the PeakMatrix, grouped by cluster_hybrid), min-max scaled to 0-1 per row. A random 1,000 of the regions (seed 42) are drawn, split
# into 8 k-means groups (row_km = 8, 25 repeats, seed 42) and ordered within each group by Ward clustering on correlation distance. Same computation and
# drawing as the nk = 8 heatmap chunk of multiome/archr/peaks_glut.qmd (commit 112f60f), which wrote
#   multiome/archr/peaks_glut_hybrid/DAR_fdr0.01_lfc0.6_nk8_scale01.pdf
# The figure uses the raster heatmap body of that PDF (the eight k-means slices); the labels, legend and "n =" were redrawn by hand.
#
# Inputs: the marker tests of peaks_glut.qmd (peaks_nonsong_glut_se.rds; FIG3_DAR_MARKERS; recomputed into FIG3_CACHE when absent) and the shared
#         ArchR project (read only; the cell-type mean accessibility computed from it is cached in FIG3_CACHE and reused).
# Output: FIG3_OUT/panel_D_DAR_fdr0.01_lfc0.6_nk8_scale01.pdf, panel_D_DAR_rows.csv (the 1,000 drawn regions: k-means group, row order, scaled values),
#         panel_D_DAR_counts.csv (DARs per comparison and in the union)
#
#   Rscript figures/figure3/panel_D_song_nonsong_dar_heatmap.R
source(here::here("figures/figure3/common.R"))
suppressMessages({library(ArchR); library(ComplexHeatmap); library(circlize)})
# ArchR attaches plyr after the tidyverse (loaded by common.R), masking these dplyr verbs
mutate <- dplyr::mutate; count <- dplyr::count; summarize <- dplyr::summarize; arrange <- dplyr::arrange; rename <- dplyr::rename

# marker tests, the ArchR cache and project helpers (DAR_MARKERS, FIG3_CACHE, in_cache(), archr_hybrid(), song_pair_markers()): common.R

FDR_THRESH <- 0.01; LFC_THRESH <- 0.6; NK <- 8
CT_ORDER <- c("Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-DACH2-HVCra", "Glut-DACH2-1", "Glut-DACH2-HVCx", "Glut-DACH2-4")
# legend title: ComplexHeatmap's automatic name in the notebook session; the figure's legend was redrawn ("row-scaled accessibility")
LEGEND_NAME <- "matrix_2"

scale_01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))  # as in colquittlab utils R/stats.R
cor_dist <- function(x) proxy::dist(x, method = "correlation")

# ---- DARs ---------------------------------------------------------------------------------------------------------------------------------------------
markers <- song_pair_markers()
dar <- map(markers, \(se) getMarkers(se, cutOff = sprintf("FDR < %s & abs(Log2FC) >= %s", FDR_THRESH, LFC_THRESH), returnGR = TRUE)) %>%
  map(as.data.frame) %>% bind_rows(.id = "comparison") %>% mutate(region = sprintf("%s:%s-%s", seqnames, start, end))

# ---- mean accessibility per cell type (getGroupSE on the shared project; cached) ---------------------------------------------------------------------
group_mean_file <- file.path(FIG3_CACHE, "peakmatrix_groupse_cluster_hybrid.rds")
if (!file.exists(group_mean_file)) saveRDS(in_cache(getGroupSE(ArchRProj = archr_hybrid(), useMatrix = "PeakMatrix", groupBy = "cluster_hybrid")), group_mean_file)
sePeak <- readRDS(group_mean_file)
rdata <- rowData(sePeak) %>% as.data.frame() %>% mutate(region = sprintf("%s:%s-%s", seqnames, start, end))
mat_avg <- assay(sePeak); rownames(mat_avg) <- rdata$region
mat_avg <- mat_avg[, CT_ORDER]

mat_filter <- mat_avg[rownames(mat_avg) %in% dar$region, ]
counts <- dar %>% count(comparison, direction = ifelse(Log2FC > 0, "up", "down")) %>%
  bind_rows(tibble(comparison = "union", direction = "both", n = nrow(mat_filter)))
print(as.data.frame(counts))
write_csv(counts, file.path(FIG3_OUT, "panel_D_DAR_counts.csv"))

# the random-number generator of the notebook: its setup called addArchRThreads(threads = 36), which switches R to L'Ecuyer-CMRG (for forked
# workers), so the sampling and the k-means below ran under that generator, not R's default Mersenne-Twister
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
idx <- sample(nrow(mat_filter), min(nrow(mat_filter), 1000))
mat_sample <- mat_filter[idx, ]
mat_sample_scale <- t(apply(mat_sample, 1, scale_01))
colnames(mat_sample_scale) <- CT_ORDER

# ---- heatmap ------------------------------------------------------------------------------------------------------------------------------------------
pal <- scales::brewer_pal(palette = "Greys")(9)
cols <- colorRamp2(c(0, max(mat_sample_scale)), colors = pal[c(1, 9)])
width <- .7; height <- 3
# the k-means slices are computed at the first draw() (seed 42); later prints reuse them, so the PDF and the row table below agree
pdf(NULL)
set.seed(42)
hm <- Heatmap(mat_sample_scale, col = cols, name = LEGEND_NAME,
              row_km = NK, row_km_repeats = 25,
              clustering_method_rows = "ward.D", clustering_distance_rows = cor_dist,
              cluster_columns = FALSE, column_names_gp = gpar(fontsize = 6), column_names_side = "top",
              show_row_names = FALSE, width = unit(width, "in"), height = unit(height, "in"),
              use_raster = TRUE, raster_quality = 5, border = TRUE)
hm <- draw(hm)
invisible(dev.off())
out_pdf <- file.path(FIG3_OUT, sprintf("panel_D_DAR_fdr%s_lfc%s_nk%s_scale01.pdf", FDR_THRESH, LFC_THRESH, NK))
pdf(out_pdf, height = height + 2, width = width + 2)
print(hm)
invisible(dev.off())

# ---- row membership -----------------------------------------------------------------------------------------------------------------------------------
ro <- row_order(hm)
rows <- imap_dfr(ro, \(i, k) tibble(slice = k, row_in_slice = seq_along(i), region = rownames(mat_sample_scale)[i])) %>%
  mutate(slice_position = match(slice, names(ro)), across(slice, as.integer)) %>%
  left_join(as_tibble(mat_sample_scale, rownames = "region") %>% mutate(across(where(is.double), ~ signif(.x, 4))), by = "region")
print(rows %>% count(slice_position, slice))
write_csv(rows, file.path(FIG3_OUT, "panel_D_DAR_rows.csv"))
message("saved ", out_pdf)
