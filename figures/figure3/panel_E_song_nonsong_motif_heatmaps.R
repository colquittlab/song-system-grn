#!/usr/bin/env Rscript
# Figure 3, panel E -- motifs enriched in the DARs more accessible in the song neurons ("up") and in their non-song neighbors ("down").
#
# For each pair (ra = RA vs C1H-1, hvc_1 = HVCra vs D2-1, hvc_2 = HVCx vs D2-4), ArchR motif enrichment (peakAnnoEnrichment, hypergeometric, MotifHS
# annotation = the CIS-BP human motifs of the project) in the DARs with FDR <= 0.01 and log2FC >= 0.6 (up) or <= -0.6 (down), against all peaks. The 25
# motifs named by an Ensembl gene ID rather than a TF are dropped. Rows: the union over the three pairs of each pair's top 25 motifs by -log10 adjusted
# p (ties kept); value: fold enrichment (sign flipped for the down set), rows clustered by ComplexHeatmap's defaults (Euclidean, complete linkage);
# viridis "C" (up) and reversed viridis "E" (down) scales between the matrix minimum and maximum. Same computation and drawing as the "Plot up and down
# separately" chunk of multiome/archr/peaks_glut.qmd (commit 112f60f), which wrote
#   multiome/archr/peaks_glut_hybrid/song-glut_motifs_up_top25_heatmap.pdf and song-glut_motifs_down_top25_heatmap.pdf
#
# Inputs: the marker tests of peaks_glut.qmd (peaks_nonsong_glut_se.rds; FIG3_DAR_MARKERS; recomputed into FIG3_CACHE when absent), the shared ArchR
#         project's MotifHS peak annotation (read only; the enrichment results are cached in FIG3_CACHE and reused) and data/lambert_tfs/tfs.csv.
# Output: FIG3_OUT/panel_E_song-glut_motifs_{up,down}_top25_heatmap.pdf, panel_E_song-glut_motifs_top25.csv (the drawn values, with each motif's
#         -log10 adjusted p and rank per pair)
#
#   Rscript figures/figure3/panel_E_song_nonsong_motif_heatmaps.R
source(here::here("figures/figure3/common.R"))
suppressMessages({library(ArchR); library(ComplexHeatmap); library(circlize)})
# ArchR attaches plyr after the tidyverse (loaded by common.R), masking these dplyr verbs
mutate <- dplyr::mutate; count <- dplyr::count; summarize <- dplyr::summarize; arrange <- dplyr::arrange; rename <- dplyr::rename

# marker tests, the ArchR cache and project helpers (DAR_MARKERS, FIG3_CACHE, in_cache(), archr_hybrid(), song_pair_markers()): common.R

TFS <- here::here("data/lambert_tfs/tfs.csv")

TOP_N <- 25
# legend titles: ComplexHeatmap's automatic names in the notebook session; the figure's legends were redrawn ("motif enrichment")
LEGEND_NAME <- c(up = "matrix_6", down = "matrix_7")

# ---- motif enrichment in the up and down DARs of each pair (cached) ----------------------------------------------------------------------------------
enrich_file <- file.path(FIG3_CACHE, "motifhs_enrichment_song_glut_dars.rds")
if (!file.exists(enrich_file)) {
  markers <- song_pair_markers(); proj <- archr_hybrid()
  enr <- list(up = "FDR <= 0.01 & Log2FC >= .6", down = "FDR <= 0.01 & Log2FC <= -0.6") %>%
    map(\(cut) map(markers, \(se) in_cache(peakAnnoEnrichment(seMarker = se, ArchRProj = proj, peakAnnotation = "MotifHS", cutOff = cut))))
  saveRDS(enr, enrich_file); rm(proj, markers); invisible(gc())
}
enr <- readRDS(enrich_file)

drop_unnamed_motifs <- function(df) df[!grepl("^ENSG", df$TF), , drop = FALSE]
to_df <- function(lst) map(lst, function(x) {
  df <- data.frame(TF = rownames(x), mlog10Padj = assay(x)[, 1], enrichment = assay(x, "Enrichment")[, 1])
  df <- drop_unnamed_motifs(df)
  df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
  df$rank <- seq_len(nrow(df))
  df
}) %>% bind_rows(.id = "pair")
motifs_df <- bind_rows(list(up = to_df(enr$up), down = to_df(enr$down) %>% mutate(enrichment = -1 * enrichment)), .id = "direction") %>%
  mutate(Name = str_extract(TF, "[A-Z0-9]+")) %>%
  left_join(read_delim(TFS, show_col_types = FALSE), by = "Name")

# ---- heatmaps -----------------------------------------------------------------------------------------------------------------------------------------
drawn <- map(c("up", "down"), function(direction_cur) {
  motifs_df_cur <- motifs_df %>% filter(direction == direction_cur)
  motifs_df_f <- motifs_df_cur %>% group_by(pair) %>% top_n(TOP_N, mlog10Padj)
  m <- motifs_df_cur %>% filter(TF %in% motifs_df_f$TF) %>% select(pair, TF, enrichment) %>%
    pivot_wider(names_from = "pair", values_from = "enrichment", values_fill = 0) %>% column_to_rownames("TF") %>% as.matrix()
  pal <- if (direction_cur == "up") scales::viridis_pal(option = "C")(9) else scales::viridis_pal(option = "E", direction = -1)(9)
  cols <- colorRamp2(c(min(m), max(m)), colors = pal[c(1, 9)])
  width <- .35; height <- 3
  hm <- Heatmap(m, col = cols, name = LEGEND_NAME[[direction_cur]], cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6),
                width = unit(width, "in"), height = unit(height, "in"))
  out_pdf <- file.path(FIG3_OUT, sprintf("panel_E_song-glut_motifs_%s_top%s_heatmap.pdf", direction_cur, TOP_N))
  pdf(out_pdf, height = height + 2, width = width + 2); print(hm); invisible(dev.off())
  message("saved ", out_pdf)
  cat(direction_cur, ":", nrow(m), "motifs; range", signif(range(m), 4), "\n")
  pdf(NULL); ro <- row_order(draw(hm)); invisible(dev.off())
  tibble(direction = direction_cur, TF = rownames(m)[ro], row = seq_along(ro))
})

# the drawn values, top to bottom, with each motif's statistics in every pair
rows <- bind_rows(drawn) %>%
  left_join(motifs_df %>% select(direction, pair, TF, enrichment, mlog10Padj, rank), by = c("direction", "TF")) %>%
  mutate(across(c(enrichment, mlog10Padj), ~ signif(.x, 4)))
write_csv(rows, file.path(FIG3_OUT, "panel_E_song-glut_motifs_top25.csv"))
