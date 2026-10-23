#!/usr/bin/env Rscript
# Figure 3, panel F -- mean chromVAR motif enrichment per TF family, song minus non-song, compared between the three song neurons.
#
# chromVAR deviation z-scores (cisbp human motifs, ArchR addDeviationsMatrix) on the glutamatergic subset (the six song / non-song clusters) of the
# multiome ArchR project. Motifs are those matching the notebook's list of TF names (`MOTIFS`); each motif's z-score is averaged over the cells of each
# cluster, motifs are pooled into families (mean of the motif means), and each song neuron's family value has its non-song partner's subtracted
# (RA - C1H-1, HVCra - D2-1, HVCx - D2-4). One scatter per pair of song neurons, with Pearson R in the top left. Same computation and drawing as the
# "Remove redundant TFs" section of multiome/archr/chromvar.qmd, which wrote
#   multiome/archr/chromvar_hybrid/points_fam_<A>_vs_<B>_norm.pdf
#
# Input: the chromVAR-annotated glutamatergic subset project that chromvar.qmd caches (proj_glut_sub; FIG3_CHROMVAR_PROJ). It is only read here. The
# deviations take ~94 min to compute and are not recomputed by this script; see README.md.
# Output: FIG3_OUT/panel_F_points_fam_<A>_vs_<B>_norm.pdf (x3), panel_F_motif_family_song_minus_nonsong.csv
#
#   Rscript figures/figure3/panel_F_motif_family_scatter.R
# ArchR attaches plyr; load it before the tidyverse (common.R) so dplyr's mutate/summarize/rename win, as in the notebook
suppressMessages({library(ArchR); library(cowplot)})
source(here::here("figures/figure3/common.R"))
source(here::here("multiome/archr/hybrid_labels.R"))
select <- dplyr::select

CHROMVAR_PROJ <- Sys.getenv("FIG3_CHROMVAR_PROJ", here::here("multiome/archr/chromvar_hybrid/proj_glut_sub"))
addArchRLogging(useLogs = FALSE)   # ArchR otherwise writes ArchRLogs/ into the working directory
addArchRVerbose(verbose = FALSE)
addArchRThreads(threads = 8)

proj_sub <- loadArchRProject(CHROMVAR_PROJ, showLogo = FALSE)
stopifnot("MotifHSMatrix" %in% getAvailableMatrices(proj_sub))
# the cache must hold the six focal clusters under the current hybrid labels (chromvar.qmd checks the same against the shared project)
cd <- getCellColData(proj_sub)
stopifnot(all(cd$cluster_hybrid %in% names(hybrid_position)),
          identical(unname(hybrid_cluster_map()[as.character(cd$cluster)]), as.character(cd$cluster_hybrid)))
print(table(cd$cluster_hybrid))

# motifs of the notebook's TF list (z-scores only)
MOTIFS <- c("SMARCC1", "FOS", "FOS2L", "JUNB", "JUN", "JUND", "BACH1", "BACH2", "JDP2", "NFE2", "NFE2L2", "BATF", "IRF4",
            "MAFK", "MAFB", "MAFF", "MAF", "MAFG", "NR3C1", "NFIC", "AR_", "CIC",
            "MEF2", "TCF7", "SIX2", "EMX", "PPAR", "RFX", "TCF", "ARID", "BHLH", "NEUROD", "OLIG", "ATOH", "NEUROG", "TAL")
marker_motifs <- getFeatures(proj_sub, select = paste(MOTIFS, collapse = "|"), useMatrix = "MotifHSMatrix")
marker_motifs <- sub("z:", "", grep("z:", marker_motifs, value = TRUE))
cat(length(marker_motifs), "motifs\n")

mat <- getMatrixFromProject(proj_sub, useMatrix = "MotifHSMatrix")
mat_z <- assay(mat[marker_motifs, ], "z")
md <- as.data.frame(getCellColData(proj_sub)) %>% rownames_to_column("barcode") %>% select(barcode, cluster_hybrid)
rm(mat); invisible(gc())

mat_z_l <- as.data.frame(as.matrix(mat_z)) %>% rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "barcode", values_to = "value") %>%
  left_join(md, by = "barcode") %>% rename(cluster = cluster_hybrid)

# per-cluster motif means -> family means (family assignment as in the notebook)
mat_z_mean <- mat_z_l %>% group_by(cluster, TF) %>% summarize(value_mean = mean(value), .groups = "drop") %>%
  mutate(family = case_when(grepl("CTCF", TF) ~ "CTCF", grepl("^FOS", TF) ~ "FOS", grepl("JUN", TF) ~ "JUN", grepl("^MAF", TF) ~ "MAF",
                            grepl("^MEF2", TF) ~ "MEF2", grepl("^NFE2", TF) ~ "NFE2", grepl("NR3C1", TF) ~ "NR3C", grepl("PPAR", TF) ~ "PPAR",
                            grepl("NFIC", TF) ~ "NFI", TRUE ~ str_extract(TF, "[A-Z]+")))
mat_z_mean_fam <- mat_z_mean %>% group_by(cluster, family) %>% summarize(value_mean = mean(value_mean), .groups = "drop")

# song minus its non-song partner
pair_of <- unlist(lapply(names(hybrid_pairs_glut), function(p) setNames(rep(p, 2), unlist(hybrid_pairs_glut[[p]]))))
rel <- mat_z_mean_fam %>% mutate(set = pair_of[cluster], song = if_else(cluster %in% hybrid_song_clusters, "song", "nonsong")) %>%
  group_by(set, family) %>% mutate(value_mean_rel = value_mean - value_mean[song == "nonsong"]) %>% ungroup()
rel_w <- rel %>% select(family, value_mean_rel, cluster) %>% pivot_wider(names_from = "cluster", values_from = "value_mean_rel")

out_csv <- rel %>% filter(song == "song") %>% select(family, cluster, value_mean_rel) %>%
  pivot_wider(names_from = cluster, values_from = value_mean_rel) %>%
  left_join(mat_z_mean %>% count(family, name = "n_motifs_x6") %>% mutate(n_motifs = n_motifs_x6 / 6) %>% select(-n_motifs_x6), by = "family") %>%
  arrange(desc(`Glut-CACNA1H-RA`))
write_csv(out_csv %>% mutate(across(where(is.double), ~ signif(.x, 4))), file.path(FIG3_OUT, "panel_F_motif_family_song_minus_nonsong.csv"))
cat(nrow(rel_w), "families\n")

# drawing (notebook: theme_bw, 6/7 pt text, cowplot::save_plot 1.5 x 1.5 in)
annotate_pearson <- function(gg, d, x, y, size_pt = 6) {
  r <- cor(d[[x]], d[[y]], use = "complete.obs")
  message(sprintf("  %s vs %s: R = %.4f (n = %d)", x, y, r, sum(complete.cases(d[, c(x, y)]))))
  gg + annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, label = sprintf("R = %.2f", r), size = size_pt / .pt,
                family = "")   # the pdf() device default (Helvetica), as in the notebook; figure_theme.R makes Arial the geom default
}
to_plot_list <- list(c("Glut-CACNA1H-RA", "Glut-DACH2-HVCra"), c("Glut-CACNA1H-RA", "Glut-DACH2-HVCx"), c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx"))
for (to_plot in to_plot_list) {
  gg <- ggplot(rel_w, aes(.data[[to_plot[1]]], .data[[to_plot[2]]])) +
    geom_hline(yintercept = 0, linetype = 1, linewidth = .25) +
    geom_vline(xintercept = 0, linetype = 1, linewidth = .25) +
    geom_point(size = 0.5) +
    theme_bw() +
    theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6))
  gg <- annotate_pearson(gg, rel_w, to_plot[1], to_plot[2])
  save_plot(file.path(FIG3_OUT, sprintf("panel_F_points_fam_%s_vs_%s_norm.pdf", to_plot[1], to_plot[2])), gg, base_height = 1.5, base_asp = 1)
}
message("saved ", FIG3_OUT, "/panel_F_points_fam_*_norm.pdf")
