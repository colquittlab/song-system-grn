#!/usr/bin/env Rscript
# Heatmap of cell-type specificity (z-score of the cluster mean across cell types) for MAFB's +/+ target genes in the main
# comparison config (config37). Rows are targets (Ward-clustered on their z profile), columns are the cell types in the project's
# hybrid order. Gene labels are shown for the genes discussed in the RA / PVALB / SST analysis; the left strip marks the targets
# found in at least half of the sweep configs (the stable core). The z matrix is saved as a small csv next to this script.
#
#   Rscript scenicplus_mafb_target_specificity_heatmap.R
suppressMessages({library(Seurat); library(qs2); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))   # hybrid_ct_order
fig_check_font()

MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
HERE <- here::here("grn/analysis")
LABEL <- c("PVALB", "KCNC1", "LOC110481261", "MAFB", "PCSK5", "PIK3R6", "ABHD12", "FAM184B", "RUNX2", "PIK3R5", "ADAMTS18", "MTNR1B",
           "ST6GALNAC5", "PLXNA4", "TENM4", "KIAA1217", "CNTNAP2", "VIPR2", "ERBB4", "ARX")   # LOC110481261 = KCNC3 (fs_gene_aliases_lonStrDom2.tsv)
ZCAP <- 3

obj <- qs_read(OBJ, nthreads = 8)
obj$cluster <- obj$cluster_hybrid
n <- table(obj$cluster)
avg <- log1p(AverageExpression(obj, assays = "SCT", layer = "data", group.by = "cluster")[[1]])
rm(obj); invisible(gc())
cl <- intersect(hybrid_ct_order, names(n)[n >= 50])
avg <- avg[, cl]
sdv <- apply(avg, 1, sd)
z <- (avg - rowMeans(avg)) / sdv

e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name))
targets <- intersect(unique(e$Gene), rownames(z)[!is.na(sdv) & sdv > 0])
cat(MAIN, "MAFB +/+ targets:", length(unique(e$Gene)), "| with expression variation:", length(targets), "| clusters:", length(cl), "\n")
zt <- z[targets, , drop = FALSE]
write_csv(as_tibble(round(zt, 3), rownames = "gene"), file.path(HERE, paste0("scenicplus_mafb_target_specificity_", MAIN, ".csv")))

ord <- rownames(zt)[hclust(dist(zt), method = "ward.D2")$order]
core <- read_csv(file.path(HERE, "scenicplus_mafb_target_membership.csv"), show_col_types = FALSE) %>% filter(n_pp >= 12) %>% pull(Gene)
long <- as_tibble(zt, rownames = "gene") %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "z") %>%
  mutate(gene = factor(gene, levels = rev(ord)), cluster = factor(cluster, levels = cl), zc = pmax(pmin(z, ZCAP), -ZCAP))
lab <- intersect(LABEL, targets)
cat("labeled:", paste(lab, collapse = ", "), "\n")
cat("not MAFB +/+ targets in", MAIN, ":", paste(setdiff(LABEL, targets), collapse = ", "), "\n")

strip <- tibble(gene = factor(ord, levels = rev(ord)), core = ord %in% core)
nx <- length(cl)
lab_df <- tibble(gene = factor(lab, levels = rev(ord)), x = nx + 0.8)
p_main <- ggplot(long, aes(cluster, gene, fill = zc)) +
  geom_raster() +
  ggrepel::geom_text_repel(data = lab_df, aes(x = x, y = gene, label = gene), inherit.aes = FALSE, hjust = 0, direction = "y",
                           xlim = c(nx + 0.8, NA), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                           segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0, box.padding = 0.12,
                           point.padding = 0, max.overlaps = Inf, seed = 1) +
  scale_fill_gradient2(low = FIG_DIV_LOW, mid = FIG_DIV_MID, high = FIG_DIV_HIGH, midpoint = 0, limits = c(-ZCAP, ZCAP),
                       name = "Specificity (z across cell types)", breaks = c(-3, 0, 3), labels = c("\u2264\u22123", "0", "\u22653")) +
  scale_x_discrete(expand = expansion(add = c(0, 5.5))) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FIG_PT_AXIS_TEXT),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        legend.position = "bottom", legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.9, "cm"),
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
p_strip <- ggplot(strip, aes(1, gene, fill = core)) +
  geom_raster() +
  scale_fill_manual(values = c(`TRUE` = FIG_INK_SECONDARY, `FALSE` = FIG_GRID), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  theme_void()
p <- (p_strip | p_main) + patchwork::plot_layout(widths = c(0.03, 1)) +
  patchwork::plot_annotation(caption = "Left bar: target present in at least 12 of 24 sweep configs (stable core)",
                             theme = theme(plot.caption = element_text(size = FIG_PT_AXIS_TEXT, family = FIG_FONT, hjust = 0, color = FIG_INK_SECONDARY)))
out <- file.path(OUT, MAIN, "mafb_target_specificity_heatmap")
fig_save(p, out, width = 4, height = 5.2)
cat("saved", out, "\n")
