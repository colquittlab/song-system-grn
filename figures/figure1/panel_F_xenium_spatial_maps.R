#!/usr/bin/env Rscript
# Figure 1, panel F -- Xenium positions of GABA-MGE-PVALB-1/2/3 and Astro-1/2/3 cells.
#
# Cells labelled by RCTD label transfer from the snRNA-seq hybrid labels (confident calls: singlet or confident doublet; >= 400 transcripts) drawn
# in black over every profiled cell of the section in grey, all ten sections in a common orientation (posterior left, dorsal up). Same as the
# "every cell call in the taxonomy (gte400)" section of xenium/label_transfer/proseg_spatial_maps_hybrid.qmd, which wrote
#   xenium/label_transfer/proseg_spatial_maps_hybrid/spatial_all_<type>_gte400.png
# The figure shows one section of each of these ten-section maps, cropped and relabelled by hand: OR52YW26_1_4 (1.2 mm from the midline) for the
# PVALB types and OR52YW26_1_7 (2.2 mm) for the astrocyte types.
#
# Input: the RCTD calls, proseg centroids and transcript counts (common.R).
# Output: FIG1_OUT/panel_F_spatial_<type>_gte400.png, one per type
#
#   Rscript figures/figure1/panel_F_xenium_spatial_maps.R
source(here::here("figures/figure1/common.R"))
suppressMessages(library(cowplot))
TYPES <- c("GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2", "GABA-MGE-PVALB-3", "Astro-1", "Astro-2", "Astro-3")

d <- read_xenium_calls() %>%
  group_by(section_id) %>% group_modify(~ orient_xy(.x, .y$section_id, "x_centroid", "y_centroid")) %>% ungroup() %>%
  mutate(section_id = factor(section_id, levels = SECTIONS), confident = spot_class %in% CONFIDENT_CLASSES)
stopifnot(all(TYPES %in% d$first_type))

theme_spatial <- theme_cowplot() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5))
section_plot <- function(sid, hl) {
  p <- ggplot() +
    geom_point(data = d %>% filter(section_id == sid), aes(x_centroid, y_centroid), colour = "grey85", size = 0.03, alpha = 0.4) +
    coord_equal() + labs(title = sid) + theme_spatial
  hl <- hl %>% filter(section_id == sid)
  if (nrow(hl) == 0) return(p)
  p + geom_point(data = hl, aes(x_centroid, y_centroid), colour = "#000000", size = 0.3, alpha = 0.85)
}

for (ty in TYPES) {
  hl <- d %>% filter(confident, nCount >= MIN_TRANSCRIPTS, first_type == ty)
  p <- plot_grid(plotlist = map(SECTIONS, section_plot, hl = hl), ncol = 5) +
    draw_figure_label(paste0(ty, "  (n = ", nrow(hl), " confident cells, nCount >= ", MIN_TRANSCRIPTS, ")"), position = "top", size = 14)
  out <- file.path(FIG1_OUT, paste0("panel_F_spatial_", str_replace_all(ty, "[/-]", "_"), "_gte400.png"))
  ggsave(out, p, width = 20, height = 9.5, dpi = 250, bg = "white")
  message("saved ", out)
}
