#!/usr/bin/env Rscript
# Bars of KCNC1 expression (SCT log-normalized) in RA, C1H-1, PVALB-1 and LAMP5 with one dot per bird: the unit is a pseudobulk of one bird in one library
# (sample_short x assignment) in one cell type, mean of the cells' values, at least MIN_CELLS cells; the bar is the mean of those per-bird values, as a standalone SVG to sit beside the four
# accessibility tracks of Fig. 4D. Written in millimeters (the figure's user unit), one row per track at the measured track-label pitch
# (fig4_scenicplus_v8.svg: label centers 15.77 px = 4.172 mm apart), rows ordered RA, C1H-1, PVALB-1, LAMP5 so it can be dropped in beside the tracks and
# aligned by eye. Text is live Arial at 5 pt (ticks) / 6 pt (title). Also writes the values (csv) with percent of cells expressing.
#
# PV_COMPARATOR=LGE1 replaces the fourth row (LAMP5) with GABA-LGE-1 and writes to <config37>/PVALB1_vs_LGE1/; the default (LAMP5) is unchanged.
#
#   Rscript kcnc1_expression_bars.R
#   PV_COMPARATOR=LGE1 Rscript kcnc1_expression_bars.R
suppressMessages({library(Seurat); library(qs2); library(tidyverse); library(here)})
OBJ <- "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/config37/")
GENE <- "KCNC1"
LGE <- identical(Sys.getenv("PV_COMPARATOR", "LAMP5"), "LGE1")
if (LGE) OUT <- file.path(OUT, "PVALB1_vs_LGE1")
ROWS <- tribble(~cluster, ~label, ~color,
                "Glut-CACNA1H-RA", "RA", "#e4332f", "Glut-CACNA1H-1", "C1H-1", "#f4a9a8", "GABA-MGE-PVALB-1", "PVALB-1", "#2a8a3a", if (LGE) "GABA-LGE-1" else "GABA-MGE-LAMP5", if (LGE) "LGE-1" else "LAMP5", "#2a8a3a")
PITCH <- 4.172; BAR_H <- 2.5; MM_PER_UNIT <- 8 / 3; MIN_CELLS <- 10; DOT_R <- 0.33; PT <- 25.4 / 72     # mm per point
o <- qs_read(OBJ, nthreads = 8); DefaultAssay(o) <- "SCT"; stopifnot(GENE %in% rownames(o))
d <- FetchData(o, vars = c(GENE, "cluster_hybrid", "sample_short", "assignment"), layer = "data") %>% rename(expr = 1) %>% filter(cluster_hybrid %in% ROWS$cluster)
bird <- d %>% group_by(cluster = cluster_hybrid, sample_short, assignment) %>% summarize(n_cells = n(), pseudobulk = mean(expr), .groups = "drop") %>% filter(n_cells >= MIN_CELLS)
res <- ROWS %>% left_join(bird %>% group_by(cluster) %>% summarize(n_birds = n(), mean_expr = mean(pseudobulk), sd_birds = sd(pseudobulk), .groups = "drop"), by = "cluster") %>%
  left_join(d %>% group_by(cluster = cluster_hybrid) %>% summarize(n_cells = n(), mean_expr_cells = mean(expr), pct_expressing = 100 * mean(expr > 0), .groups = "drop"), by = "cluster")
print(as.data.frame(res %>% select(-color)))
write_csv(res %>% select(-color) %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(OUT, "kcnc1_expression_bars.csv"))
write_csv(bird %>% mutate(pseudobulk = round(pseudobulk, 3)), file.path(OUT, "kcnc1_expression_bars_per_bird.csv"))
XMAX <- ceiling(max(c(res$mean_expr, bird$pseudobulk)))
BAR_MAX_MM <- XMAX * MM_PER_UNIT
x0 <- 0.4                                   # axis origin (mm from left)
y0 <- 7.6                                   # center of the first row (room above for the title)
W <- x0 + BAR_MAX_MM + 5; H <- y0 + 3 * PITCH + BAR_H / 2 + 5.4
yc <- y0 + (0:3) * PITCH
bars <- sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="%s"/>', x0, yc - BAR_H / 2, res$mean_expr * MM_PER_UNIT, BAR_H, res$color)
set.seed(1)
bird <- bird %>% mutate(row = match(cluster, ROWS$cluster), yd = yc[row] + runif(n(), -0.45, 0.45) * BAR_H)
dots <- sprintf('<circle cx="%.3f" cy="%.3f" r="%.2f" fill="#ffffff" fill-opacity="0.85" stroke="#222" stroke-width="0.13"/>', x0 + bird$pseudobulk * MM_PER_UNIT, bird$yd, DOT_R)
ya <- yc[4] + BAR_H / 2 + 0.5               # axis line below the last bar
ticks <- 0:XMAX; xt <- x0 + ticks * MM_PER_UNIT
axis <- c(sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', x0, ya, x0 + BAR_MAX_MM, ya),
          sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', x0, yc[1] - BAR_H / 2 - 0.5, x0, ya),
          sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', xt, ya, xt, ya + 0.5),
          sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" text-anchor="middle" fill="#222">%d</text>', xt, ya + 0.5 + 5 * PT * 0.9, 5 * PT, ticks))
title <- c(sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" font-style="italic" fill="#222">%s</text>', x0, 2.0, 6 * PT, GENE),
           sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" fill="#222">expression</text>', x0, 4.2, 6 * PT))
make_svg <- function(with_points) c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%.3fmm" height="%.3fmm" viewBox="0 0 %.3f %.3f">', W, H, W, H),
                                    '<g id="kcnc1_expression_bars">', title, bars, if (with_points) dots, axis, '</g>', '</svg>')
writeLines(make_svg(TRUE), file.path(OUT, "kcnc1_expression_bars.svg"))
writeLines(make_svg(FALSE), file.path(OUT, "kcnc1_expression_bars_no_points.svg"))     # same bars and axis, no per-bird dots
cat(sprintf("saved %s (%.1f x %.1f mm; row pitch %.3f mm; axis 0-%d = %.0f mm)\n", file.path(OUT, "kcnc1_expression_bars.svg"), W, H, PITCH, XMAX, BAR_MAX_MM))
