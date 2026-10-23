#!/usr/bin/env Rscript
# Figure 3, panel C -- SIX2 expression beside the SIX2-locus genome-browser tracks. The tracks themselves are a genome-browser export (pseudobulk ATAC
# coverage per cell type); see README.md.
#
# One bar per track row (RA, C1H-1, HVCra, D2-1, HVCx, D2-4): SCT log-normalized expression (natural log, log1p of corrected counts), averaged per bird
# (one bird in one library, >= MIN_CELLS cells) and then across birds; optional dots = the per-bird values. Same construction as the KCNC1 bars of
# Figure 4E. Written as an SVG in millimeters with one row per track at the track pitch of the figure (2.56 mm), so it can be aligned beside the tracks.
#
# Inputs: the hybrid-labeled Seurat object (common.R).
# Output: FIG3_OUT/panel_C_six2_expression_bars_no_points.svg (the version in the figure), _bars.svg (+ per-bird dots), _labeled.svg (row and gene
#         labels, for use away from the tracks), + csv of the values
#
#   Rscript figures/figure3/panel_C_six2_expression_bars.R
source(here::here("figures/figure3/common.R"))
suppressMessages({library(Seurat); library(qs2)})

GENE <- "SIX2"; MIN_CELLS <- 10
PITCH <- 2.5598; BAR_H <- 1.6; BAR_MAX_MM <- 9; DOT_R <- 0.25; PT <- 25.4 / 72; LABEL_W <- 6.2

o <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8)
ex <- gene_expression_by_bird(o, GENE, SONG_PAIRS$cluster, MIN_CELLS); rm(o); invisible(gc())
res <- SONG_PAIRS %>% left_join(ex$summary, by = "cluster")
bird <- ex$per_bird
print(as.data.frame(res %>% select(label, n_birds, mean_expr, sd_birds, n_cells, pct_expressing)))
write_csv(res %>% select(-color) %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(FIG3_OUT, "panel_C_six2_expression_bars.csv"))

# The bars stand directly to the right of the track rows, so in the figure they carry no row labels (the track labels name the rows) and the x-axis title is
# "log expr." (one line under the ticks). The _labeled version names each row and the gene, for use away from the tracks.
y0 <- BAR_H / 2 + 0.5
yc <- y0 + (0:5) * PITCH
ya <- yc[6] + BAR_H / 2 + 0.4
txt <- function(x, y, s, pt, anchor = "start", extra = "") sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" text-anchor="%s" fill="#222"%s>%s</text>', x, y, pt * PT, anchor, extra, s)
set.seed(1)
bird <- bird %>% mutate(row = match(cluster, SONG_PAIRS$cluster), yd = yc[row] + runif(n(), -0.45, 0.45) * BAR_H)
# x scale: BAR_MAX_MM spans the largest value drawn (the bars, plus the per-bird dots when shown); at most two tick intervals so the 5 pt labels fit
make_svg <- function(with_points, labeled = FALSE) {
  x0 <- if (labeled) LABEL_W else 0.4
  W <- x0 + BAR_MAX_MM + if (labeled) 4.5 else 1.5; H <- ya + 5.4
  labels <- if (labeled) txt(x0 - 0.6, yc + 5 * PT * 0.35, res$label, 5, "end")
  title <- sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" text-anchor="middle" fill="#222">%s</text>', x0 + BAR_MAX_MM / 2,
                   ya + 0.4 + 5 * PT * 0.85 + 6 * PT * 1.05, 6 * PT, if (labeled) sprintf('<tspan font-style="italic">%s</tspan> log expr.', GENE) else "log expr.")
  vmax <- max(c(res$mean_expr, if (with_points) bird$pseudobulk))
  steps <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1, 2, 5)
  step <- steps[which(vmax / steps <= 2)[1]]
  xmax <- step * ceiling(vmax / step); mm <- BAR_MAX_MM / xmax
  ticks <- seq(0, xmax, by = step); xt <- x0 + ticks * mm
  bars <- sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="%s"/>', x0, yc - BAR_H / 2, res$mean_expr * mm, BAR_H, res$color)
  dots <- sprintf('<circle cx="%.3f" cy="%.3f" r="%.2f" fill="#ffffff" fill-opacity="0.85" stroke="#222" stroke-width="0.1"/>', x0 + bird$pseudobulk * mm, bird$yd, DOT_R)
  axis <- c(sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.1"/>', x0, ya, x0 + BAR_MAX_MM, ya),
            sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.1"/>', x0, yc[1] - BAR_H / 2 - 0.4, x0, ya),
            sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.1"/>', xt, ya, xt, ya + 0.4),
            txt(xt, ya + 0.4 + 5 * PT * 0.85, format(ticks, drop0trailing = TRUE, trim = TRUE), 5, "middle"))
  c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%.3fmm" height="%.3fmm" viewBox="0 0 %.3f %.3f">', W, H, W, H),
    '<g id="six2_expression_bars">', title, labels, bars, if (with_points) dots, axis, '</g>', '</svg>')
}
writeLines(make_svg(TRUE), file.path(FIG3_OUT, "panel_C_six2_expression_bars.svg"))
writeLines(make_svg(FALSE), file.path(FIG3_OUT, "panel_C_six2_expression_bars_no_points.svg"))
writeLines(make_svg(FALSE, labeled = TRUE), file.path(FIG3_OUT, "panel_C_six2_expression_bars_labeled.svg"))
message("saved ", file.path(FIG3_OUT, "panel_C_six2_expression_bars_no_points.svg"))
