#!/usr/bin/env Rscript
# Figure 4, panel E -- the parts of the KCNC1 locus panel that are made by code. The accessibility tracks themselves are a genome-browser export
# (pseudobulk ATAC coverage of each cell type at the KCNC1 locus); see README.md.
#
# 1. MAFB cistrome track (BED, for the genome browser): the regions of the MAFB +/+ eRegulon (config37), one line per region, named by the linked target
#    gene(s); score = 250 + 750 * sqrt(region-to-gene importance / max), i.e. compressed so that weak links stay visible (UCSC useScore shading).
#    Coordinates are lonStrDom2 (UCSC chromosome names), 0-based BED.
# 2. KCNC1 expression next to the tracks, one bar per track row (RA, C1H-1, PVALB-1, LGE-1): SCT log-normalized expression (natural log, log1p of corrected
#    counts), averaged per bird (one bird in one library, >= MIN_CELLS cells) and then across birds; optional dots = the per-bird values. Written as an SVG in
#    millimeters with one row per track at the track-label pitch of the figure (4.172 mm), so it can be aligned beside the tracks.
# Same output as the one-off cistrome export (MAFB_cistrome_regions_config37.bed) and grn/analysis/kcnc1_expression_bars.R with PV_COMPARATOR=LGE1
# (<scenicplus output>/config37/PVALB1_vs_LGE1/kcnc1_expression_bars_no_points.svg).
#
# Inputs: the config37 SCENIC+ MuData and the hybrid-labeled Seurat object (common.R).
# Output: FIG4_OUT/panel_E_MAFB_cistrome_regions.bed, panel_E_kcnc1_expression_bars[_no_points].svg (+ csv of the values)
#
#   Rscript figures/figure4/panel_E_kcnc1_track_annotations.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(Seurat); library(qs2)})

## 1. MAFB cistrome BED
e <- read_eregulons() %>% filter(TF == "MAFB")
cat("MAFB eRegulons:", paste(unique(e$eRegulon_name), collapse = ", "), "\n")
# the importance is first expressed on UCSC's 0-1000 score scale (rounded), then compressed, as in the BED used for the figure
bed <- e %>% mutate(importance_R2G = as.numeric(importance_R2G)) %>% group_by(Region) %>%
  summarize(genes = paste(sort(unique(Gene), method = "radix"), collapse = "/"), imp = max(importance_R2G), .groups = "drop") %>%
  extract(Region, c("chrom", "start", "end"), "^(.+):(\\d+)-(\\d+)$", convert = TRUE) %>%
  mutate(score1000 = round(1000 * imp / max(imp)), score = round(250 + 750 * sqrt(score1000 / 1000)), name = paste0("MAFB:", genes), strand = ".") %>% arrange(chrom, start)
f <- file.path(FIG4_OUT, "panel_E_MAFB_cistrome_regions.bed")
writeLines('track name="MAFB_cistrome" description="MAFB eRegulon regions (config37, +/+); name = linked genes; score = 250 + 750*sqrt(max R2G importance / max), compressed so weak links stay visible" useScore=1', f)
write_tsv(bed %>% select(chrom, start, end, name, score, strand), f, append = TRUE, col_names = FALSE)
cat(nrow(bed), "regions ->", f, "\n")

## 2. KCNC1 expression bars
GENE <- "KCNC1"; MIN_CELLS <- 10
ROWS <- tribble(~cluster, ~label, ~color, RA, "RA", "#e4332f", C1H1, "C1H-1", "#f4a9a8", PV1, "PVALB-1", "#2a8a3a", LGE1, "LGE-1", "#2a8a3a")
PITCH <- 4.172; BAR_H <- 2.5; MM_PER_UNIT <- 8 / 3; DOT_R <- 0.33; PT <- 25.4 / 72
o <- qs_read(SEURAT_HYBRID_OBJ, nthreads = 8); DefaultAssay(o) <- "SCT"; stopifnot(GENE %in% rownames(o))
d <- FetchData(o, vars = c(GENE, "cluster_hybrid", "sample_short", "assignment"), layer = "data") %>% rename(expr = 1) %>% filter(cluster_hybrid %in% ROWS$cluster)
bird <- d %>% group_by(cluster = cluster_hybrid, sample_short, assignment) %>% summarize(n_cells = n(), pseudobulk = mean(expr), .groups = "drop") %>% filter(n_cells >= MIN_CELLS)
res <- ROWS %>% left_join(bird %>% group_by(cluster) %>% summarize(n_birds = n(), mean_expr = mean(pseudobulk), sd_birds = sd(pseudobulk), .groups = "drop"), by = "cluster") %>%
  left_join(d %>% group_by(cluster = cluster_hybrid) %>% summarize(n_cells = n(), mean_expr_cells = mean(expr), pct_expressing = 100 * mean(expr > 0), .groups = "drop"), by = "cluster")
print(as.data.frame(res %>% select(label, n_birds, mean_expr, sd_birds, n_cells)))
write_csv(res %>% select(-color) %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(FIG4_OUT, "panel_E_kcnc1_expression_bars.csv"))
XMAX <- ceiling(max(c(res$mean_expr, bird$pseudobulk))); BAR_MAX_MM <- XMAX * MM_PER_UNIT
x0 <- 0.4; y0 <- 7.6
W <- x0 + BAR_MAX_MM + 5; H <- y0 + 3 * PITCH + BAR_H / 2 + 5.4
yc <- y0 + (0:3) * PITCH
bars <- sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="%s"/>', x0, yc - BAR_H / 2, res$mean_expr * MM_PER_UNIT, BAR_H, res$color)
set.seed(1)
bird <- bird %>% mutate(row = match(cluster, ROWS$cluster), yd = yc[row] + runif(n(), -0.45, 0.45) * BAR_H)
dots <- sprintf('<circle cx="%.3f" cy="%.3f" r="%.2f" fill="#ffffff" fill-opacity="0.85" stroke="#222" stroke-width="0.13"/>', x0 + bird$pseudobulk * MM_PER_UNIT, bird$yd, DOT_R)
ya <- yc[4] + BAR_H / 2 + 0.5
ticks <- 0:XMAX; xt <- x0 + ticks * MM_PER_UNIT
axis <- c(sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', x0, ya, x0 + BAR_MAX_MM, ya),
          sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', x0, yc[1] - BAR_H / 2 - 0.5, x0, ya),
          sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#222" stroke-width="0.12"/>', xt, ya, xt, ya + 0.5),
          sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" text-anchor="middle" fill="#222">%d</text>', xt, ya + 0.5 + 5 * PT * 0.9, 5 * PT, ticks))
title <- c(sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" font-style="italic" fill="#222">%s</text>', x0, 2.0, 6 * PT, GENE),
           sprintf('<text x="%.3f" y="%.3f" font-family="Arial" font-size="%.3f" fill="#222">expression</text>', x0, 4.2, 6 * PT))
make_svg <- function(with_points) c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%.3fmm" height="%.3fmm" viewBox="0 0 %.3f %.3f">', W, H, W, H),
                                    '<g id="kcnc1_expression_bars">', title, bars, if (with_points) dots, axis, '</g>', '</svg>')
writeLines(make_svg(TRUE), file.path(FIG4_OUT, "panel_E_kcnc1_expression_bars.svg"))
writeLines(make_svg(FALSE), file.path(FIG4_OUT, "panel_E_kcnc1_expression_bars_no_points.svg"))
message("saved ", file.path(FIG4_OUT, "panel_E_kcnc1_expression_bars_no_points.svg"))
