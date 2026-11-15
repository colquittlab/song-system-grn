## Whole-tree vs oscine-lineage constraint per accessible element, RA and G4 side by side.
##
## Reads results/mafb_domain/element_clusters_k4.tsv (element_profiles.py, plot_element_heatmaps.R):
## one point per lobe of a SnapATAC2 peak overlapping an own-contrast ArchR DAR in the DHX35-TOP1
## domain; x = mean whole-tree phyloP, y = mean oscine-lineage (Anc222) phyloP, both over the
## accessibility summit +/-100 bp. Panels share fixed axis limits.
##
##   Rscript plot_element_scatter.R [mean|pct1|block]
##   mean   mean phyloP; pct1  % of bases with phyloP > 1; block  % of bases in 11-bp blocks with
##   mean phyloP > 1 (continuity). Panels share axis limits within a measure.
## Writes results/mafb_domain/fig_element_scatter_whole_vs_oscine[_<measure>].{pdf,png}.
suppressMessages({library(ggplot2); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
args <- commandArgs(trailingOnly = TRUE)
ms <- if (length(args)) args[1] else "mean"
MS <- list(mean = c("whole_mean", "osc_mean", "mean phyloP", ""), pct1 = c("whole_pct1", "osc_pct1", "% bases phyloP > 1", "_pct1"),
           block = c("whole_block", "osc_block", "% bases in 11-bp blocks > 1", "_block"))[[ms]]
d <- fread(file.path(res, "element_clusters_k4.tsv"))
d[, `:=`(xv = get(MS[1]), yv = get(MS[2]))]
d[, cls := factor(cls, levels = c("RA", "G4"), labels = c("Glut-RA DARs", "GABA-4 DARs"))]
pad <- function(r) r + c(-1, 1) * 0.04 * diff(r)
lim_x <- pad(range(d$xv)); lim_y <- pad(range(d$yv))
if (ms != "mean") { lim_x[1] <- min(lim_x[1], 0); lim_y[1] <- min(lim_y[1], 0) }
lab <- d[, .(n = .N, rho = cor(xv, yv, method = "spearman")), by = cls]
## n and rho go in the panel titles, so no point can sit under them
ttl <- setNames(sprintf("%s  (n = %d, rho = %.2f)", lab$cls, lab$n, lab$rho), as.character(lab$cls))

p <- ggplot(d, aes(xv, yv, colour = cls)) +
  geom_hline(yintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = c("Glut-RA DARs" = FIG_PAL[[2]], "GABA-4 DARs" = FIG_PAL[[1]]), guide = "none") +
  scale_x_continuous(limits = lim_x, expand = fig_expand()) +
  scale_y_continuous(limits = lim_y, expand = fig_expand()) +
  facet_wrap(~cls, nrow = 1, labeller = as_labeller(ttl)) +
  labs(x = sprintf("whole tree: %s (summit ±100 bp)", MS[3]), y = sprintf("oscine lineage: %s", MS[3])) +
  theme_fig() +
  theme(strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
        panel.spacing = unit(4, "mm"), aspect.ratio = 1)
stem <- paste0("fig_element_scatter_whole_vs_oscine", MS[4])
fig_save(p, file.path(res, stem), width = 3.6, height = 2.0)
message("wrote ", file.path(res, stem), ".pdf/.png")
print(lab[, .(cls, n, rho = round(rho, 2))])
