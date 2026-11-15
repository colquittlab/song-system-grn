## Whole-tree vs oscine-lineage constraint per accessible element, RA and G4 side by side.
##
## Reads results/mafb_domain/element_clusters_k4.tsv (element_profiles.py, plot_element_heatmaps.R):
## one point per lobe of a SnapATAC2 peak overlapping an own-contrast ArchR DAR in the DHX35-TOP1
## domain; x = mean whole-tree phyloP, y = mean oscine-lineage (Anc222) phyloP, both over the
## accessibility summit +/-100 bp. Panels share fixed axis limits.
##
##   Rscript plot_element_scatter.R
## Writes results/mafb_domain/fig_element_scatter_whole_vs_oscine.{pdf,png}.
suppressMessages({library(ggplot2); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
d <- fread(file.path(res, "element_clusters_k4.tsv"))
d[, cls := factor(cls, levels = c("RA", "G4"), labels = c("Glut-RA DARs", "GABA-4 DARs"))]
lim_x <- range(d$whole_mean) + c(-0.15, 0.15); lim_y <- range(d$osc_mean) + c(-0.1, 0.1)
lab <- d[, .(n = .N, rho = cor(whole_mean, osc_mean, method = "spearman")), by = cls]
lab[, txt := sprintf("n = %d\nrho = %.2f", n, rho)]

p <- ggplot(d, aes(whole_mean, osc_mean, colour = cls)) +
  geom_hline(yintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.8, stroke = 0) +
  geom_text(data = lab, aes(x = lim_x[2], y = lim_y[2], label = txt), hjust = 1, vjust = 1,
            size = fig_pt(5), colour = FIG_INK_SECONDARY, lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(values = c("Glut-RA DARs" = FIG_PAL[[2]], "GABA-4 DARs" = FIG_PAL[[1]]), guide = "none") +
  scale_x_continuous(limits = lim_x, expand = fig_expand()) +
  scale_y_continuous(limits = lim_y, expand = fig_expand()) +
  facet_wrap(~cls, nrow = 1) +
  labs(x = "whole-tree phyloP (mean, summit ±100 bp)", y = "oscine-lineage phyloP (mean)") +
  theme_fig() +
  theme(strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
        panel.spacing = unit(4, "mm"), aspect.ratio = 1)
fig_save(p, file.path(res, "fig_element_scatter_whole_vs_oscine"), width = 3.6, height = 2.0)
message("wrote ", file.path(res, "fig_element_scatter_whole_vs_oscine"), ".pdf/.png")
print(lab[, .(cls, n, rho = round(rho, 2))])
