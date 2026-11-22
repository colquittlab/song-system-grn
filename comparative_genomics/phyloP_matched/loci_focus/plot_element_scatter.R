## Whole-tree vs oscine-lineage constraint per accessible element, RA, G4 and both side by side.
##
## Reads results/mafb_domain/element_clusters_k4.tsv (element_profiles.py, plot_element_heatmaps.R):
## one point per lobe of a SnapATAC2 peak overlapping an own-contrast ArchR DAR in the DHX35-TOP1
## domain; x = mean whole-tree phyloP, y = mean oscine-lineage (Anc222) phyloP, both over the
## accessibility summit +/-100 bp. Panels share fixed axis limits. Points are coloured by the
## element's peak-to-gene correlation with MAFB: the highest r among its own class's ArchR p2g links (either
## class's for "both" elements, RA and G4 peaks overlapping each other)
## (hybrid labels; ra links for RA, gaba links for G4; r >= 0.2 file) whose linked peak, resized to
## 500 bp about its centre, overlaps the element's SnapATAC peak. Grey = no MAFB link at r >= 0.2.
##
##   Rscript plot_element_scatter.R [mean|pct1|block|ext] [depth=T]
##   depth=T keeps lobes with every 50-bp window aligned in >= T of oscine and non-oscine genomes
##        (element_profiles.py); axis limits stay those of all lobes; outputs carry _depth<T>
##   ext  the rule-class scores over each lobe's full extent, accessibility-weighted (x whole-tree
##        block %, y oscine % > 1), with the element_rule_classes.py thresholds dashed
##   mean   mean phyloP; pct1  % of bases with phyloP > 1; block  % of bases in 11-bp blocks with
##   mean phyloP > 1 (continuity). Panels share axis limits within a measure.
## Writes results/mafb_domain/fig_element_scatter_whole_vs_oscine[_<measure>]_p2g.{pdf,png} and element_p2g_mafb.tsv.
suppressMessages({library(ggplot2); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
args <- commandArgs(trailingOnly = TRUE)
ms <- if (length(args) && !grepl("=", args[1])) args[1] else "mean"
dmin <- sub("^depth=", "", grep("^depth=", args, value = TRUE)); dmin <- if (length(dmin)) as.numeric(dmin) else NA
MS <- list(mean = c("whole_mean", "osc_mean", "mean phyloP", ""), pct1 = c("whole_pct1", "osc_pct1", "% bases phyloP > 1", "_pct1"),
           block = c("whole_block", "osc_block", "% bases in 11-bp blocks > 1", "_block"),
           ext = c("whole_ext_block", "osc_ext_pct1", "", "_ext"))[[ms]]
d <- fread(file.path(res, "element_clusters_k4.tsv"))
d <- merge(d, fread(file.path(res, "element_rule_classes.tsv"))[, .(element, osc_depth_min50, nonosc_depth_min50)], by = "element")
d[, `:=`(xv = get(MS[1]), yv = get(MS[2]))]
## peak-to-gene correlation with MAFB, per element, from its own class's link set
p2g_dir <- "/private/home/bcolquit/group/song-system-grn/multiome/archr/peaks2geneslinks_hybrid"
links <- rbindlist(lapply(c(RA = "ra", G4 = "gaba"), function(x) {
  l <- fread(file.path(p2g_dir, sprintf("peak2geneLinks_%s_corCutOff-0.2.bed", x)), header = FALSE)[V12 == "MAFB"]
  m <- regmatches(l$V17, regexec("^(.+)_([0-9]+)-([0-9]+)$", l$V17))
  mid <- (as.numeric(sapply(m, `[`, 3)) + as.numeric(sapply(m, `[`, 4))) %/% 2
  data.table(chrom = sapply(m, `[`, 2), s = mid - 250, e = mid + 250, r = l$V6)
}), idcol = "cls")
d[, p2g_r := mapply(function(k, c, a, b) { x <- links[(cls == k | k == "both") & chrom == c & s < b & e > a, r]; if (length(x)) max(x) else NA_real_ },
                    cls, chrom, peak_start, peak_end)]
d[, cls := factor(cls, levels = c("RA", "G4", "both"), labels = c("Glut-RA DARs", "GABA-4 DARs", "shared DARs"))]
pad <- function(r) r + c(-1, 1) * 0.04 * diff(r)
lim_x <- pad(range(d$xv)); lim_y <- pad(range(d$yv))
if (ms != "mean") { lim_x[1] <- min(lim_x[1], 0); lim_y[1] <- min(lim_y[1], 0) }
if (!is.na(dmin)) {                                  # axis limits stay those of all lobes
  n0 <- nrow(d); d <- d[osc_depth_min50 >= dmin & nonosc_depth_min50 >= dmin]
  message(sprintf("depth >= %g across the lobe (every 50-bp window, both clades) keeps %d of %d lobes", dmin, nrow(d), n0))
}
lab <- d[, .(n = .N, rho = cor(xv, yv, method = "spearman")), by = cls]
## n and rho go in the panel titles, so no point can sit under them
ttl <- setNames(sprintf("%s  (n = %d, rho = %.2f)", lab$cls, lab$n, lab$rho), as.character(lab$cls))

p <- ggplot(d[order(!is.na(p2g_r), p2g_r)], aes(xv, yv, colour = p2g_r)) +
  geom_hline(yintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = FIG_AXIS, linewidth = 0.3) +
  {if (ms == "ext") list(geom_vline(xintercept = d$thr_whole[1], colour = FIG_INK_SECONDARY, linewidth = 0.3, linetype = 2),
                         geom_hline(yintercept = d$thr_oscine[1], colour = FIG_INK_SECONDARY, linewidth = 0.3, linetype = 2))} +
  geom_point(size = 1.0, alpha = 0.9, stroke = 0) +
  scale_colour_gradient(low = "#cfe0f3", high = "#0b3a73", limits = c(0.2, 0.9), na.value = FIG_NEUTRAL_FILL,
                        name = "MAFB p2g r", breaks = c(0.2, 0.5, 0.8)) +
  scale_x_continuous(limits = lim_x, expand = fig_expand()) +
  scale_y_continuous(limits = lim_y, expand = fig_expand()) +
  facet_wrap(~cls, nrow = 1, labeller = as_labeller(ttl)) +
  labs(x = if (ms == "ext") "whole tree: % bases in 11-bp blocks > 1 (full extent)" else sprintf("whole tree: %s (summit ±100 bp)", MS[3]),
       y = if (ms == "ext") "oscine lineage: % bases > 1 (full extent)" else sprintf("oscine lineage: %s", MS[3])) +
  theme_fig() +
  theme(strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
        panel.spacing = unit(4, "mm"), aspect.ratio = 1,
        legend.key.width = unit(2, "mm"), legend.key.height = unit(6, "mm"),
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
stem <- paste0("fig_element_scatter_whole_vs_oscine", MS[4], "_p2g", if (is.na(dmin)) "" else sprintf("_depth%g", dmin))
fig_save(p, file.path(res, stem), width = 5.8, height = 2.0)
message("wrote ", file.path(res, stem), ".pdf/.png")
print(lab[, .(cls, n, rho = round(rho, 2))])
print(d[, .(lobes = .N, linked = sum(!is.na(p2g_r)), median_r = round(median(p2g_r, na.rm = TRUE), 2)), by = cls])
if (is.na(dmin)) fwrite(d[, .(element, cls, cluster, rule, p2g_r)], file.path(res, "element_p2g_mafb.tsv"), sep = "\t")
