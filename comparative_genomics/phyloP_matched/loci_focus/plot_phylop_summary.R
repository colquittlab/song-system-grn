## Constraint at high-accessibility-difference bases, Glut-RA vs GABA-4 activity, by span.
##
## Reads results/compare_mafb_spans/phylop_summary_q<Q>.tsv (phylop_summary_table.py) and draws,
## per track, the percentage of high-difference bases with phyloP > Q (top) and the mean phyloP
## of those bases (bottom), for the DHX35-TOP1 domain, MAFB +/-1.5 Mb and genome-wide. Points are
## pooled over each class's bases with 95% peak-cluster bootstrap intervals; filled = all
## elements, open = Log2FC-matched.
##
##   Rscript plot_phylop_summary.R [Q=1] [window="500 bp"]
## Writes results/compare_mafb_spans/fig_phylop_summary_{3model,perbase}_<window>_q<Q>.{pdf,png}.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
## the repo's renv may be unavailable off its home machine; find the theme relative to this script
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
.theme <- if (length(.self)) {
  normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R"))
} else here::here("config/figure_theme.R")
source(.theme)

args <- commandArgs(trailingOnly = TRUE)
q <- if (length(args) >= 1) args[1] else "1"
win <- if (length(args) >= 2) args[2] else "500 bp"
res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/compare_mafb_spans"

d <- read.delim(file.path(res, sprintf("phylop_summary_q%s.tsv", q)), check.names = FALSE) %>%
  filter(window == win) %>%
  mutate(span = factor(span, levels = c("DHX35–TOP1 domain", "MAFB ±1.5 Mb", "genome-wide"),
                       labels = c("domain", "±1.5 Mb", "genome-wide")),
         cls = factor(cls, levels = c("RA", "G4"), labels = c("Glut-RA", "GABA-4")),
         set = factor(set, levels = c("all", "Log2FC-matched")))

long <- bind_rows(
  d %>% transmute(tracks, track, set, span, cls, value = pct, lo = pct_lo, hi = pct_hi,
                  measure = sprintf("%% bases phyloP > %s", q)),
  d %>% transmute(tracks, track, set, span, cls, value = mean, lo = mean_lo, hi = mean_hi,
                  measure = sprintf("mean phyloP, bases > %s", q))
) %>% mutate(measure = factor(measure, levels = unique(measure)))

cls_col <- c("Glut-RA" = FIG_PAL[[2]], "GABA-4" = FIG_PAL[[1]])
dodge <- position_dodge(width = 0.75)

for (ts in unique(long$tracks)) {
  x <- long %>% filter(tracks == ts) %>% mutate(grp = interaction(cls, set))
  trk <- unique(x$track)
  x$track <- factor(x$track, levels = trk,
                    labels = sub(" \\(matched\\)", "", sub(", birdAnc222", " (Anc222)", trk)))
  p <- ggplot(x, aes(span, value, colour = cls, group = grp)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.35, position = dodge) +
    geom_point(aes(shape = set), size = 1.1, stroke = 0.45, fill = FIG_SURFACE, position = dodge) +
    scale_colour_manual(values = cls_col, name = NULL) +
    scale_shape_manual(values = c("all" = 16, "Log2FC-matched" = 21), name = NULL) +
    scale_y_continuous(expand = fig_expand(0.06)) +
    facet_grid(measure ~ track, scales = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_fig() +
    theme(strip.placement = "outside",
          strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
          legend.position = "bottom", legend.key.size = unit(6, "pt"),
          legend.text = element_text(size = FIG_PT_AXIS_TEXT),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "pt"),
          axis.text.x = element_text(angle = 30, hjust = 1))
  stem <- sprintf("fig_phylop_summary_%s_%s_q%s", ifelse(ts == "3-model", "3model", "perbase"),
                  gsub(" ", "", win), q)
  fig_save(p, file.path(res, stem), width = 1.15 * length(trk) + 0.5, height = 2.6)
  message("wrote ", file.path(res, stem), ".pdf/.png")
}
