## ---------------------------------------------------------------------------
## Class-level summary of cross-species integration strength,
## finch x chicken vs finch x mouse.
##
## Tests the write-up claim that all finch cell classes integrate well with
## chicken, while with mouse non-neuronal cells integrate strongly, GABAergic
## neurons moderately and glutamatergic neurons weakly. Reads the two full-suite
## hybrid-label composites (GSI + SAMap + CCA + SATURN):
##   composite_scoring/results/gg_adult_hybrid/composite_calls.csv
##   composite_scoring/results/yao_adult_hybrid/composite_calls.csv
##
## Two per-cluster readouts, summarised by finch cell class (Non-neuronal /
## GABAergic / Glutamatergic, from the celltype_hybrid name prefix):
##   top row     composite confidence -- the Zaremba magnitude channel at the
##               rank-aggregated winner; "is the best match strong?"
##   bottom row  method agreement -- how many of the 4 methods independently pick
##               the same reference label as their own top call.
## Points are finch clusters; bars are the class median (top) or mean (bottom).
## Within-pair class differences: one-sided Mann-Whitney U (Non-neuronal > GABA
## > Glut), normal approximation with continuity correction. Chicken-vs-mouse
## drop per class: paired one-sided Wilcoxon over shared finch clusters.
##
## Writes class_summary_table.csv / class_summary_tests.csv (tracked) and the
## figure (PDF + PNG, gitignored) to composite_scoring/results/class_summary/.
## ---------------------------------------------------------------------------

suppressMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
source(here::here("config/figure_theme.R"))

RES <- here::here("snrna/integration/composite_scoring/results")
OUT <- file.path(RES, "class_summary"); dir.create(OUT, showWarnings = FALSE)

PAIRS <- c(gg_adult_hybrid = "Finch × chicken", yao_adult_hybrid = "Finch × mouse")
CLASSES <- c("Non-neuronal", "GABAergic", "Glutamatergic")
CLASS_COL <- setNames(FIG_PAL[1:3], CLASSES)   # fixed slots, all-pairs safe
METRICS <- c(confidence = "Composite confidence\n(match strength at winner)",
             agree_frac = "Method agreement\n(fraction of 4 methods)")

cell_class <- function(q) ifelse(grepl("^GABA", q), "GABAergic",
                          ifelse(grepl("^Glut", q), "Glutamatergic", "Non-neuronal"))

D <- bind_rows(lapply(names(PAIRS), function(tag) {
  read.csv(file.path(RES, tag, "composite_calls.csv"), check.names = FALSE) %>%
    transmute(pair = PAIRS[[tag]], tag, query, confidence, confidence_tier,
              methods_agreeing_on_own_top, n_methods,
              agree_frac = methods_agreeing_on_own_top / n_methods,
              cls = factor(cell_class(query), CLASSES))
}))

## --- Tables -------------------------------------------------------------------
summ <- D %>% group_by(pair, cell_class = cls) %>%
  summarise(n_clusters = n(),
            confidence_median = median(confidence), confidence_mean = mean(confidence),
            frac_high_tier = mean(confidence_tier == "high"),
            frac_low_tier = mean(confidence_tier == "low"),
            agree_frac_mean = mean(agree_frac),
            frac_all_methods_agree = mean(methods_agreeing_on_own_top == n_methods),
            frac_le1_method_agree = mean(methods_agreeing_on_own_top <= 1), .groups = "drop")
write.csv(summ, file.path(OUT, "class_summary_table.csv"), row.names = FALSE)

mwu <- function(a, b) wilcox.test(a, b, alternative = "greater", exact = FALSE, correct = TRUE)$p.value
COMPS <- list(c("Non-neuronal", "GABAergic"), c("GABAergic", "Glutamatergic"),
              c("Non-neuronal", "Glutamatergic"))
tests <- bind_rows(
  lapply(names(PAIRS), function(tag) bind_rows(lapply(names(METRICS), function(met) bind_rows(lapply(COMPS, function(cp) {
    d <- D[D$tag == tag, ]
    tibble(pair = PAIRS[[tag]], metric = met, test = sprintf("%s > %s", cp[1], cp[2]),
           p = mwu(d[[met]][d$cls == cp[1]], d[[met]][d$cls == cp[2]]))
  }))))),
  lapply(CLASSES, function(c) bind_rows(lapply(names(METRICS), function(met) {
    G <- D[D$tag == "gg_adult_hybrid" & D$cls == c, ]; Y <- D[D$tag == "yao_adult_hybrid", ]
    idx <- G$query[G$query %in% Y$query]
    p <- tryCatch(wilcox.test(G[[met]][match(idx, G$query)], Y[[met]][match(idx, Y$query)],
                              paired = TRUE, alternative = "greater")$p.value,
                  error = function(e) NA_real_)
    tibble(pair = "chicken vs mouse (paired)", metric = met, test = sprintf("%s: chicken > mouse", c), p = p)
  })))
)
write.csv(tests, file.path(OUT, "class_summary_tests.csv"), row.names = FALSE)
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
print(as.data.frame(tests %>% mutate(p = signif(p, 3))), row.names = FALSE)

## --- Figure ---------------------------------------------------------------------
## Full p-value rather than a significance star, per direct instruction --
## scientific notation below 0.001 (where a fixed-decimal value would round to
## 0.000), three decimal places otherwise.
pfmt <- function(p) ifelse(p < 0.001, sprintf("p=%.1e", p), sprintf("p=%.3f", p))
long <- D %>% pivot_longer(c(confidence, agree_frac), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, names(METRICS)), pair = factor(pair, PAIRS), x = as.integer(cls))

## point positions: jitter for the continuous metric, deterministic fan-out in rows
## of six for the discrete agreement fractions (0.25 steps)
set.seed(0)
pts <- long %>% group_by(pair, metric, cls, value) %>%
  mutate(k = row_number() - 1L, n = n(),
         row_i = k %/% 6L, m = pmin(6L, n - 6L * row_i), j = k %% 6L,
         xp = ifelse(metric == "agree_frac", x + (j - (m - 1) / 2) * 0.075, x + runif(n(), -0.18, 0.18)),
         yp = ifelse(metric == "agree_frac", value + 0.045 * row_i, value)) %>% ungroup()
bars <- long %>% group_by(pair, metric, cls, x) %>%
  summarise(stat = if (first(metric) == "confidence") median(value) else mean(value), .groups = "drop")
ns <- D %>% count(pair, cls) %>% mutate(pair = factor(pair, PAIRS))
xlab_df <- ns %>% mutate(lab = sprintf("%s\n(n=%d)", c("Non-\nneuronal", "GABA", "Glut")[as.integer(cls)], n))
brk <- tests %>% filter(pair %in% PAIRS, test %in% c("Non-neuronal > GABAergic", "GABAergic > Glutamatergic")) %>%
  mutate(pair = factor(pair, PAIRS), metric = factor(metric, names(METRICS)),
         k = ifelse(test == "Non-neuronal > GABAergic", 1, 2), lab = pfmt(p), y = 1.11)

## facet_grid only ever draws the y-axis on its left column and the x-axis on
## its bottom row, even with scales = "free" (verified: ggplot2 does not repeat
## axes to interior/opposite panels in a grid facet -- that needs ggh4x or
## lemon, neither installed here). Built as four independent, unfaceted
## ggplots instead and assembled with patchwork, so every one draws its own
## full x/y axis, per direct instruction.
panel_plot <- function(met, pr, show_title, show_ylab) {
  b <- bars %>% filter(metric == met, pair == pr)
  pt <- pts %>% filter(metric == met, pair == pr)
  bk <- brk %>% filter(metric == met, pair == pr)
  g <- ggplot()
  if (met == "confidence")
    g <- g + geom_hline(aes(yintercept = 0.5), colour = FIG_INK_MUTED, linewidth = 0.3, linetype = "22")
  g <- g +
    geom_col(data = b, aes(x, stat, fill = cls), width = 0.62, alpha = 0.28) +
    geom_segment(data = b, aes(x = x - 0.31, xend = x + 0.31, y = stat, yend = stat, colour = cls), linewidth = 1.1) +
    geom_point(data = pt, aes(xp, yp, fill = cls), shape = 21, size = 1.1, colour = "white", stroke = 0.25) +
    geom_segment(data = bk, aes(x = k + 0.05, xend = k + 0.05, y = y - 0.02, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
    geom_segment(data = bk, aes(x = k + 0.95, xend = k + 0.95, y = y - 0.02, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
    geom_segment(data = bk, aes(x = k + 0.05, xend = k + 0.95, y = y, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
    ## Rotated: at this panel width two side-by-side horizontal p-value strings
    ## no longer fit within their own bracket's span and visually merge --
    ## vertical text only needs the bracket-span WIDTH for its stroke width,
    ## not for its string length, so it fits regardless.
    geom_text(data = bk, aes(k + 0.5, y + 0.03, label = lab), size = fig_pt(5.5), colour = FIG_INK_PRIMARY,
             angle = 90, hjust = 0, vjust = 0.5) +
    scale_fill_manual(values = CLASS_COL, guide = "none") +
    scale_colour_manual(values = CLASS_COL, guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = c("Non-\nneuronal", "GABA", "Glut"), expand = expansion(add = 0.55)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 2.0), expand = expansion(0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = if (show_ylab) METRICS[[met]] else NULL,
        title = if (show_title) pr else NULL) +
    theme_fig(base_size = 6.5) +
    theme(panel.grid.major.x = element_blank(), axis.line.y = element_line(colour = FIG_AXIS, linewidth = 0.3),
          plot.title = element_text(size = FIG_PT_AXIS_TITLE, hjust = 0.5))
  g
}

p <- (panel_plot("confidence", PAIRS[["gg_adult_hybrid"]], TRUE, TRUE) +
      panel_plot("confidence", PAIRS[["yao_adult_hybrid"]], TRUE, FALSE) +
      panel_plot("agree_frac", PAIRS[["gg_adult_hybrid"]], FALSE, TRUE) +
      panel_plot("agree_frac", PAIRS[["yao_adult_hybrid"]], FALSE, FALSE)) +
  plot_layout(ncol = 2) +
  plot_annotation(
    ## Short lines throughout -- at this panel width a sentence-per-line
    ## caption overran the right edge.
    caption = paste0("Points: finch clusters\n(non-neuronal n=", ns$n[1], ", GABA n=", ns$n[2], ", Glut n=", ns$n[3], ").\n",
                     "Bars: class median (top)\nor mean (bottom).\n",
                     "Dashed: high-confidence\ntier (0.5).\n",
                     "Brackets: one-sided\nMann-Whitney U,\nleft > right, p shown."),
    theme = theme(plot.caption = element_text(colour = FIG_INK_MUTED, size = 6.5 - 2.4, hjust = 0,
                                              family = FIG_FONT)))

fig_save(p, file.path(OUT, "class_integration_summary"), width = 2.66, height = 4.9)
message("wrote ", OUT)
