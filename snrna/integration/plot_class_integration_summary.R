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
pstar <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
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
         k = ifelse(test == "Non-neuronal > GABAergic", 1, 2), lab = pstar(p), y = 1.11)

facet_labs <- function(pair, met) sprintf("%s|%s", pair, met)
p <- ggplot() +
  geom_hline(data = data.frame(metric = factor("confidence", names(METRICS)), y = 0.5),
             aes(yintercept = y), colour = FIG_INK_MUTED, linewidth = 0.3, linetype = "22") +
  geom_col(data = bars, aes(x, stat, fill = cls), width = 0.62, alpha = 0.28) +
  geom_segment(data = bars, aes(x = x - 0.31, xend = x + 0.31, y = stat, yend = stat, colour = cls), linewidth = 1.1) +
  geom_point(data = pts, aes(xp, yp, fill = cls), shape = 21, size = 1.1, colour = "white", stroke = 0.25) +
  geom_segment(data = brk, aes(x = k + 0.05, xend = k + 0.05, y = y - 0.02, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
  geom_segment(data = brk, aes(x = k + 0.95, xend = k + 0.95, y = y - 0.02, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
  geom_segment(data = brk, aes(x = k + 0.05, xend = k + 0.95, y = y, yend = y), linewidth = 0.3, colour = FIG_INK_PRIMARY) +
  geom_text(data = brk, aes(k + 0.5, y + 0.02, label = lab), size = fig_pt(5.5), colour = FIG_INK_PRIMARY, vjust = 0) +
  facet_grid(metric ~ pair, switch = "y", labeller = labeller(metric = METRICS)) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("Non-\nneuronal", "GABA", "Glut"), expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.18), expand = expansion(0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL,
       caption = paste0("Points: finch clusters (", paste(sprintf("%s n=%d", c("non-neuronal", "GABA", "Glut"), ns$n[1:3]), collapse = ", "),
                        "). Bars: class median (top) or mean (bottom).\nDashed: high-confidence tier (0.5). ",
                        "Brackets: one-sided Mann-Whitney U, left > right;\n* p<0.05, ** p<0.01, *** p<0.001, ns p>0.05.")) +
  theme_fig(base_size = 6.5) +
  theme(panel.grid.major.x = element_blank(), axis.line.y = element_line(colour = FIG_AXIS, linewidth = 0.3),
        strip.placement = "outside", strip.text = element_text(colour = FIG_INK_SECONDARY, size = FIG_PT_AXIS_TITLE),
        strip.text.y.left = element_text(angle = 90), panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(4, "mm"))

fig_save(p, file.path(OUT, "class_integration_summary"), width = 3.7, height = 3.5)
message("wrote ", OUT)
