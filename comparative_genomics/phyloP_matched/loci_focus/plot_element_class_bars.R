## Number of MAFB-domain accessible elements (lobes) in each phyloP threshold class, by DAR class.
##
## Reads results/mafb_domain/element_threshold_classes.tsv (element_threshold_classes.py): each lobe's
## threshold class from its accessibility-weighted whole-tree and oscine mean phyloP over its full
## extent (deep+oscine, strong avian, moderate avian, oscine, weak), and its DAR class (RA, G4, shared = "both"). Panels: all
## lobes, and lobes with every 50-bp window aligned in >= 0.8 / >= 0.9 of oscine and non-oscine
## genomes. Positions are DAR classes (RA, G4, shared), bars within each coloured by phyloP class (as
## in the threshold heatmaps); counts are printed above each bar; y axis shared.
##
##   Rscript plot_element_class_bars.R [oscany] [pct | stack | horiz | mini]
## mini: the stacked bars for depth >= 0.8 only, in a 1.25 x 1.25 in panel (5 pt text), with the legend
## written to its own file (_stack_depth0.8_mini[_legend]).
## horiz draws the stacked bars horizontally, panels one above another; outputs carry _hstack.
## pct scales each bar to the % of its DAR class's lobes (per panel), keeping the count above the bar;
## outputs carry _pct. stack draws one 100% stacked bar per DAR class (counts inside segments of >= 6%,
## totals above); outputs carry _stack. Stacked plots carry brackets with the pairwise Holm-adjusted p.
## oscany combines deep+oscine and oscine into "any oscine" (any oscine-lineage constraint), next to
## "deep only" and "weak", and reports per DAR class the share with any oscine constraint (RA vs G4,
## Fisher's exact test); outputs carry _oscany. Without it, the class mix is tested RA vs G4 (and with
## shared) by Fisher's exact test (simulated p, 1e6) and each class RA vs G4 by Fisher's exact test.
## Writes results/mafb_domain/fig_element_class_bars[_oscany].{pdf,png}.
suppressMessages({library(ggplot2); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
oscany <- "oscany" %in% commandArgs(trailingOnly = TRUE)
pct <- "pct" %in% commandArgs(trailingOnly = TRUE)    # bar height = % of the DAR class's lobes; count above
stack <- "stack" %in% commandArgs(trailingOnly = TRUE) # one stacked 100% bar per DAR class
horiz <- "horiz" %in% commandArgs(trailingOnly = TRUE) # horizontal stacked bars, panels stacked vertically
if (horiz) stack <- TRUE
mini <- "mini" %in% commandArgs(trailingOnly = TRUE)   # depth >= 0.8 only, 1.25 x 1.25 in; legend saved separately
if (mini) { stack <- TRUE; horiz <- FALSE }
if (stack) pct <- TRUE
d <- fread(file.path(res, "element_threshold_classes.tsv"))
CLS <- c(RA = "Glut-RA", G4 = "GABA-4", both = "shared")
## phyloP classes kept clear of the DAR-class colours (Glut-RA red, GABA-4 blue): avian constraint in the heatmap's
## purples (dark = strong), oscine green, both ochre, weak light grey
PCOL <- c("deep+oscine" = "#d4a017", "strong avian" = "#3f007d", "moderate avian" = "#8c6bb1", oscine = "#1baf7a", weak = "#bdbdbd")
DARK_FILL <- c("strong avian", "moderate avian", "deep only")   # white counts on these, dark ink on the rest
d[, dar := factor(CLS[cls], levels = CLS)]
d[, phy := factor(thr_class, levels = names(PCOL))]
if (oscany) {
  d[, phy := factor(fcase(thr_class %in% c("deep+oscine", "oscine"), "any oscine", thr_class %in% c("strong avian", "moderate avian"), "deep only", default = "weak"),
                    levels = c("any oscine", "deep only", "weak"))]
  PCOL <- c("any oscine" = "#1baf7a", "deep only" = "#5e3c99", weak = "#bdbdbd")
}
sets <- list("all lobes" = d, "depth ≥ 0.8" = d[osc_depth_min50 >= 0.8 & nonosc_depth_min50 >= 0.8],
             "depth ≥ 0.9" = d[osc_depth_min50 >= 0.9 & nonosc_depth_min50 >= 0.9])
if (mini) sets <- sets["depth ≥ 0.8"]
n <- rbindlist(lapply(names(sets), function(s) sets[[s]][, .N, by = .(phy, dar)][, set := s]))
n <- n[CJ(phy = levels(d$phy), dar = levels(d$dar), set = names(sets)), on = .(phy, dar, set)][is.na(N), N := 0L]
n[, set := factor(set, levels = names(sets), labels = sprintf("%s (n = %d)", names(sets), sapply(sets, nrow)))]
n[, phy := factor(phy, levels = levels(d$phy))]; n[, dar := factor(dar, levels = CLS)]
n[, P := if (sum(N)) 100 * N / sum(N) else 0, by = .(set, dar)]
print(dcast(n, set + phy ~ dar, value.var = "N"))
if (!oscany) {                                       # exact tests of the class mix, and per class RA vs G4
  set.seed(1)
  for (s in names(sets)) {
    x <- sets[[s]]; t2 <- table(x[cls %in% c("RA", "G4"), .(phy, cls = factor(cls, c("RA", "G4")))])
    t3 <- table(x[, .(phy, cls = factor(cls, c("RA", "G4", "both")))])
    message(sprintf("%s (n = %d): class mix RA vs G4 Fisher p = %.3g; RA vs G4 vs shared p = %.3g (simulated, 1e6)", s, nrow(x),
                    fisher.test(t2, simulate.p.value = TRUE, B = 1e6)$p.value, fisher.test(t3, simulate.p.value = TRUE, B = 1e6)$p.value))
    for (k in levels(x$phy)) {
      a <- t2[k, "RA"]; na <- sum(t2[, "RA"]); b <- t2[k, "G4"]; nb <- sum(t2[, "G4"])
      message(sprintf("    %-12s RA %2d/%-2d (%3.0f%%)  G4 %2d/%-2d (%3.0f%%)  Fisher p = %.3g", k, a, na, 100 * a / na, b, nb, 100 * b / nb,
                      fisher.test(matrix(c(a, na - a, b, nb - b), 2))$p.value))
    }
  }
}
if (oscany) for (s in names(sets)) {                 # share with any oscine constraint, RA vs G4
  x <- sets[[s]]; f <- function(k) x[cls == k, c(sum(phy == "any oscine"), .N)]
  r <- f("RA"); g <- f("G4"); b <- f("both")
  ft <- fisher.test(matrix(c(r[1], r[2] - r[1], g[1], g[2] - g[1]), 2))
  message(sprintf("%-12s any oscine: RA %d/%d (%.0f%%), G4 %d/%d (%.0f%%), shared %d/%d (%.0f%%); RA vs G4 Fisher p = %.2g",
                  s, r[1], r[2], 100 * r[1] / r[2], g[1], g[2], 100 * g[1] / g[2], b[1], b[2], 100 * b[1] / b[2], ft$p.value))
}

pd <- position_dodge(width = 0.8)
tot <- n[, .(N = sum(N)), by = .(set, dar)]
## pairwise tests of the class distribution (RA vs G4, G4 vs shared, RA vs shared): Fisher's exact test (simulated
## p, 1e6), Holm-adjusted across the three pairs within each panel; drawn as brackets on the stacked plots
if (stack) {
  set.seed(1)
  PR <- list(CLS[c("RA", "G4")], CLS[c("G4", "both")], CLS[c("RA", "both")])
  xpos <- if (horiz) rev(unname(CLS)) else unname(CLS)
  br <- rbindlist(lapply(names(sets), function(s) {
    x <- sets[[s]]
    ps <- sapply(PR, function(pr) { t <- table(droplevels(x[dar %in% pr, .(phy, dar = factor(dar, pr))])); t <- t[rowSums(t) > 0, , drop = FALSE]
                                    fisher.test(t, simulate.p.value = TRUE, B = 1e6)$p.value })
    data.table(set = s, x1 = sapply(PR, function(pr) match(pr[1], xpos)), x2 = sapply(PR, function(pr) match(pr[2], xpos)),
               lev = 1:3, p = p.adjust(ps, "holm"))
  }))
  br[, set := factor(set, levels = names(sets), labels = levels(n$set))]
  br[, y := (if (mini) 116 else 113) + (lev - 1) * (if (mini) 10 else 8)]
  br[, lab := ifelse(p < 0.001, "p < 0.001", sprintf("p = %s", formatC(signif(p, 2), format = "fg", digits = 2)))]
  print(br[, .(set, pair = sapply(PR, paste, collapse = " vs ")[lev], p_holm = signif(p, 3))])
}
p <- if (stack) {
  ## counts inside segments (white; omitted below 6%), class totals above the bars
  ggplot(n, aes(dar, P, fill = phy)) +
    geom_col(width = 0.7, colour = FIG_SURFACE, linewidth = 0.4) +
    geom_text(aes(label = ifelse(P >= 6, N, ""), colour = phy %in% DARK_FILL, group = phy), position = position_stack(vjust = 0.5),
              size = FIG_PT_AXIS_TEXT / .pt, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#171B22"), guide = "none") +
    geom_text(data = tot, aes(dar, 100, label = sprintf("n = %d", N)), inherit.aes = FALSE,
              vjust = if (horiz) 0.5 else -0.5, hjust = if (horiz) -0.15 else 0.5,
              size = FIG_PT_AXIS_TEXT / .pt, colour = FIG_INK_SECONDARY)
} else {
  ggplot(n, aes(dar, if (pct) P else N, fill = phy)) +
    geom_col(position = pd, width = 0.72, colour = FIG_SURFACE, linewidth = 0.25) +
    geom_text(aes(label = N), position = pd, vjust = -0.4, size = FIG_PT_AXIS_TEXT / .pt, colour = FIG_INK_SECONDARY)
}
p <- p +
  scale_fill_manual(values = PCOL, name = "phyloP class") +
  (if (stack) scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 100, 25),
                                 limits = c(0, if (horiz) 150 else max(br$y) + 7))
   else scale_y_continuous(expand = expansion(mult = c(0, 0.12)))) +
  facet_wrap(~set, nrow = if (horiz) 3 else 1) +
  labs(x = "DAR class", y = if (pct) "% of the DAR class's elements (lobes)" else "elements (lobes)") +
  theme_fig() +
  theme(strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
        panel.grid.major.x = element_blank(), panel.spacing = unit(4, "mm"),
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        legend.key.size = unit(2.5, "mm"))
if (stack && !horiz) {                               # p-value brackets (ticks point back towards the bars)
  bs <- if (mini) 4 / .pt else FIG_PT_AXIS_TEXT / .pt
  seg <- rbind(br[, .(set, x = x1, xend = x2, y, yend = y)], br[, .(set, x = x1, xend = x1, y, yend = y - 2)],
               br[, .(set, x = x2, xend = x2, y, yend = y - 2)])
  p <- p + geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                        linewidth = 0.25, colour = FIG_INK_SECONDARY) +
    geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + 0.8, label = lab), inherit.aes = FALSE,
              vjust = 0, hjust = 0.5, size = bs, colour = FIG_INK_SECONDARY)
}
if (horiz) {                                         # horizontal: the three pairwise p as a list right of the bars
  pl <- br[, .(lab = paste(sprintf("%s: %s", c("RA vs G4", "G4 vs shared", "RA vs shared")[lev], sub("^p ", "p ", lab)), collapse = "\n")), by = set]
  p <- p + geom_text(data = pl, aes(x = 2, y = 113, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0.5,
                     lineheight = 0.95, size = FIG_PT_AXIS_TEXT / .pt, colour = FIG_INK_SECONDARY)
}
if (horiz) {                                         # RA at the top; the stack reads left to right in legend order
  p <- p + coord_flip() + scale_x_discrete(limits = unname(rev(CLS))) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25), panel.grid.major.y = element_blank(),
          legend.position = "bottom")
}
if (mini) {                                          # small panel: no facet strip, 5 pt text, legend in its own file
  TXT <- 5
  p <- p + facet_null() + labs(x = NULL, y = "% of elements", title = NULL) +
    theme(legend.position = "none", axis.text = element_text(size = TXT), axis.title.y = element_text(size = TXT),
          axis.text.x = element_text(size = TXT, colour = FIG_INK_SECONDARY), plot.margin = margin(2, 2, 1, 1, "pt"))
  p$layers[[2]]$aes_params$size <- 4 / .pt; p$layers[[3]]$aes_params$size <- 4 / .pt
  p <- p + scale_x_discrete(labels = c("Glut-RA" = "RA", "GABA-4" = "G4", shared = "shared"))
  stem <- paste0("fig_element_class_bars", if (oscany) "_oscany" else "", "_stack_depth0.8_mini")
  fig_save(p, file.path(res, stem), width = 1.25, height = 1.25)
  leg <- cowplot::get_legend(p + theme(legend.position = "right", legend.title = element_text(size = TXT),
                                       legend.text = element_text(size = TXT), legend.key.size = unit(2.2, "mm")))
  fig_save(cowplot::ggdraw(leg), file.path(res, paste0(stem, "_legend")), width = 0.9, height = 0.8)
  message("wrote ", file.path(res, stem), ".pdf/.png and _legend")
  quit(save = "no")
}
stem <- paste0("fig_element_class_bars", if (oscany) "_oscany" else "", if (horiz) "_hstack" else if (stack) "_stack" else if (pct) "_pct" else "")
fig_save(p, file.path(res, stem), width = if (horiz) 3.6 else 5.6, height = if (horiz) 3.4 else 2.1)
message("wrote ", file.path(res, stem), ".pdf/.png")
