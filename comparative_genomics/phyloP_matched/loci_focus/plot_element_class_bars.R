## Number of MAFB-domain accessible elements (lobes) in each phyloP threshold class, by DAR class.
##
## Reads results/mafb_domain/element_threshold_classes.tsv (element_threshold_classes.py): each lobe's
## threshold class from its accessibility-weighted whole-tree and oscine mean phyloP over its full
## extent (deep+oscine, deep, oscine, weak), and its DAR class (RA, G4, shared = "both"). Panels: all
## lobes, and lobes with every 50-bp window aligned in >= 0.8 / >= 0.9 of oscine and non-oscine
## genomes. Positions are DAR classes (RA, G4, shared), bars within each coloured by phyloP class (as
## in the threshold heatmaps); counts are printed above each bar; y axis shared.
##
##   Rscript plot_element_class_bars.R [oscany]
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
d <- fread(file.path(res, "element_threshold_classes.tsv"))
CLS <- c(RA = "Glut-RA", G4 = "GABA-4", both = "shared")
PCOL <- c("deep+oscine" = "#7a3ca0", deep = "#c8161a", oscine = "#1baf7a", weak = "#5f5f5f")
d[, dar := factor(CLS[cls], levels = CLS)]
d[, phy := factor(thr_class, levels = c("deep+oscine", "deep", "oscine", "weak"))]
if (oscany) {
  d[, phy := factor(fcase(thr_class %in% c("deep+oscine", "oscine"), "any oscine", thr_class == "deep", "deep only", default = "weak"),
                    levels = c("any oscine", "deep only", "weak"))]
  PCOL <- c("any oscine" = "#1baf7a", "deep only" = "#c8161a", weak = "#5f5f5f")
}
sets <- list("all lobes" = d, "depth ≥ 0.8" = d[osc_depth_min50 >= 0.8 & nonosc_depth_min50 >= 0.8],
             "depth ≥ 0.9" = d[osc_depth_min50 >= 0.9 & nonosc_depth_min50 >= 0.9])
n <- rbindlist(lapply(names(sets), function(s) sets[[s]][, .N, by = .(phy, dar)][, set := s]))
n <- n[CJ(phy = levels(d$phy), dar = levels(d$dar), set = names(sets)), on = .(phy, dar, set)][is.na(N), N := 0L]
n[, set := factor(set, levels = names(sets), labels = sprintf("%s (n = %d)", names(sets), sapply(sets, nrow)))]
n[, phy := factor(phy, levels = levels(d$phy))]; n[, dar := factor(dar, levels = CLS)]
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
p <- ggplot(n, aes(dar, N, fill = phy)) +
  geom_col(position = pd, width = 0.72, colour = FIG_SURFACE, linewidth = 0.25) +
  geom_text(aes(label = N), position = pd, vjust = -0.4, size = FIG_PT_AXIS_TEXT / .pt, colour = FIG_INK_SECONDARY) +
  scale_fill_manual(values = PCOL, name = "phyloP class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~set, nrow = 1) +
  labs(x = "DAR class", y = "elements (lobes)") +
  theme_fig() +
  theme(strip.text = element_text(size = FIG_PT_AXIS_TITLE, colour = FIG_INK_SECONDARY),
        panel.grid.major.x = element_blank(), panel.spacing = unit(4, "mm"),
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        legend.key.size = unit(2.5, "mm"))
stem <- paste0("fig_element_class_bars", if (oscany) "_oscany" else "")
fig_save(p, file.path(res, stem), width = 5.6, height = 2.1)
message("wrote ", file.path(res, stem), ".pdf/.png")
