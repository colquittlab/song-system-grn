## ---------------------------------------------------------------------------
## Figure: Glut-DACH2-HVCra-Int is transcriptionally intermediate between
## Glut-DACH2-HVCra and the surrounding nidopallial Glut-DACH2-1.
##
## Three panels, three independent levels of evidence for the same claim:
##   A  per cell   -- every Glut cluster's position on the DACH2-1 -> HVCra axis,
##                    so HVCra-Int's position is read against the whole dataset
##                    rather than against its two anchors alone. Note what this
##                    panel does NOT claim: other song-nucleus types (HVCx,
##                    LMANco, CACNA1H-RA) also sit off both anchors. Mid-axis
##                    position alone does not make something an intermediate --
##                    those are discrete types, with 421 private markers in
##                    HVCx's case against HVCra-Int's 9 (see the QC section in
##                    snrna/naming/hybrid_division_naming.qmd). What is specific
##                    to HVCra-Int is holding a mid-axis position while having
##                    almost nothing of its own.
##   B  per gene set -- the genes that define the axis, shown graded across the
##                    three clusters.
##   C  per gene   -- each axis gene's shift in HVCra-Int against its shift in
##                    HVCra. Slope is the headline number: 0 = indistinguishable
##                    from DACH2-1, 1 = indistinguishable from HVCra.
##
## Reads snrna/naming/hybrid_division_naming.qmd's output read-only and writes
## only its own figures, by the same convention used throughout this project.
##
## Colours are the validated three-slot categorical palette (blue/orange/aqua),
## which clears the all-pairs CVD and normal-vision floors in both modes -- panel
## C is a scatter, so the all-pairs pairlist is the one that applies. Aqua sits
## below 3:1 on the light surface, so every cluster is also directly labelled
## rather than identified by colour alone.
## ---------------------------------------------------------------------------

suppressMessages({
  library(Seurat)
  library(qs2)
  library(tidyverse)
  library(patchwork)
})

obj_fname <- here::here("snrna/naming/hybrid_division_naming", "obj_hybrid_labels.qs2")
out_dir <- here::here("snrna/deg", "hvcra_int_intermediacy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

LOW <- "Glut-DACH2-1"           # the surrounding nidopallium
MID <- "Glut-DACH2-HVCra-Int"   # the cluster under test
HIGH <- "Glut-DACH2-HVCra"      # the song-nucleus projection neurons
N_AXIS_GENES <- 100             # per direction, for the axis score
N_HEATMAP_GENES <- 20           # per direction, for panel B

## Style --------------------------------------------------------------------
## Type, palette and theme come from the shared project standard rather than
## being re-declared here, so a change to the standard reaches every figure.
source(here::here("config/figure_theme.R"))

## Local aliases, so the plotting code below reads in its own terms.
pal <- FIG_PAL[seq_len(FIG_PAL_ALLPAIRS_MAX)]        # slots 1-3, fixed order
names(pal) <- c(LOW, MID, HIGH)
ink_primary <- FIG_INK_PRIMARY; ink_secondary <- FIG_INK_SECONDARY
ink_muted <- FIG_INK_MUTED; grid_col <- FIG_GRID; axis_col <- FIG_AXIS
surface <- FIG_SURFACE; neutral_fill <- FIG_NEUTRAL_FILL
div_low <- FIG_DIV_LOW; div_mid <- FIG_DIV_MID; div_high <- FIG_DIV_HIGH
FONT <- FIG_FONT
theme_viz <- theme_fig

## Data ---------------------------------------------------------------------
obj <- qs_read(obj_fname, nthreads = 8)
DefaultAssay(obj) <- "SCT"
obj <- obj[, !is.na(obj$celltype_hybrid)]
Idents(obj) <- obj$celltype_hybrid

## The axis genes: differential expression between the two anchors only. Nothing
## about HVCra-Int enters the definition of the axis, so its position on that
## axis is a prediction the data can refuse, not a fit.
de_anchor <- FindMarkers(obj, ident.1 = HIGH, ident.2 = LOW,
                         min.pct = 0.2, logfc.threshold = 0.5) %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.01)
high_genes <- de_anchor %>% arrange(desc(avg_log2FC)) %>% slice_head(n = N_AXIS_GENES) %>% pull(gene)
low_genes <- de_anchor %>% arrange(avg_log2FC) %>% slice_head(n = N_AXIS_GENES) %>% pull(gene)

scored <- AddModuleScore(obj, features = list(high_genes, low_genes),
                         name = "axis", assay = "SCT", seed = 1)
anchor_med <- tapply(scored$axis1 - scored$axis2, obj$celltype_hybrid, median)
rescale_axis <- function(x) (x - anchor_med[[LOW]]) / (anchor_med[[HIGH]] - anchor_med[[LOW]])

## `assignment` is souporcell's per-library index: the same numeral means a
## different bird in a different library, so a bird is (sample_name, assignment)
## and bird identity is only comparable WITHIN a library. Carried through here so
## the bird-level test below can key on it.
cells <- tibble(celltype = as.character(obj$celltype_hybrid),
                score = rescale_axis(scored$axis1 - scored$axis2),
                library = as.character(obj$sample_name),
                bird = paste0(obj$sample_name, "|", obj$assignment)) %>%
  filter(str_starts(celltype, "Glut"))

## Panel A ------------------------------------------------------------------
cell_order <- cells %>% group_by(celltype) %>% summarise(m = median(score)) %>%
  arrange(m) %>% pull(celltype)
cells <- cells %>%
  mutate(celltype = factor(celltype, levels = cell_order),
         role = if_else(celltype %in% names(pal), as.character(celltype), "other"))

panel_a <- ggplot(cells, aes(score, celltype, fill = role)) +
  geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.9) +
  stat_summary(fun = median, geom = "point", size = 1.4,
               colour = surface, stroke = 0, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", size = 0.9,
               colour = ink_primary, show.legend = FALSE) +
  scale_fill_manual(values = c(pal, other = neutral_fill), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("0\nDACH2-1", "0.5", "1\nHVCra")) +
  geom_vline(xintercept = c(0, 1), colour = axis_col, linewidth = 0.3, linetype = "22") +
  labs(title = "Every glutamatergic cluster on the DACH2-1 → HVCra axis",
       subtitle = paste0("Per-cell score on 200 genes defined by the two anchors alone,\n",
                         "rescaled so the anchor medians are 0 and 1. Other song-nucleus\n",
                         "types sit off the anchors too — see panel C for what is specific here."),
       x = "Position on axis", y = NULL, tag = "A") +
  theme_viz() +
  theme(panel.grid.major.y = element_blank())

## Direct labels for the three named clusters (relief rule: identity is never
## colour-alone).
panel_a <- panel_a +
  theme(axis.text.y = element_text(
    colour = if_else(cell_order %in% names(pal), ink_primary, ink_muted),
    face = if_else(cell_order %in% names(pal), "bold", "plain")))

## Panel B ------------------------------------------------------------------
hm_genes <- c(de_anchor %>% arrange(avg_log2FC) %>% slice_head(n = N_HEATMAP_GENES) %>% pull(gene),
              de_anchor %>% arrange(desc(avg_log2FC)) %>% slice_head(n = N_HEATMAP_GENES) %>% pull(gene))
expr <- LayerData(obj, assay = "SCT", layer = "data")[hm_genes, , drop = FALSE]
ct <- as.character(obj$celltype_hybrid)
pb <- sapply(c(LOW, MID, HIGH), function(k) Matrix::rowMeans(expr[, ct == k, drop = FALSE]))
pbz <- t(scale(t(pb)))

hm <- as_tibble(pbz, rownames = "gene") %>%
  pivot_longer(-gene, names_to = "celltype", values_to = "z") %>%
  mutate(celltype = factor(celltype, levels = c(LOW, MID, HIGH)),
         gene = factor(gene, levels = rev(hm_genes)))

panel_b <- ggplot(hm, aes(celltype, gene, fill = z)) +
  geom_tile(colour = surface, linewidth = 0.4) +
  scale_fill_gradient2(low = div_low, mid = div_mid, high = div_high, midpoint = 0,
                       name = "z-score", breaks = c(-1, 0, 1)) +
  scale_x_discrete(labels = function(x) str_remove(x, "^Glut-")) +
  labs(title = "The axis genes themselves",
       subtitle = paste0("Top ", N_HEATMAP_GENES, " per direction; mean expression, z-scored across the three.\n",
                         "The z-score does not force a middle column — HVCra-Int could sit at either pole."),
       x = NULL, y = NULL, tag = "B") +
  theme_viz() +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 5.2),
        axis.text.x = element_text(colour = ink_primary, face = "bold"),
        legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.5, "cm"))

## Panel C ------------------------------------------------------------------
## Per gene: how far has HVCra-Int moved, against how far HVCra moved, both
## measured from DACH2-1. A precursor or a distinct type has no reason to lie on
## a straight line through the origin; an intermediate does, and its slope is the
## fraction of the way it sits.
de_mid <- FindMarkers(obj, ident.1 = MID, ident.2 = LOW,
                      min.pct = 0.2, logfc.threshold = 0) %>%
  rownames_to_column("gene")
shift <- de_anchor %>%
  select(gene, hvcra = avg_log2FC) %>%
  inner_join(de_mid %>% select(gene, int = avg_log2FC), by = "gene")
fit <- coef(lm(int ~ 0 + hvcra, data = shift))[["hvcra"]]

panel_c <- ggplot(shift, aes(hvcra, int)) +
  geom_hline(yintercept = 0, colour = axis_col, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = axis_col,
              linewidth = 0.4, linetype = "22") +
  geom_point(colour = ink_secondary, alpha = 0.35, size = 0.9, stroke = 0) +
  geom_abline(slope = fit, intercept = 0, colour = pal[[MID]], linewidth = 0.9) +
  annotate("text", x = max(shift$hvcra), y = max(shift$hvcra) * 1.02,
           label = "identical to HVCra", hjust = 1, vjust = 0,
           size = 2.5, colour = ink_muted) +
  annotate("text", x = max(shift$hvcra), y = 0.15,
           label = "identical to DACH2-1", hjust = 1, vjust = 0,
           size = 2.5, colour = ink_muted) +
  annotate("text", x = max(shift$hvcra) * 0.98, y = max(shift$hvcra) * fit,
           label = sprintf("HVCra-Int: %.0f%% of the way", 100 * fit),
           hjust = 1, vjust = -0.6, size = 2.9, colour = pal[[MID]], fontface = "bold") +
  labs(title = "Each axis gene, shifted part-way",
       subtitle = "log2 fold-change from Glut-DACH2-1, per gene",
       x = "HVCra vs DACH2-1", y = "HVCra-Int vs DACH2-1", tag = "C") +
  theme_viz()

## ---------------------------------------------------------------------------
## Compact standalone version of panel C, for use on its own at figure scale.
## Same data and same fit as the panel above -- only the text budget and the
## mark sizes are re-cut for a 2.8 in measure, so the two cannot disagree.
## ---------------------------------------------------------------------------
sx_max <- max(shift$hvcra)
sx_diag <- min(sx_max, max(shift$int)) * 0.92
scatter_compact <- ggplot(shift, aes(hvcra, int)) +
  geom_hline(yintercept = 0, colour = axis_col, linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, colour = axis_col,
              linewidth = 0.3, linetype = "22") +
  geom_point(colour = ink_secondary, alpha = 0.3, size = 0.45, stroke = 0) +
  geom_abline(slope = fit, intercept = 0, colour = pal[[MID]], linewidth = 0.7) +
  ## Anchor the reference labels inside the panel. The y = x line leaves the top
  ## of the panel well before the x data end, so labelling it at max(hvcra)
  ## places the text outside the drawing area and it is silently cropped.
  annotate("text", x = sx_diag, y = sx_diag, label = "= HVCra",
           hjust = 1, vjust = -0.5, size = 5 / .pt, colour = ink_muted) +
  annotate("text", x = sx_max, y = 0, label = "= DACH2-1",
           hjust = 1, vjust = -0.7, size = 5 / .pt, colour = ink_muted) +
  annotate("text", x = sx_max, y = sx_max * fit,
           label = sprintf("%.0f%% of the way", 100 * fit),
           hjust = 1, vjust = -0.8, size = 6 / .pt,
           colour = pal[[MID]], fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Each axis gene shifts part-way",
       subtitle = "log2 fold-change from DACH2-1, per gene",
       x = "HVCra vs DACH2-1", y = "HVCra-Int vs DACH2-1") +
  theme_viz(base_size = 7) +
  theme(plot.margin = margin(2, 4, 2, 2))

ggsave(file.path(out_dir, "hvcra_int_gene_shift_scatter.pdf"), scatter_compact,
       width = 2.8, height = 2.2, device = cairo_pdf, family = FONT)
ggsave(file.path(out_dir, "hvcra_int_gene_shift_scatter.png"), scatter_compact,
       width = 2.8, height = 2.2, dpi = 600, bg = surface)

## Compose ------------------------------------------------------------------
fig <- panel_a | (panel_b / panel_c + plot_layout(heights = c(1.15, 1)))
fig <- fig + plot_layout(widths = c(1.1, 1)) &
  theme(plot.background = element_rect(fill = surface, colour = NA))

ggsave(file.path(out_dir, "hvcra_int_intermediacy.pdf"), fig,
       width = 11, height = 7.5, device = cairo_pdf, family = FONT)
ggsave(file.path(out_dir, "hvcra_int_intermediacy.png"), fig,
       width = 11, height = 7.5, dpi = 300, bg = surface)

## The numbers behind the figure, so the claim is checkable without re-running it.
cells %>%
  group_by(celltype) %>%
  summarise(n_cells = n(), median_axis_position = round(median(score), 3)) %>%
  arrange(median_axis_position) %>%
  write_csv(file.path(out_dir, "axis_position_by_celltype.csv"))
write_csv(shift, file.path(out_dir, "per_gene_shift.csv"))

## ---------------------------------------------------------------------------
## Focused three-way violin: the hypothesis on its own terms
##
## Panel A above deliberately shows all 23 glutamatergic clusters, because the
## honest version of the claim has to survive the whole dataset being in frame.
## This is the same axis restricted to the three clusters the hypothesis is
## actually about, for when that is the question being asked.
##
## The medians carry a bootstrap 95% CI so this reads as a test rather than a
## picture: the question "is HVCra-Int strictly between the anchors" is answered
## by whether its interval excludes 0 and 1, not by whether the violin looks
## like it sits in the middle. Note the axis is defined by DE between the two
## anchors alone -- HVCra-Int contributes nothing to it, so its position is a
## prediction that could have come back at either end.
## ---------------------------------------------------------------------------

boot_median_ci <- function(x, n_boot = 2000, seed = 1) {
  set.seed(seed)
  b <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
  c(lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)))
}

three <- cells %>%
  filter(celltype %in% c(LOW, MID, HIGH)) %>%
  mutate(celltype = factor(as.character(celltype), levels = c(LOW, MID, HIGH)))

three_stats <- three %>%
  group_by(celltype) %>%
  summarise(n_cells = n(),
            median = median(score),
            lo = boot_median_ci(score)[["lo"]],
            hi = boot_median_ci(score)[["hi"]],
            q25 = quantile(score, 0.25),
            q75 = quantile(score, 0.75),
            iqr_width = q75 - q25,
            .groups = "drop")

## ---------------------------------------------------------------------------
## Bird-level test: the same comparison with birds, not cells, as the replicate.
##
## The cell-level Mann-Whitney above is inflated -- cells within a bird are not
## independent draws -- so this repeats it on per-bird medians.
##
## What the design actually allows, checked rather than assumed:
## Glut-DACH2-HVCra exists ONLY in the hvc_run1 library, and bird identity is not
## comparable across libraries, so the only design in which all three clusters
## share the same birds is within hvc_run1. That is n = 3 birds, paired. The
## DACH2-1 arm is therefore its hvc-library cells (34/121/170 per bird), i.e. the
## surrounding nidopallium captured in the HVC dissection, not the whole cluster.
##
## n = 3 is the honest limit here, and it is a hard one: a paired rank test on
## three pairs cannot return a two-sided p below 0.25 no matter how large the
## effect, because there are only 2^3 = 8 sign assignments. Both tests are
## reported below so that floor is visible rather than hidden behind whichever
## one happens to look better. The paired t-test on three bird medians can clear
## 0.05, but it buys that by assuming normality of three differences, which three
## points cannot evidence. Treat the per-bird points on the figure -- 3/3 birds
## ordered the same way, tightly grouped -- as the real bird-level evidence, and
## these p-values as a formality.
BIRD_LIB <- "hvc_run1"
BIRD_MIN_CELLS <- 20

bird_cells <- three %>% filter(library == BIRD_LIB)
bird_med <- bird_cells %>%
  group_by(celltype, bird) %>%
  summarise(n_cells = n(), median = median(score), .groups = "drop") %>%
  filter(n_cells >= BIRD_MIN_CELLS)

## Keep only birds seen in all three clusters, so the test is genuinely paired.
paired_birds <- bird_med %>% count(bird) %>% filter(n == 3) %>% pull(bird)
bird_med <- bird_med %>% filter(bird %in% paired_birds) %>%
  mutate(celltype = factor(as.character(celltype), levels = c(LOW, MID, HIGH)))
message("bird-level replicates (", BIRD_LIB, ", >= ", BIRD_MIN_CELLS,
        " cells in all three clusters): n = ", length(paired_birds))

bird_wide <- bird_med %>%
  select(celltype, bird, median) %>%
  pivot_wider(names_from = celltype, values_from = median)

bird_test <- map_dfr(list(c(LOW, MID), c(MID, HIGH), c(LOW, HIGH)), function(p) {
  x <- bird_wide[[p[1]]]; y <- bird_wide[[p[2]]]
  tt <- t.test(x, y, paired = TRUE)
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = TRUE))
  tibble(group1 = p[1], group2 = p[2], n_birds = length(x),
         mean_diff = mean(y - x),
         t_p = tt$p.value, wilcox_p = wt$p.value,
         same_direction_in_all_birds = all((y - x) > 0) || all((y - x) < 0))
}) %>%
  mutate(t_p_adj = p.adjust(t_p, method = "BH"),
         wilcox_p_adj = p.adjust(wilcox_p, method = "BH"))
write_csv(bird_test, file.path(out_dir, "three_way_bird_level_test.csv"))
write_csv(bird_med, file.path(out_dir, "three_way_bird_medians.csv"))

message("bird-level paired tests:")
print(as.data.frame(bird_test %>% mutate(across(where(is.numeric), ~signif(.x, 3)))))

## Only HVCra-Int carries a numeric label. The two anchors are 0 and 1 *by
## construction* of the rescaling, so printing a CI on them would invite reading
## them as independent measurements of something -- they are the definition of
## the scale, and their bootstrap interval is just sampling noise around a fixed
## point. HVCra-Int's interval is the only one that is a result.
three_labels <- three_stats %>%
  mutate(label = if_else(celltype == MID,
                         sprintf("%.2f [%.2f, %.2f]", median, lo, hi),
                         NA_character_))

## Mann-Whitney U (Wilcoxon rank-sum) between the three groups, all three
## pairs, Benjamini-Hochberg corrected across them.
##
## Read the effect size, not the p-value. Each group is hundreds to thousands of
## *cells* drawn from three birds in a single library per dissection, so cells
## are not independent samples: n is inflated relative to the number of animals
## and any real difference is guaranteed to come back astronomically significant.
## The p-value is a statement about cells, not about birds. `prob_superiority`
## (U / n1n2 -- the chance a random cell from group 2 scores above a random cell
## from group 1) is the number that is comparable across pairs and unaffected by
## that inflation.
##
## Note also what this test can and cannot establish. It shows HVCra-Int differs
## from BOTH anchors, which rules out "it is really just one of them". It does
## NOT establish intermediacy -- a cluster off one end of the axis would give the
## same two rejections. Intermediacy is carried by the position estimate (0.42,
## CI [0.40, 0.45]), not by these p-values.
mwu_pairs <- list(c(LOW, MID), c(MID, HIGH), c(LOW, HIGH))
mwu <- map_dfr(mwu_pairs, function(p) {
  x <- three$score[three$celltype == p[1]]
  y <- three$score[three$celltype == p[2]]
  w <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided", exact = FALSE))
  tibble(group1 = p[1], group2 = p[2], n1 = length(x), n2 = length(y),
         U = unname(w$statistic), p_value = w$p.value,
         prob_superiority = 1 - unname(w$statistic) / (length(x) * length(y)),
         median_diff = median(y) - median(x))
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         rank_biserial = 2 * prob_superiority - 1)
write_csv(mwu, file.path(out_dir, "three_way_mannwhitney.csv"))

fmt_p <- function(p) if_else(p < 2.2e-16, "< 2.2e-16", sprintf("= %.2g", p))
## %.2f rounded 5.7e-04 to "0.00"; keep two significant figures instead.
fmt_p_sig <- function(p) sprintf("%.2g", p)
mwu_caption <- sprintf(
  "Anchors 0 and 1 by definition. Bar = bootstrap 95%% CI; box = IQR.\nPoints above = per-bird medians (%2$s, n = %3$d birds, paired).\nPaired t on bird medians p = %4$s / %5$s; cell-wise MWU p %1$s but inflated.\nA rank test on 3 pairs cannot go below p = 0.25 — 3/3 ordering is the evidence.",
  fmt_p(max(mwu$p_adj[mwu$group1 == MID | mwu$group2 == MID])),
  BIRD_LIB, length(paired_birds),
  fmt_p_sig(bird_test$t_p_adj[bird_test$group1 == LOW & bird_test$group2 == MID]),
  fmt_p_sig(bird_test$t_p_adj[bird_test$group1 == MID & bird_test$group2 == HIGH]))

BIRD_OFFSET <- 0.38
three_violin <- ggplot(three, aes(score, celltype, fill = celltype)) +
  geom_vline(xintercept = c(0, 1), colour = axis_col, linewidth = 0.35, linetype = "22") +
  geom_violin(scale = "width", width = 0.7, colour = NA, alpha = 0.9) +
  geom_boxplot(width = 0.09, outlier.shape = NA, colour = ink_primary,
               fill = surface, linewidth = 0.3) +
  geom_linerange(data = three_stats,
                 aes(y = celltype, xmin = lo, xmax = hi),
                 colour = ink_primary, linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(data = three_stats, aes(median, celltype),
             size = 2.6, colour = surface, inherit.aes = FALSE) +
  geom_point(data = three_stats, aes(median, celltype),
             size = 1.7, colour = ink_primary, inherit.aes = FALSE) +
  geom_text(data = three_labels %>% filter(!is.na(label)),
            aes(median, celltype, label = label),
            vjust = 2.6, size = 6 / .pt, colour = ink_primary,
            fontface = "bold", inherit.aes = FALSE) +
  ## Per-bird medians, offset above each violin. With n = 3 the spread of these
  ## points is the bird-level evidence -- see three_way_bird_level_test.csv for
  ## the paired tests, and note a rank test on 3 pairs floors at p = 0.25.
  geom_point(data = bird_med,
             aes(x = median, y = as.numeric(celltype) + BIRD_OFFSET),
             size = 1.5, colour = surface, inherit.aes = FALSE) +
  geom_point(data = bird_med,
             aes(x = median, y = as.numeric(celltype) + BIRD_OFFSET),
             size = 0.85, colour = ink_primary, inherit.aes = FALSE) +
  scale_fill_manual(values = pal, guide = "none") +
  ## Explicit short labels rather than a prefix strip -- stripping "Glut-DACH2-"
  ## turns Glut-DACH2-1 into a bare "1", which reads as a number, not a cluster.
  scale_y_discrete(labels = function(x) {
    short <- c("DACH2-1", "HVCra-Int", "HVCra")
    names(short) <- c(LOW, MID, HIGH)
    n <- three_stats$n_cells[match(x, as.character(three_stats$celltype))]
    paste0(short[x], "\n(", n, ")")
  }) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = 0.02)) +
  coord_cartesian(ylim = c(0.55, 3.62), clip = "off") +
  labs(title = "HVCra-Int is intermediate",
       subtitle = paste0("Per-cell score, ", 2 * N_AXIS_GENES, " anchor-differential genes"),
       x = "Position on the DACH2-1 → HVCra axis", y = NULL,
       caption = mwu_caption) +
  theme_viz(base_size = 7) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary, face = "bold", size = 6),
        plot.caption = element_text(colour = ink_muted, size = 4.6, hjust = 0,
                                    margin = margin(t = 1.5)),
        plot.margin = margin(2, 4, 2, 2))

## 40% off both dimensions of the previous 5.0 x 3.4 in.
ggsave(file.path(out_dir, "hvcra_int_three_way_violin.pdf"), three_violin,
       width = 3, height = 2.5, device = cairo_pdf, family = FONT)
ggsave(file.path(out_dir, "hvcra_int_three_way_violin.png"), three_violin,
       width = 3, height = 2.5, dpi = 600, bg = surface)
write_csv(three_stats %>% mutate(across(median:iqr_width, ~round(.x, 4))),
          file.path(out_dir, "three_way_axis_medians.csv"))

message("Mann-Whitney U:")
print(as.data.frame(mwu %>% mutate(across(c(prob_superiority, median_diff, rank_biserial), ~round(.x, 3)))))
message("three-way medians [95% CI]:")
print(as.data.frame(three_stats %>% mutate(across(median:iqr_width, ~round(.x, 3)))))

message("slope (fraction of the way HVCra-Int sits): ", round(fit, 3))
message("Wrote ", out_dir)
