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

## Palette ------------------------------------------------------------------
pal <- c("#2a78d6", "#eb6834", "#1baf7a")           # slots 1-3, fixed order
names(pal) <- c(LOW, MID, HIGH)
ink_primary <- "#0b0b0b"; ink_secondary <- "#52514e"; ink_muted <- "#898781"
grid_col <- "#e1e0d9"; axis_col <- "#c3c2b7"; surface <- "#fcfcfb"
neutral_fill <- "#c9c8c2"                            # context clusters
div_low <- "#2a78d6"; div_mid <- "#f0efec"; div_high <- "#e34948"

theme_viz <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = surface, colour = NA),
      panel.background = element_rect(fill = surface, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.25),
      axis.line = element_line(colour = axis_col, linewidth = 0.3),
      axis.ticks = element_line(colour = axis_col, linewidth = 0.3),
      axis.text = element_text(colour = ink_muted),
      axis.title = element_text(colour = ink_secondary),
      plot.title = element_text(colour = ink_primary, face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = ink_secondary, size = base_size - 0.5),
      plot.tag = element_text(colour = ink_primary, face = "bold", size = base_size + 3),
      legend.text = element_text(colour = ink_secondary),
      legend.title = element_text(colour = ink_secondary)
    )
}

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

cells <- tibble(celltype = as.character(obj$celltype_hybrid),
                score = rescale_axis(scored$axis1 - scored$axis2)) %>%
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

## Compose ------------------------------------------------------------------
fig <- panel_a | (panel_b / panel_c + plot_layout(heights = c(1.15, 1)))
fig <- fig + plot_layout(widths = c(1.1, 1)) &
  theme(plot.background = element_rect(fill = surface, colour = NA))

ggsave(file.path(out_dir, "hvcra_int_intermediacy.pdf"), fig,
       width = 11, height = 7.5, device = cairo_pdf)
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

## Only HVCra-Int carries a numeric label. The two anchors are 0 and 1 *by
## construction* of the rescaling, so printing a CI on them would invite reading
## them as independent measurements of something -- they are the definition of
## the scale, and their bootstrap interval is just sampling noise around a fixed
## point. HVCra-Int's interval is the only one that is a result.
three_labels <- three_stats %>%
  mutate(label = if_else(celltype == MID,
                         sprintf("median %.2f   95%% CI [%.2f, %.2f]", median, lo, hi),
                         NA_character_))

three_violin <- ggplot(three, aes(score, celltype, fill = celltype)) +
  geom_vline(xintercept = c(0, 1), colour = axis_col, linewidth = 0.35, linetype = "22") +
  geom_violin(scale = "width", width = 0.8, colour = NA, alpha = 0.9) +
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
            vjust = -2.0, size = 3.2, colour = ink_primary,
            fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_discrete(labels = function(x) {
    n <- three_stats$n_cells[match(x, as.character(three_stats$celltype))]
    paste0(str_remove(x, "^Glut-"), "\n(n = ", n, ")")
  }) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  annotate("text", x = 0, y = 3.58, label = "DACH2-1 anchor (0)", hjust = 0.5, vjust = 0,
           size = 2.8, colour = ink_muted) +
  annotate("text", x = 1, y = 3.58, label = "HVCra anchor (1)", hjust = 0.5, vjust = 0,
           size = 2.8, colour = ink_muted) +
  coord_cartesian(ylim = c(0.62, 3.78), clip = "off") +
  labs(title = "Glut-DACH2-HVCra-Int sits between its two anchors",
       subtitle = paste0("Per-cell score on ", 2 * N_AXIS_GENES,
                         " genes differential between Glut-DACH2-1 and Glut-DACH2-HVCra alone,\n",
                         "rescaled so the anchor medians are 0 and 1 — so the anchors are fixed by definition and\n",
                         "only HVCra-Int's position is a measurement. Point and bar are its median with a bootstrap\n",
                         "95% CI (2000 resamples); box is the IQR."),
       x = "Position on the DACH2-1 → HVCra axis", y = NULL) +
  theme_viz(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(colour = ink_primary, face = "bold", size = 9.5),
        plot.margin = margin(5.5, 12, 5.5, 5.5))

ggsave(file.path(out_dir, "hvcra_int_three_way_violin.pdf"), three_violin,
       width = 7.4, height = 4.2, device = cairo_pdf)
ggsave(file.path(out_dir, "hvcra_int_three_way_violin.png"), three_violin,
       width = 7.4, height = 4.2, dpi = 300, bg = surface)
write_csv(three_stats %>% mutate(across(median:iqr_width, ~round(.x, 4))),
          file.path(out_dir, "three_way_axis_medians.csv"))

message("three-way medians [95% CI]:")
print(as.data.frame(three_stats %>% mutate(across(median:iqr_width, ~round(.x, 3)))))

message("slope (fraction of the way HVCra-Int sits): ", round(fit, 3))
message("Wrote ", out_dir)
