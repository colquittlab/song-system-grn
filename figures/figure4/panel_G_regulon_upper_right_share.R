#!/usr/bin/env Rscript
# Figure 4, panel G -- share of each eRegulon's targets that are higher in both MAFB on/off contrasts ("upper-right quadrant"), across eRegulons.
#
# Contrasts: RA vs C1H-1 and PVALB-1 vs LGE-1. Genes: pseudobulk DESeq2 log2 fold change (as in panel F). Regions: ArchR accessibility log2 fold change
# (bias-matched Wilcoxon on the consensus peak matrix; markers_<FG>_VS_<BG>.rds). For every +/+ eRegulon (config37) with >= MIN_REGIONS linked regions:
#   share   fraction of its target genes (resp. linked regions) with log2FC > 0 in both contrasts, among those tested in both
#   p       one-sided permutation p-value of the share: each eRegulon against B random sets of the same size drawn from the pooled targets (resp. regions)
#           of all these eRegulons (10,000 for genes, 2,000 for regions; seed 20260 per kind)
# Plot (compact, flipped): one point per eRegulon, share on x, genes and regions as two strips; point area = number of targets (log scale), fill = -log10 p
# (single-hue); MAFB with a thick outline; eRegulons with a larger share than MAFB labeled. 2.6 x 1.6 in.
# Same computation as grn/analysis/scenicplus_regulon_quadrant_fractions.R + scenicplus_regulon_quadrant_permutation.R +
# scenicplus_regulon_quadrant_plot_size_n_color_p_flipped.R (all with PV_COMPARATOR=LGE1), which wrote
#   <scenicplus output>/config37/regulon_scatter_plots/PVALB1_vs_LGE1/quadrant_percentages_by_regulon_size_n_color_p_flipped.pdf
# To reproduce that file exactly, gene log2FCs are rounded to 3 decimals and the share and p-value to 4 significant digits, as they were in the
# intermediate tables it was drawn from.
#
# Inputs: the hybrid-labeled Seurat object, the config37 SCENIC+ MuData and the ArchR marker tests (common.R).
# Output: FIG4_OUT/panel_G_quadrant_share_size_n_color_p_flipped.{pdf,png,svg} and the per-eRegulon table (csv)
#
#   Rscript figures/figure4/panel_G_regulon_upper_right_share.R
source(here::here("figures/figure4/common.R"))
MIN_REGIONS <- 100; B_GENE <- 10000; B_REGION <- 2000
fig_check_font()

## inputs: log2 fold changes in the two contrasts
rna <- mafb_onoff_rna_contrasts() %>% filter(!is.na(lfc_RA), !is.na(lfc_PV1)) %>% mutate(lfc_RA = round(lfc_RA, 3), lfc_PV1 = round(lfc_PV1, 3))
gx <- setNames(rna$lfc_RA, rna$gene); gy <- setNames(rna$lfc_PV1, rna$gene)
LFC_RA <- archr_lfc(RA, C1H1); LFC_PV <- archr_lfc(PV1, LGE1)
e <- read_eregulons() %>% distinct(TF, Gene, Region) %>% filter(Gene %in% rna$gene, Region %in% names(LFC_RA), Region %in% names(LFC_PV))
tfs <- e %>% group_by(TF) %>% filter(n_distinct(Region) >= MIN_REGIONS) %>% pull(TF) %>% unique() %>% sort()
gsets <- lapply(setNames(tfs, tfs), function(tf) unique(e$Gene[e$TF == tf])); rsets <- lapply(setNames(tfs, tfs), function(tf) unique(e$Region[e$TF == tf]))
gpool <- sort(unique(unlist(gsets))); rpool <- sort(unique(unlist(rsets)))
cat(sprintf("%d eRegulons; gene pool %d, region pool %d\n", length(tfs), length(gpool), length(rpool)))

## share in the upper-right quadrant and its permutation p-value
share_of <- function(a, b) mean(a & b)
perm_kind <- function(sets, xa, yb, B, kind) {
  pool <- names(xa); A <- xa > 0; Bp <- yb > 0
  set.seed(20260)
  map_dfr(names(sets), function(tf) {
    idx <- match(sets[[tf]], pool); n <- length(idx)
    obs <- share_of(A[idx], Bp[idx])
    nul <- vapply(seq_len(B), function(i) { s <- sample.int(length(pool), n); share_of(A[s], Bp[s]) }, numeric(1))
    tibble(TF = tf, kind = kind, n = n, frac_UR = obs, null_frac_mean = mean(nul), p_frac = (1 + sum(nul >= obs)) / (1 + B))
  })
}
R <- bind_rows(perm_kind(gsets, gx[gpool], gy[gpool], B_GENE, "gene_expression"), perm_kind(rsets, LFC_RA[rpool], LFC_PV[rpool], B_REGION, "region_accessibility")) %>%
  mutate(across(c(frac_UR, null_frac_mean, p_frac), ~ signif(.x, 4)))
R <- R %>% arrange(match(TF, tfs), kind)      # row order of the original table (per eRegulon: genes, then regions); the jitter below depends on it
write_csv(R, file.path(FIG4_OUT, "panel_G_quadrant_share.csv"))
m <- R %>% filter(TF == "MAFB")
cat(sprintf("MAFB: genes %.1f%% (p = %.2g), regions %.1f%% (p = %.2g)\n", 100 * m$frac_UR[1], m$p_frac[1], 100 * m$frac_UR[2], m$p_frac[2]))

## plot
KIND <- c(gene_expression = "genes", region_accessibility = "regions")
YPOS <- c(genes = 1.55, regions = 1)
R <- R %>% select(TF, kind, n, frac_UR, p_frac) %>% mutate(pct = 100 * frac_UR, nlp = -log10(p_frac), kind = factor(KIND[kind], levels = KIND), yi = YPOS[as.character(kind)]) %>%
  group_by(kind) %>% mutate(rank = rank(-pct, ties.method = "min")) %>% ungroup()
set.seed(2026)
R <- R %>% mutate(yj = ifelse(TF == "MAFB", yi, yi + runif(n(), -0.16, 0.16)))
M <- R %>% filter(TF == "MAFB") %>% mutate(lab = sprintf("MAFB #%d, %.1f%%", rank, pct))
A <- R %>% left_join(M %>% select(kind, pct_mafb = pct), by = "kind") %>% filter(TF != "MAFB", pct > pct_mafb)
# labeled dots (eRegulons above MAFB, then MAFB) staggered in y in share order so no leader line crosses another labeled dot
ST <- bind_rows(A, M %>% select(any_of(names(A)))) %>% group_by(kind) %>% arrange(desc(pct), .by_group = TRUE) %>% mutate(yj = yi + 0.17 - 0.055 * (row_number() - 1)) %>% ungroup()
A <- ST %>% filter(TF != "MAFB"); M <- M %>% select(-yj) %>% left_join(ST %>% filter(TF == "MAFB") %>% select(kind, yj), by = "kind")
LX <- 56; GAP <- 0.095
LAB <- bind_rows(A %>% mutate(lab = TF), M) %>% arrange(kind, desc(pct)) %>% group_by(kind) %>%
  mutate(ly = { y <- numeric(n()); for (k in seq_len(n())) y[k] <- if (k == 1) yj[k] + (if (kind[1] == "genes") 0.06 else -0.05) else min(yj[k], y[k - 1] - GAP); y }) %>% ungroup()
OUTLINE <- 0.12
p <- ggplot(R, aes(pct, yj)) +
  geom_point(data = R %>% filter(TF != "MAFB", !paste(TF, kind) %in% paste(A$TF, A$kind)), aes(size = n, fill = nlp), shape = 21, color = FIG_INK_SECONDARY, stroke = OUTLINE, alpha = 0.9) +
  geom_segment(data = LAB, aes(x = LX - 0.8, xend = pct, y = ly, yend = yj), linewidth = 0.2, color = FIG_INK_PRIMARY) +
  geom_point(data = A, aes(size = n, fill = nlp), shape = 21, color = FIG_INK_SECONDARY, stroke = OUTLINE) +
  geom_point(data = M, aes(size = n, fill = nlp), shape = 21, color = FIG_INK_PRIMARY, stroke = 0.5) +
  geom_text(data = LAB, aes(x = LX, y = ly, label = lab), hjust = 0, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY) +
  scale_size_continuous(name = "targets", range = c(0.4, 2.2), trans = "log10", breaks = c(100, 1000, 10000)) +
  scale_fill_gradient(name = expression(-log[10]~p), low = "#e6eef9", high = "#0f3c78", limits = c(0, NA), breaks = c(0, 2, 4)) +
  guides(size = guide_legend(order = 2, title.position = "left", nrow = 1, override.aes = list(fill = FIG_NEUTRAL_FILL, color = FIG_INK_SECONDARY, stroke = OUTLINE)),
         fill = guide_colorbar(order = 1, title.position = "left", direction = "horizontal", barwidth = unit(0.5, "in"), barheight = unit(0.07, "in"), ticks.linewidth = 0.2)) +
  scale_x_continuous(breaks = seq(10, 50, 10), limits = c(5, 77), expand = c(0, 0)) +
  scale_y_continuous(breaks = unname(YPOS), labels = names(YPOS), limits = c(0.82, 1.86), expand = c(0, 0)) +
  labs(x = "% in upper-right quadrant", y = NULL) + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.11, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.03, "in"), legend.spacing.x = unit(0.04, "in"),
        legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal",
        legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
save_gg(p, "panel_G_quadrant_share_size_n_color_p_flipped", width = 2.6, height = 1.6)
