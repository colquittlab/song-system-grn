#!/usr/bin/env Rscript
# Figure 4, panel B -- expansion of panel A: GRNs with differential activity between song neurons and their matched non-song neurons, and the MAFB module.
#
# 1. Activity: per-TF positive (+/+) eRegulon gene-based AUC, as in panel A (direct and extended pooled), averaged per cell type and scaled 0-1 across cell types
#    for display; and averaged per bird (one bird in one library, sample_short x assignment) and cell type, for the contrasts.
# 2. Three matched contrasts: HVCra vs D2-1, HVCx vs D2-4, RA vs C1H-1. Per regulon and contrast, the per-bird difference in AUC between the two cell types on
#    the panel-A scale (AUC / the regulon's range across the cell-type means), over the birds with >= MIN_CELLS cells of both types. A contrast is called when
#      mean per-bird difference >= DIFF_MIN, positive in >= POS_FRAC of the paired birds (>= MIN_BIRDS), and the song neuron's scaled AUC >= FG_ACTIVE.
#    (An effect-size and replication rule, not a significance test: with 5-7 paired birds the smallest attainable p-value cannot survive correction over 123 TFs.)
#    Regulons calling at least one contrast are grouped by the pattern of contrasts they call.
# 3. MAFB module: MAFB's branch of an average-linkage tree on 1 - Pearson r of the 0-1 cell-type profiles (all 27 types), cut at MAFB_CUT.
# Drawing: AUC heatmap WIDTH_IN wide with square cells, the three per-bird differences alongside, group titles above each block; no legends (the gray scale is
# panel A's; the difference scale is written separately as *_legend). Cropped to the drawing with scripts/crop_svg.py; PDF and PNG are made from the SVG.
# Same computation as grn/analysis/scenicplus_gene_auc_cluster_means.R + scenicplus_gene_auc_contrast_clusters.R, which wrote
#   <scenicplus output>/config37/gene_auc_expanded/gene_auc_contrast_modules.pdf
#
# Input: the config37 SCENIC+ MuData (common.R). Needs Inkscape and python3 (lxml) for the crop.
# Output: FIG4_OUT/panel_B_gene_auc_contrast_modules.{svg,pdf,png}, panel_B_..._legend.*, panel_B_gene_auc_contrast_modules.csv
#
#   Rscript figures/figure4/panel_B_grn_contrast_modules.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(ComplexHeatmap); library(circlize)})
DIFF_MIN <- 0.15; POS_FRAC <- 0.8; MIN_BIRDS <- 3; MIN_CELLS <- 10; FG_ACTIVE <- 0.4; MAFB_CUT <- 0.25; WIDTH_IN <- 1.9
PAIRS <- list(`HVCra vs D2-1` = c("Glut-DACH2-HVCra", "Glut-DACH2-1"), `HVCx vs D2-4` = c("Glut-DACH2-HVCx", "Glut-DACH2-4"), `RA vs C1H-1` = c(RA, C1H1))
FGT <- sapply(PAIRS, `[`, 1)
PATTERN <- c(`111` = "all three contrasts", `110` = "HVCra & HVCx", `101` = "HVCra & RA", `011` = "HVCx & RA", `100` = "HVCra only", `010` = "HVCx only", `001` = "RA only")

## 1. activity
auc <- read_gene_auc()
long <- positive_auc_long(auc, by = c("cluster", "sample_short", "assignment"))
Mraw <- long %>% group_by(tf, cluster) %>% summarize(v = mean(value), .groups = "drop")
rng <- Mraw %>% group_by(tf) %>% summarize(lo = min(v), rng = max(v) - min(v), .groups = "drop")
S <- Mraw %>% left_join(rng, by = "tf") %>% mutate(s = (v - lo) / rng) %>% select(tf, cluster, s) %>% pivot_wider(names_from = cluster, values_from = s)
A <- as.matrix(S[, -1]); rownames(A) <- S$tf
ORDER <- celltype_order(colnames(A))
ncell <- auc %>% count(cluster, sample_short, assignment, name = "n_cells")
bird <- long %>% group_by(tf, cluster, sample_short, assignment) %>% summarize(value_mean = mean(value), .groups = "drop") %>%
  left_join(ncell, by = c("cluster", "sample_short", "assignment")) %>% filter(n_cells >= MIN_CELLS) %>%
  left_join(rng, by = "tf") %>% mutate(bird = paste(sample_short, assignment, sep = "|"), v = value_mean / rng)

## 2. matched contrasts
ct <- imap_dfr(PAIRS, function(p, k) {
  x <- bird %>% filter(cluster == p[1]) %>% select(tf, bird, va = v); y <- bird %>% filter(cluster == p[2]) %>% select(tf, bird, vc = v)
  inner_join(x, y, by = c("tf", "bird")) %>% mutate(d = va - vc) %>% group_by(tf) %>% summarize(n_birds = n(), mean_diff = mean(d), frac_pos = mean(d > 0), .groups = "drop") %>% mutate(contrast = k)
}) %>% mutate(fg_active = A[cbind(tf, FGT[contrast])], pass = n_birds >= MIN_BIRDS & mean_diff >= DIFF_MIN & frac_pos >= POS_FRAC & fg_active >= FG_ACTIVE)
W <- ct %>% select(tf, contrast, mean_diff, frac_pos, n_birds, pass) %>% pivot_wider(names_from = contrast, values_from = c(mean_diff, frac_pos, n_birds, pass), names_glue = "{.value}.{contrast}")
passmat <- as.matrix(W[, paste0("pass.", names(PAIRS))]); rownames(passmat) <- W$tf
pat <- apply(passmat, 1, function(r) paste0(as.integer(r), collapse = ""))
sc <- W %>% mutate(n_pass = rowSums(passmat), pattern = pat[tf], module = ifelse(n_pass >= 1, unname(PATTERN[pat[tf]]), NA))

## 3. MAFB module
cm <- cutree(hclust(as.dist(1 - cor(t(A))), "average"), h = MAFB_CUT)
sc$mafb_module <- sc$tf %in% names(cm)[cm == cm["MAFB"]]
write_csv(sc %>% mutate(across(where(is.double), ~ round(.x, 3))) %>% arrange(desc(n_pass)), file.path(FIG4_OUT, "panel_B_gene_auc_contrast_modules.csv"))
for (m in PATTERN) if (any(sc$module == m, na.rm = TRUE)) cat(sprintf("  %-20s %s\n", m, paste(sc$tf[which(sc$module == m)], collapse = ", ")))
cat("  MAFB module          ", paste(sc$tf[sc$mafb_module], collapse = ", "), "\n")

## drawing
grp_levels <- c(unname(PATTERN[c("111", "110", "101", "011", "100", "010", "001")]), "MAFB module")
keep <- sc %>% filter(!is.na(module) | mafb_module) %>% mutate(group = factor(ifelse(mafb_module, "MAFB module", module), levels = grp_levels)) %>% droplevels() %>%
  arrange(group, desc(rowSums(across(starts_with("mean_diff")))))
mat <- A[keep$tf, ORDER]
D <- as.matrix(keep[, paste0("mean_diff.", names(PAIRS))]); colnames(D) <- names(PAIRS); rownames(D) <- keep$tf
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)
CELL_IN <- WIDTH_IN / ncol(mat); GAP <- 0.09; nsp <- nlevels(keep$group)
body_h <- unit(nrow(mat) * CELL_IN + (nsp - 1) * GAP, "in")
cols <- colorRamp2(seq(0, 1, length.out = 9), scales::brewer_pal(palette = "Greys", direction = 1)(9))
DC <- 0.6; dcol <- colorRamp2(c(-DC, 0, DC), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH))
ht <- Heatmap(mat, name = "target gene enrichment", col = cols, cluster_rows = FALSE, cluster_columns = FALSE, row_split = keep$group, row_title = NULL,
              column_names_gp = gp5, row_gap = unit(GAP, "in"), show_row_names = FALSE, show_heatmap_legend = FALSE,
              width = unit(ncol(mat) * CELL_IN, "in"), height = body_h, border = FALSE)
htd <- Heatmap(D, name = "per-bird difference", col = dcol, cluster_rows = FALSE, cluster_columns = FALSE, row_split = keep$group, row_title = NULL,
               row_names_gp = gp5, column_names_gp = gp5, row_gap = unit(GAP, "in"), show_row_names = TRUE, row_names_side = "right", show_heatmap_legend = FALSE,
               width = unit(ncol(D) * CELL_IN, "in"), height = body_h, border = FALSE)
draw_all <- function() {
  draw(ht + htd, padding = unit(c(1, 1, 1, 1), "mm"), ht_gap = unit(0.04, "in"))
  for (i in seq_len(nsp)) decorate_heatmap_body("target gene enrichment", grid.text(levels(keep$group)[i], x = unit(0, "npc"), y = unit(1, "npc") + unit(1.5, "pt"), just = c("left", "bottom"), gp = gp5), slice = i)
}
export <- function(draw_fun, stem, w = 6, h = 6) {          # draw on a generous transparent canvas, crop the SVG to the drawing, make PDF/PNG from it
  raw <- tempfile(fileext = ".svg")
  svglite::svglite(raw, width = w, height = h, bg = "transparent", system_fonts = list(sans = FIG_FONT)); draw_fun(); invisible(dev.off())
  size <- system2("python3", c(here::here("scripts/crop_svg.py"), raw, paste0(stem, ".svg")), stdout = TRUE)
  system2("inkscape", c(paste0(stem, ".svg"), "--export-type=pdf", paste0("--export-filename=", stem, ".pdf")), stdout = FALSE, stderr = FALSE)
  system2("inkscape", c(paste0(stem, ".svg"), "--export-type=png", "--export-dpi=600", "--export-background-opacity=0", paste0("--export-filename=", stem, ".png")), stdout = FALSE, stderr = FALSE)
  message("saved ", stem, " (", size, ")")
}
stem <- file.path(FIG4_OUT, "panel_B_gene_auc_contrast_modules")
export(draw_all, stem)
lgd <- Legend(title = "Difference in\nscaled AUC\n(per-bird mean)", col_fun = dcol, at = c(-DC, 0, DC), labels = c("≤−0.6", "0", "≥0.6"), title_gp = gp6, labels_gp = gp5,
              legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in"))
export(function() draw(lgd), paste0(stem, "_legend"), 2, 2)
