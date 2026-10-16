#!/usr/bin/env Rscript
# A defined way to pull regulon sets out of the panel-A heatmap (positive gene-AUC, config37), instead of cutting the row dendrogram by eye.
#
# Question: regulons whose target-gene expression is enriched in HVCra, HVCx and/or RA relative to D2-1, D2-4 and C1H-1, taken as the three matched contrasts
# used for the song-pair DARs:  HVCra vs D2-1,  HVCx vs D2-4,  RA vs C1H-1.
#
# 1. Per contrast, per regulon, with birds as the replicates (one bird in one library = sample_short x assignment, >= MIN_CELLS cells of the type): the per-bird
#    difference in AUC between the two types, on the heatmap's 0-1 scale (a regulon's AUC divided by its range across the 27 cell-type means). A contrast is
#    "differential" for a regulon when
#      mean per-bird difference >= DIFF_MIN,  the difference is > 0 in >= POS_FRAC of the paired birds (>= MIN_BIRDS of them),  and the foreground type is active
#      (its heatmap value, the 0-1 cluster mean, is >= FG_ACTIVE: not merely less active in the reference).
#    This is an effect size plus replication across birds, not a test: with 5-7 paired birds per contrast the smallest attainable p-value (2^-n) cannot survive
#    correction over 123 regulons. (A permutation of FG/REF labels across glutamatergic types was also tried and has no power by construction.)
# 2. A regulon is selected if at least one contrast is differential; its module is the pattern of contrasts it passes (all three; HVCra & HVCx; RA only; ...).
# 3. The MAFB module: regulons in MAFB's branch of an average-linkage tree on 1 - Pearson r of the AUC profiles over ALL 27 cell types, cut at MAFB_CUT.
# Output: a csv with the per-contrast numbers, pass/fail and module of every regulon, and the expanded heatmap (cells square) with the three differences alongside.
#
# Input: scenicplus_gene_auc_cluster_means.R (cell-type and per-bird AUC means).
#   Rscript scenicplus_gene_auc_contrast_clusters.R
suppressMessages({library(tidyverse); library(ComplexHeatmap); library(circlize); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
DIFF_MIN <- 0.15; POS_FRAC <- 0.8; MIN_BIRDS <- 3; MIN_CELLS <- 10; FG_ACTIVE <- 0.4; MAFB_CUT <- 0.25; CELL_IN <- 0.12
PAIRS <- list(`HVCra vs D2-1` = c("Glut-DACH2-HVCra", "Glut-DACH2-1"), `HVCx vs D2-4` = c("Glut-DACH2-HVCx", "Glut-DACH2-4"), `RA vs C1H-1` = c("Glut-CACNA1H-RA", "Glut-CACNA1H-1"))
FGT <- sapply(PAIRS, `[`, 1); REFT <- sapply(PAIRS, `[`, 2)
PATTERN <- c(`111` = "all three contrasts", `110` = "HVCra & HVCx", `101` = "HVCra & RA", `011` = "HVCx & RA", `100` = "HVCra only", `010` = "HVCx only", `001` = "RA only")
S <- read_csv(file.path(HERE, paste0("scenicplus_gene_auc_cluster_means_scaled_", MAIN, ".csv")), show_col_types = FALSE)
A <- as.matrix(S[, -1]); rownames(A) <- S$tf
M <- read_csv(file.path(HERE, paste0("scenicplus_gene_auc_cluster_means_", MAIN, ".csv")), show_col_types = FALSE)
rng <- tibble(tf = M$tf, rng = apply(as.matrix(M[, -1]), 1, function(x) max(x) - min(x)))
# column order of the panel-A heatmap (labels_order in scenicplus_hybrid_all_template.qmd)
cl <- colnames(A)
ORDER <- unique(c(sort(grep("HVC", cl, value = TRUE)), sort(grep("DACH2-[0-9]", cl, value = TRUE)), sort(grep("CACNA1H-RA", cl, value = TRUE)), sort(grep("CACNA1H-[0-9]", cl, value = TRUE)),
                  intersect(c("Glut-Im", "Glut-NB", "Glut-NSC"), cl), sort(setdiff(grep("GABA", cl, value = TRUE), "GABA-Im")), intersect("GABA-Im", cl), sort(grep("Astro", cl, value = TRUE)),
                  sort(grep("Epen", cl, value = TRUE)), sort(grep("Oligo", cl, value = TRUE)), sort(grep("OPC", cl, value = TRUE)), sort(grep("Micro", cl, value = TRUE)), sort(grep("Endo", cl, value = TRUE))))
stopifnot(setequal(ORDER, cl))

## 1. per-bird contrasts
b <- read_csv(file.path(OUT, MAIN, "gene_auc_expanded", "scenicplus_gene_auc_per_bird.csv.gz"), show_col_types = FALSE) %>% filter(n_cells >= MIN_CELLS) %>%
  left_join(rng, by = "tf") %>% mutate(bird = paste(sample_short, assignment, sep = "|"), v = value_mean / rng)
ct <- imap_dfr(PAIRS, function(p, k) {
  x <- b %>% filter(cluster == p[1]) %>% select(tf, bird, va = v); y <- b %>% filter(cluster == p[2]) %>% select(tf, bird, vc = v)
  inner_join(x, y, by = c("tf", "bird")) %>% mutate(d = va - vc) %>% group_by(tf) %>% summarize(n_birds = n(), mean_diff = mean(d), frac_pos = mean(d > 0), .groups = "drop") %>% mutate(contrast = k)
}) %>% mutate(fg_active = A[cbind(tf, FGT[contrast])], pass = n_birds >= MIN_BIRDS & mean_diff >= DIFF_MIN & frac_pos >= POS_FRAC & fg_active >= FG_ACTIVE)
W <- ct %>% select(tf, contrast, mean_diff, frac_pos, n_birds, pass) %>% pivot_wider(names_from = contrast, values_from = c(mean_diff, frac_pos, n_birds, pass), names_glue = "{.value}.{contrast}")
passmat <- as.matrix(W[, paste0("pass.", names(PAIRS))]); rownames(passmat) <- W$tf
pat <- apply(passmat, 1, function(r) paste0(as.integer(r), collapse = ""))
sc <- W %>% mutate(n_pass = rowSums(passmat), pattern = pat[tf], module = ifelse(n_pass >= 1, unname(PATTERN[pat[tf]]), NA)) %>% arrange(desc(n_pass), desc(rowMeans(across(starts_with("mean_diff")))))
## 3. the MAFB module
hm <- hclust(as.dist(1 - cor(t(A))), "average")
cm <- cutree(hm, h = MAFB_CUT)
mafb_set <- names(cm)[cm == cm["MAFB"]]
sc$mafb_module <- sc$tf %in% mafb_set
write_csv(sc %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_gene_auc_contrast_modules_", MAIN, ".csv")))
cat(sprintf("differential in >= 1 contrast: %d regulons (DIFF_MIN %g, >= %g%% of paired birds positive, FG active >= %g)\n", sum(sc$n_pass >= 1), DIFF_MIN, 100 * POS_FRAC, FG_ACTIVE))
for (m in PATTERN) if (any(sc$module == m, na.rm = TRUE)) cat(sprintf("  %-20s %s\n", m, paste(sc$tf[which(sc$module == m)], collapse = ", ")))
cat("MAFB branch at distance", MAFB_CUT, ":", paste(mafb_set, collapse = ", "), "\n")
chk <- c("AR", "FOSL2", "KLF6", "ETV1", "NR4A2", "LHX9", "TBR1")
cat("checks:\n"); print(as.data.frame(sc %>% filter(tf %in% chk) %>% transmute(tf, d_HVCra = round(`mean_diff.HVCra vs D2-1`, 2), d_HVCx = round(`mean_diff.HVCx vs D2-4`, 2), d_RA = round(`mean_diff.RA vs C1H-1`, 2), pattern, module)), row.names = FALSE)

## expanded heatmap: rows grouped by module, cells square, the three differences alongside
grp_levels <- c(unname(PATTERN[c("111", "110", "101", "011", "100", "010", "001")]), "MAFB module")
keep <- sc %>% filter(!is.na(module) | mafb_module) %>% mutate(group = ifelse(mafb_module, "MAFB module", module)) %>% filter(group %in% grp_levels)
keep <- keep %>% mutate(group = factor(group, levels = grp_levels)) %>% droplevels() %>% arrange(group, desc(rowSums(across(starts_with("mean_diff")))))
mat <- A[keep$tf, ORDER]
D <- as.matrix(keep[, paste0("mean_diff.", names(PAIRS))]); colnames(D) <- names(PAIRS); rownames(D) <- keep$tf
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)
GAP <- 0.05
nsp <- nlevels(keep$group)
body_h <- unit(nrow(mat) * CELL_IN + (nsp - 1) * GAP, "in")
cols <- colorRamp2(seq(0, 1, length.out = 9), scales::brewer_pal(palette = "Greys", direction = 1)(9))
fgmark <- HeatmapAnnotation(type = anno_simple(ifelse(colnames(mat) %in% FGT, "foreground", ifelse(colnames(mat) %in% REFT, "reference", "other")),
                                               col = c(foreground = FIG_PAL[2], reference = FIG_PAL[1], other = "white"), height = unit(0.06, "in")), show_annotation_name = FALSE, show_legend = FALSE)
ht <- Heatmap(mat, name = "target gene enrichment", col = cols, cluster_rows = FALSE, cluster_columns = FALSE, row_split = keep$group, row_title_gp = gp6, row_title_rot = 0,
              row_names_gp = gp5, column_names_gp = gp5, row_gap = unit(GAP, "in"), top_annotation = fgmark, show_row_names = FALSE,
              width = unit(ncol(mat) * CELL_IN, "in"), height = body_h, border = FALSE,
              heatmap_legend_param = list(title = "Target gene\nenrichment", title_gp = gp6, labels_gp = gp5, legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in"), at = c(0, 1), labels = c("0", "1")))
DC <- 0.6
htd <- Heatmap(D, name = "per-bird difference", col = colorRamp2(c(-DC, 0, DC), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH)), cluster_rows = FALSE, cluster_columns = FALSE, row_split = keep$group, row_title = NULL,
               row_names_gp = gp5, column_names_gp = gp5, row_gap = unit(GAP, "in"), show_row_names = TRUE, row_names_side = "right", width = unit(ncol(D) * CELL_IN, "in"), height = body_h, border = FALSE,
               cell_fun = function(j, i, x, y, w, h, fill) if (passmat[rownames(D)[i], j]) grid.points(x, y, pch = 16, size = unit(0.025, "in"), gp = gpar(col = FIG_INK_PRIMARY)),
               heatmap_legend_param = list(title = "Difference in\nscaled AUC\n(per-bird mean)", title_gp = gp6, labels_gp = gp5, legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in"), at = c(-DC, 0, DC), labels = c("≤−0.6", "0", "≥0.6")))
dir.create(file.path(OUT, MAIN, "gene_auc_expanded"), showWarnings = FALSE)
stem <- file.path(OUT, MAIN, "gene_auc_expanded", "gene_auc_contrast_modules")
Wd <- 7.4; Ht <- nrow(mat) * CELL_IN + nsp * GAP + 2.2
cairo_pdf(paste0(stem, ".pdf"), width = Wd, height = Ht, family = FIG_FONT); draw(ht + htd, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"), ht_gap = unit(0.06, "in")); invisible(dev.off())
png(paste0(stem, ".png"), width = Wd, height = Ht, units = "in", res = 300, type = "cairo", family = FIG_FONT); draw(ht + htd, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"), ht_gap = unit(0.06, "in")); invisible(dev.off())
cat("saved", stem, "\n")
