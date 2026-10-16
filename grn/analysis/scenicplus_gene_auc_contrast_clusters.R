#!/usr/bin/env Rscript
# A defined way to pull regulon sets out of the panel-A heatmap (positive gene-AUC, config37), instead of cutting the row dendrogram by eye.
#
# 1. Contrast score (supervised). Regulons whose target-gene expression is enriched in HVCra, HVCx and/or RA relative to D2-1, D2-4 and C1H-1:
#      Delta = max over FG (HVCra, HVCx, RA) of the row-scaled AUC (0-1 across cell types, as in the heatmap) - mean over REF (D2-1, D2-4, C1H-1).
#    A regulon is selected if Delta >= DELTA_MIN and its best FG value is >= FG_MIN (it has to be active, not merely less active in REF).
#    This is an effect-size criterion, not a test: an exact permutation of FG/REF labels over the 12 glutamatergic types (18,480 assignments) was tried and has no
#    power by construction (a regulon specific to one type is in the FG in 25% of assignments; AR, Delta = 0.78, gets q = 0.98), so thresholds are reported, not p-values.
# 2. Modules among the selected regulons (unsupervised, on the columns that define the question): Ward clustering of the three FG values (HVCra, HVCx, RA; 0-1), cut into
#    K_MODULES and named by the FG types where the module mean is within 25% of its maximum. Clustering the whole heatmap on all 27 cell types is dominated by the
#    GABAergic / non-neuronal columns, which is why these sets end up apart (and mixed with other glutamatergic regulons) in the full dendrogram.
# 3. The MAFB module: regulons in MAFB's branch of an average-linkage tree on 1 - Pearson r over ALL 27 cell types, cut at distance MAFB_CUT.
# Output: table of scores, module membership (csv), and the expanded heatmap of the three groups (config37/gene_auc_expanded/).
#
#   Rscript scenicplus_gene_auc_contrast_clusters.R
suppressMessages({library(tidyverse); library(ComplexHeatmap); library(circlize); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
DELTA_MIN <- 0.2; FG_MIN <- 0.6; K_MODULES <- 3; MAFB_CUT <- 0.25
FG <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-CACNA1H-RA"); REF <- c("Glut-DACH2-1", "Glut-DACH2-4", "Glut-CACNA1H-1")
GLUT <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4", "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2",
          "Glut-Im", "Glut-NB", "Glut-NSC")
S <- read_csv(file.path(HERE, paste0("scenicplus_gene_auc_cluster_means_scaled_", MAIN, ".csv")), show_col_types = FALSE)
A <- as.matrix(S[, -1]); rownames(A) <- S$tf
# column order of the panel-A heatmap (labels_order in scenicplus_hybrid_all_template.qmd)
cl <- colnames(A)
ORDER <- unique(c(sort(grep("HVC", cl, value = TRUE)), sort(grep("DACH2-[0-9]", cl, value = TRUE)), sort(grep("CACNA1H-RA", cl, value = TRUE)), sort(grep("CACNA1H-[0-9]", cl, value = TRUE)),
                  intersect(c("Glut-Im", "Glut-NB", "Glut-NSC"), cl), sort(setdiff(grep("GABA", cl, value = TRUE), "GABA-Im")), intersect("GABA-Im", cl), sort(grep("Astro", cl, value = TRUE)),
                  sort(grep("Epen", cl, value = TRUE)), sort(grep("Oligo", cl, value = TRUE)), sort(grep("OPC", cl, value = TRUE)), sort(grep("Micro", cl, value = TRUE)), sort(grep("Endo", cl, value = TRUE))))
stopifnot(setequal(ORDER, cl))

## 1. contrast score
sc <- tibble(tf = rownames(A), delta = apply(A, 1, function(x) max(x[FG]) - mean(x[REF])), fg_max = apply(A[, FG], 1, max), ref_mean = rowMeans(A[, REF]),
             best_fg = sub("Glut-", "", FG[apply(A[, FG], 1, which.max)])) %>% mutate(selected = delta >= DELTA_MIN & fg_max >= FG_MIN) %>% arrange(desc(delta))
## 2. modules among the selected regulons
sel <- sc$tf[sc$selected]
hc <- hclust(dist(A[sel, FG]), "ward.D2")
mod <- cutree(hc, K_MODULES)
short <- function(x) sub("Glut-", "", sub("DACH2-", "", sub("CACNA1H-", "", x)))
mod_name <- sapply(split(names(mod), mod), function(m) { v <- colMeans(A[m, FG, drop = FALSE]); paste0("enriched in ", paste(short(names(v)[v >= 0.75 * max(v)]), collapse = "/")) })
sc$module <- ifelse(sc$tf %in% sel, paste0("module ", mod[sc$tf], ": ", mod_name[as.character(mod[sc$tf])]), NA)
## 3. the MAFB module
hm <- hclust(as.dist(1 - cor(t(A))), "average")
cm <- cutree(hm, h = MAFB_CUT)
mafb_set <- names(cm)[cm == cm["MAFB"]]
sc$module <- ifelse(sc$tf %in% mafb_set, "MAFB module", sc$module)
write_csv(sc %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_gene_auc_contrast_modules_", MAIN, ".csv")))
cat(sprintf("selected: %d regulons (Delta >= %g, best FG >= %g); modules (K = %d):\n", length(sel), DELTA_MIN, FG_MIN, K_MODULES))
for (m in sort(unique(na.omit(sc$module)))) cat(sprintf("  %s: %s\n", m, paste(sc$tf[which(sc$module == m)], collapse = ", ")))
r <- cor(t(A))["MAFB", ]; cat("MAFB branch at distance", MAFB_CUT, ":", paste(mafb_set, collapse = ", "), "| profile r with MAFB:", paste(round(r[setdiff(mafb_set, "MAFB")], 2), collapse = " "), "\n")

## expanded heatmap
keep <- sc %>% filter(!is.na(module)) %>% mutate(module = factor(module, levels = c(sort(grep("^module", unique(module), value = TRUE)), "MAFB module")))
ord <- unlist(lapply(levels(keep$module), function(m) { t <- keep$tf[keep$module == m]; if (length(t) < 3) t else t[hclust(as.dist(1 - cor(t(A[t, ]))), "average")$order] }))
mat <- A[ord, ORDER]
split_f <- factor(keep$module[match(ord, keep$tf)], levels = levels(keep$module))
cols <- colorRamp2(seq(0, 1, length.out = 9), scales::brewer_pal(palette = "Greys", direction = 1)(9))
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)
fgmark <- HeatmapAnnotation(`FG` = anno_simple(ifelse(colnames(mat) %in% FG, "FG", ifelse(colnames(mat) %in% REF, "REF", "other")), col = c(FG = FIG_PAL[2], REF = FIG_PAL[1], other = "white"), height = unit(0.06, "in")),
                            annotation_name_gp = gp5, show_annotation_name = FALSE)
ht <- Heatmap(mat, name = "target gene\nenrichment", col = cols, cluster_rows = FALSE, cluster_columns = FALSE, row_split = split_f, row_title_gp = gp6, row_title_rot = 0,
              row_names_gp = gp5, column_names_gp = gp5, row_gap = unit(0.05, "in"), column_gap = unit(0, "in"), top_annotation = fgmark,
              width = unit(3.1, "in"), height = unit(0.085 * nrow(mat), "in"), border = FALSE,
              heatmap_legend_param = list(title_gp = gp6, labels_gp = gp5, legend_height = unit(0.6, "in"), grid_width = unit(0.09, "in"), at = c(0, 1), labels = c("0", "1")))
dir.create(file.path(OUT, MAIN, "gene_auc_expanded"), showWarnings = FALSE)
stem <- file.path(OUT, MAIN, "gene_auc_expanded", "gene_auc_contrast_modules")
H <- 0.085 * nrow(mat) + 2.4
cairo_pdf(paste0(stem, ".pdf"), width = 6.2, height = H, family = FIG_FONT); draw(ht, padding = unit(c(2, 2, 2, 2), "mm")); invisible(dev.off())
png(paste0(stem, ".png"), width = 6.2, height = H, units = "in", res = 300, type = "cairo", family = FIG_FONT); draw(ht, padding = unit(c(2, 2, 2, 2), "mm")); invisible(dev.off())
cat("saved", stem, "\n")
