#!/usr/bin/env Rscript
# Differential accessibility (ArchR log2FC) of the regions behind the genes that are up in both RA and PVALB-1, as a two-column heatmap:
# RA vs C1H-1 and PVALB-1 vs LAMP5. All regions are shown whatever their significance; the left tracks mark which pass the DAR cutoff
# (FDR <= 0.05 and log2FC >= 0.585) in each contrast. Values come from the saved ArchR marker results of make_song_pair_dars_archr.R and
# make_interneuron_dars_archr.R (bias-matched Wilcoxon on the consensus peak matrix of the project copy).
#
# Rows are clustered with Euclidean distance (Ward). With two columns the correlation distance between rows is degenerate (every pair
# correlates at +1 or -1), so the correlation-based clustering used for the other heatmaps is not offered here.
#
#   Rscript scenicplus_mafb_region_differential_accessibility_heatmap.R
suppressMessages({library(SummarizedExperiment); library(ComplexHeatmap); library(circlize); library(tidyverse); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
STORE <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
d <- read_csv(here::here("grn/analysis/scenicplus_mafb_up_both_regions_accessibility_config37.csv"), show_col_types = FALSE) %>% select(region, genes, DAR_RA, DAR_PV1)

get_lfc <- function(rds) {
  se <- readRDS(file.path(STORE, rds))
  rr <- rowData(se)
  nm <- paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end)       # BED-style names, as in the consensus bed and the eRegulon tables
  tibble(region = nm, lfc = assay(se, "Log2FC")[, 1], fdr = assay(se, "FDR")[, 1])
}
ra <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds") %>% rename(lfc_RA = lfc, fdr_RA = fdr)
pv <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.rds") %>% rename(lfc_PV1 = lfc, fdr_PV1 = fdr)
x <- d %>% left_join(ra, by = "region") %>% left_join(pv, by = "region")
cat("regions:", nrow(x), "| not in the ArchR peak set (scaffolds outside its genome):", sum(is.na(x$lfc_RA) | is.na(x$lfc_PV1)), "\n")
x <- x %>% filter(!is.na(lfc_RA), !is.na(lfc_PV1))
write_csv(x %>% mutate(across(c(lfc_RA, lfc_PV1), ~ round(.x, 3)), across(c(fdr_RA, fdr_PV1), ~ signif(.x, 3))),
          here::here("grn/analysis/scenicplus_mafb_up_both_regions_differential_accessibility_config37.csv"))
cat(sprintf("Spearman rho of log2FC across regions: %.2f | log2FC > 0 in both: %d | > 0.585 in both: %d | DAR in RA: %d, in PVALB-1: %d, both: %d\n",
            cor(x$lfc_RA, x$lfc_PV1, method = "spearman"), sum(x$lfc_RA > 0 & x$lfc_PV1 > 0), sum(x$lfc_RA > 0.585 & x$lfc_PV1 > 0.585),
            sum(x$DAR_RA), sum(x$DAR_PV1), sum(x$DAR_RA & x$DAR_PV1)))
m <- cbind(`RA vs C1H-1` = x$lfc_RA, `PVALB-1 vs LAMP5` = x$lfc_PV1); rownames(m) <- x$region

FOCAL <- c("KCNC1", "PVALB", "ERBB4")
focal <- factor(sapply(strsplit(x$genes, ";"), function(g) { h <- FOCAL[FOCAL %in% g]; if (length(h)) h[1] else "other" }), levels = c(FOCAL, "other"))
dar_col <- c(`DAR` = FIG_INK_SECONDARY, `not DAR` = FIG_GRID)
focal_col <- setNames(c(FIG_PAL[1:3], FIG_NEUTRAL_FILL), levels(focal))
gp5 <- gpar(fontsize = FIG_PT_AXIS_TEXT, fontfamily = FIG_FONT); gp6 <- gpar(fontsize = FIG_PT_AXIS_TITLE, fontfamily = FIG_FONT)
LC <- 2
col_fun <- colorRamp2(c(-LC, 0, LC), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH))
ha <- rowAnnotation(`RA DAR` = factor(ifelse(x$DAR_RA == 1, "DAR", "not DAR"), levels = names(dar_col)),
                    `PVALB-1 DAR` = factor(ifelse(x$DAR_PV1 == 1, "DAR", "not DAR"), levels = names(dar_col)), `focal gene` = focal,
                    col = list(`RA DAR` = dar_col, `PVALB-1 DAR` = dar_col, `focal gene` = focal_col), annotation_name_gp = gp5, simple_anno_size = unit(0.09, "in"),
                    annotation_legend_param = list(title_gp = gp6, labels_gp = gp5, grid_height = unit(0.09, "in"), grid_width = unit(0.09, "in")))
hm <- Heatmap(m, name = "log2FC", col = col_fun, cluster_rows = hclust(dist(m), "ward.D"), cluster_columns = FALSE, show_row_names = FALSE,
              column_names_gp = gp5, row_dend_width = unit(0.35, "in"), left_annotation = ha, width = unit(0.5, "in"), height = unit(4, "in"),
              heatmap_legend_param = list(title = "Differential\naccessibility\n(log2FC)", title_gp = gp6, labels_gp = gp5, legend_height = unit(0.7, "in"),
                                          grid_width = unit(0.09, "in"), at = c(-LC, 0, LC), labels = c("≤−2", "0", "≥2")))
stem <- file.path(OUT, MAIN, "mafb_up_both_regions_differential_accessibility_heatmap")
cairo_pdf(paste0(stem, ".pdf"), width = 3.9, height = 6.4, family = FIG_FONT)
draw(hm, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
invisible(dev.off())
png(paste0(stem, ".png"), width = 3.9, height = 6.4, units = "in", res = 600, type = "cairo", family = FIG_FONT)
draw(hm, merge_legend = TRUE, padding = unit(c(2, 2, 2, 2), "mm"))
invisible(dev.off())
cat("saved", stem, "\n")
