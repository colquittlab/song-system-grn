#!/usr/bin/env Rscript
# The MAFB gene and region differential scatters with PVALB-1 vs GABA-LGE-1 as the interneuron contrast (instead of PVALB-1 vs LAMP5), for comparison with
# scenicplus_mafb_onoff_scatter.R and scenicplus_mafb_region_dar_scatter.R, whose LAMP5-based outputs are left untouched; everything here goes to
# <config37>/PVALB1_vs_LGE1/ with the same file names. Same plot style, labeling rules and viridis TF2G-importance color as the originals.
#   x = RA vs C1H-1 log2FC     y = PVALB-1 vs LGE-1 log2FC
# Gene scatter:   one point per MAFB +/+ target (config37) tested in both contrasts (pseudobulk DESeq2: scenicplus_all_genes_onoff_contrasts_config37.csv.gz for
#                 RA vs C1H-1, scenicplus_pv1_vs_lge1_contrast.R for PVALB-1 vs LGE-1); the original variants: all genes, and genes with |log2FC| > 0.5 / > 1
#                 in at least one contrast.
# Region scatter: the regions linked to the targets with log2FC > UP_LFC in BOTH contrasts (directional rule, no padj, the same rule that defined the
#                 gene set behind the original region scatter), ArchR accessibility log2FC for RA vs C1H-1 and PVALB-1 vs LGE-1 (make_song_pair_dars_archr.R and
#                 make_interneuron_dars_archr.R), whatever their significance.
#
#   Rscript scenicplus_mafb_pvalb1_vs_lge1_scatters.R
suppressMessages({library(SummarizedExperiment); library(tidyverse); library(ggrepel); library(egg); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
ARCHR <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
DIR <- file.path(OUT, MAIN, "PVALB1_vs_LGE1"); dir.create(DIR, showWarnings = FALSE)
PANEL_IN <- 1.6; UP_LFC <- 1; N_LABEL <- 10
LABEL_KEEP <- c("KCNC1", "ERBB4", "PVALB")
LIT_GENES <- c("KCNC1", "KCNC2", "PVALB", "ERBB4", "MAF", "MAFB", "SOX6", "CNTNAP2")     # literature support for fast-spiking / PV-class roles (pv_fs_literature_genes.csv)
LABEL_MAP <- c(LOC110468166 = "NTM", LOC110470685 = "TLE4", LOC110474504 = "VIPR1", LOC116183613 = "ncRNA-1", LOC110482360 = "LAMA3")
nm <- function(g) ifelse(g %in% names(LABEL_MAP), LABEL_MAP[g], g)
YLAB_G <- "log2 fold change\nPVALB-1 vs LGE-1"; YLAB_R <- "log2FC accessibility\nPVALB-1 vs LGE-1"

# ---- targets and importance -------------------------------------------------------------------------------------------------------------------------
E <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c", importance_TF2G = "d"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name))
imp <- E %>% group_by(gene = Gene) %>% summarise(importance_TF2G = max(importance_TF2G), .groups = "drop") %>% mutate(importance_log_TF2G = log(importance_TF2G))
ra <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE) %>% select(gene, lfc_RA_vs_C1H1)
pv <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_pv1_vs_lge1_", MAIN, ".csv.gz")), show_col_types = FALSE) %>% select(gene, lfc_PV1_vs_LGE1)
W <- imp %>% inner_join(ra, by = "gene") %>% inner_join(pv, by = "gene") %>% filter(!is.na(lfc_RA_vs_C1H1), !is.na(lfc_PV1_vs_LGE1))
cat("MAFB +/+ targets:", n_distinct(E$Gene), "| tested in both contrasts:", nrow(W), "\n")
write_csv(W %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_mafb_pvalb1_vs_lge1_targets_", MAIN, ".csv")))

theme_sc <- function() theme_fig() + theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
                                           legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
save_sc <- function(p, stem) {
  p <- egg::set_panel_size(p, width = unit(PANEL_IN, "in"), height = unit(PANEL_IN, "in"))
  Wd <- PANEL_IN + 1.75; Ht <- PANEL_IN + 0.85
  fig_save(p, stem, width = Wd, height = Ht)
  ggsave(paste0(stem, ".png"), p, width = Wd, height = Ht, dpi = 600, bg = "transparent")
  ggsave(paste0(stem, ".svg"), p, width = Wd, height = Ht, bg = "transparent",
         device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
}
repel_layer <- function(lab) geom_text_repel(data = lab, aes(label = label), xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                                             segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, point.padding = 0.05, force = 3, max.time = 15,
                                             max.iter = 300000, max.overlaps = Inf, seed = 1)

# ---- gene scatter -----------------------------------------------------------------------------------------------------------------------------------
gene_scatter <- function(min_abs = NULL) {
  d <- W %>% transmute(gene, x = lfc_RA_vs_C1H1, y = lfc_PV1_vs_LGE1, importance_log_TF2G, both = x > UP_LFC & y > UP_LFC)
  tag <- "PVALB1_vs_LGE1"
  if (!is.null(min_abs)) { d <- d %>% filter(abs(x) > min_abs | abs(y) > min_abs); tag <- paste0(tag, "_absLFC", min_abs) }
  r <- cor(d$x, d$y); rho <- cor(d$x, d$y, method = "spearman")
  cand <- d %>% filter(both) %>% mutate(weak = pmin(x, y)) %>% arrange(desc(weak))
  labeled <- union(union(head(cand$gene, N_LABEL), intersect(LABEL_KEEP, d$gene)), intersect(LIT_GENES, d$gene))
  d <- d %>% mutate(label = nm(gene))
  cat(sprintf("%s: n=%d, Pearson r=%.3f, Spearman rho=%.3f | up in both (log2FC > %g): %d\n  labeled: %s\n", tag, nrow(d), r, rho, UP_LFC, nrow(cand), paste(nm(labeled), collapse = ", ")))
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = importance_log_TF2G), size = 1, alpha = 0.9, stroke = 0) +
    repel_layer(d %>% filter(gene %in% labeled)) +
    annotate("text", x = min(d$x), y = max(d$y), label = paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(d), " genes"), hjust = 0, vjust = 1,
             size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
    scale_color_viridis_c(name = "log TF2G\nimportance") + scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
    labs(x = "log2 fold change\nRA vs C1H-1", y = YLAB_G) + theme_sc()
  save_sc(p, file.path(DIR, paste0("mafb_onoff_scatter_", tag)))
  invisible(cand$gene)
}
up_both <- gene_scatter(); gene_scatter(0.5); gene_scatter(1)

# ---- region scatter ---------------------------------------------------------------------------------------------------------------------------------
get_lfc <- function(rds) { se <- readRDS(file.path(ARCHR, rds)); rr <- rowData(se); tibble(region = paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end), lfc = as.numeric(assay(se, "Log2FC")[, 1])) }
LFC_RA <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds") %>% rename(lfc_RA = lfc)
LFC_PV <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-LGE-1.rds") %>% rename(lfc_PV1 = lfc)
gimp <- setNames(imp$importance_log_TF2G, imp$gene)
R <- E %>% filter(Gene %in% up_both) %>% distinct(Gene, Region) %>% group_by(region = Region) %>%
  summarise(genes = paste(sort(unique(Gene)), collapse = ";"), imp = max(gimp[Gene]), .groups = "drop") %>%
  left_join(LFC_RA, by = "region") %>% left_join(LFC_PV, by = "region")
cat("up-in-both genes:", length(up_both), "| linked regions:", nrow(R), "| missing from the ArchR peak set:", sum(is.na(R$lfc_RA) | is.na(R$lfc_PV1)), "\n")
R <- R %>% filter(!is.na(lfc_RA), !is.na(lfc_PV1)) %>% mutate(label = sapply(strsplit(genes, ";"), function(g) paste(nm(g), collapse = "/")), weak = pmin(lfc_RA, lfc_PV1))
write_csv(R %>% mutate(across(c(lfc_RA, lfc_PV1, imp), ~ round(.x, 3))), file.path(HERE, paste0("scenicplus_mafb_up_both_regions_lge1_differential_accessibility_", MAIN, ".csv")))
r <- cor(R$lfc_RA, R$lfc_PV1); rho <- cor(R$lfc_RA, R$lfc_PV1, method = "spearman")
cat(sprintf("regions: %d | Pearson r = %.3f, Spearman rho = %.3f\n", nrow(R), r, rho))
top <- R %>% filter(lfc_RA > UP_LFC, lfc_PV1 > UP_LFC) %>% arrange(desc(weak)) %>% head(N_LABEL)
keep <- R %>% filter(sapply(strsplit(genes, ";"), function(g) any(g %in% LABEL_KEEP))) %>% mutate(fg = sapply(strsplit(genes, ";"), function(g) LABEL_KEEP[LABEL_KEEP %in% g][1])) %>%
  group_by(fg) %>% slice_max(weak, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(label = fg)
lab <- bind_rows(top, keep) %>% distinct(region, .keep_all = TRUE)
cat("labeled regions:", paste(lab$label, collapse = ", "), "\n")
p <- ggplot(R, aes(lfc_RA, lfc_PV1)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(aes(color = imp), size = 1, alpha = 0.9, stroke = 0) + repel_layer(lab) +
  annotate("text", x = min(R$lfc_RA), y = max(R$lfc_PV1), label = paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(R), " regions"), hjust = 0, vjust = 1,
           size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
  scale_color_viridis_c(name = "log TF2G\nimportance") + scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
  labs(x = "log2FC accessibility\nRA vs C1H-1", y = YLAB_R) + theme_sc()
save_sc(p, file.path(DIR, "mafb_up_both_regions_dar_scatter"))
cat("saved to", DIR, "\n")
