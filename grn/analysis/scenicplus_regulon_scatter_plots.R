#!/usr/bin/env Rscript
# One scatter per regulon in the style of the MAFB DEG / region scatters (scenicplus_mafb_onoff_scatter.R, scenicplus_mafb_region_dar_scatter.R):
# x = log2 fold change RA vs C1H-1, y = log2 fold change PVALB-1 vs LAMP5, for the config37 +/+ regulons of the cross-regulon comparisons (TFs with
# >= MIN_REGIONS regions among tested genes). Two versions, in separate subfolders:
#   gene_expression/       one point per target gene, RNA log2FC (pseudobulk DESeq2; scenicplus_mafb_pvalb_vs_mafb_low_mge.R)
#   region_accessibility/  one point per linked region, ArchR accessibility log2FC (make_song_pair_dars_archr.R, make_interneuron_dars_archr.R)
# Color: log TF2G importance of the TF -> gene link (the highest of the linked genes for a region), viridis. Labels: up to N_LABEL points with the
# largest smaller-of-the-two log2FC, both above UP_LFC. Annotated with Pearson r, Spearman rho and n. 1.6 x 1.6 in panel, transparent background;
# each plot is written as PDF and as SVG (one text element per label).
#
#   Rscript scenicplus_regulon_scatter_plots.R
suppressMessages({library(SummarizedExperiment); library(tidyverse); library(ggrepel); library(egg); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
MIN_REGIONS <- 100; UP_LFC <- 1; N_LABEL <- 5; PANEL_IN <- 1.6
ARCHR <- "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/archr_consensus"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
D <- file.path(OUT, MAIN, "regulon_scatter_plots")
for (v in c("gene_expression", "region_accessibility")) dir.create(file.path(D, v), showWarnings = FALSE, recursive = TRUE)
LABEL_MAP <- c(LOC110468166 = "NTM", LOC110470685 = "TLE4", LOC110474504 = "VIPR1", LOC116183613 = "ncRNA-1")
nm <- function(x) ifelse(x %in% names(LABEL_MAP), LABEL_MAP[x], x)

get_lfc <- function(rds) { se <- readRDS(file.path(ARCHR, rds)); rr <- rowData(se); setNames(as.numeric(assay(se, "Log2FC")[, 1]), paste0(rr$seqnames, ":", rr$start - 1L, "-", rr$end)) }
LFC_RA <- get_lfc("markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds"); LFC_PV <- get_lfc("markers_GABA-MGE-PVALB-1_VS_GABA-MGE-LAMP5.rds")
g <- read_csv(file.path(HERE, paste0("scenicplus_all_genes_onoff_contrasts_", MAIN, ".csv.gz")), show_col_types = FALSE) %>% filter(!is.na(lfc_RA_vs_C1H1), !is.na(lfc_PV1_vs_low))
e <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", Region = "c", eRegulon_name = "c", importance_TF2G = "d"), progress = FALSE) %>%
  filter(grepl("\\+/\\+$", eRegulon_name)) %>% select(TF, Gene, Region, importance_TF2G) %>% filter(Gene %in% g$gene, Region %in% names(LFC_RA), Region %in% names(LFC_PV))

plot_one <- function(d, tf, kind, W = PANEL_IN + 1.75, H = PANEL_IN + 0.95) {
  d <- d %>% filter(is.finite(x), is.finite(y))
  r <- cor(d$x, d$y); rho <- cor(d$x, d$y, method = "spearman")
  lab <- d %>% filter(x > UP_LFC, y > UP_LFC) %>% mutate(w = pmin(x, y)) %>% arrange(desc(w)) %>% head(N_LABEL)
  big <- nrow(d) > 2000
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = imp), size = if (big) 0.5 else 1, alpha = if (big) 0.6 else 0.9, stroke = 0) +
    geom_text_repel(data = lab, aes(label = label), xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                    segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, point.padding = 0.05, force = 3, max.time = 5, max.overlaps = Inf, seed = 1) +
    annotate("text", x = min(d$x), y = max(d$y), label = paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(d), if (kind == "g") " genes" else " regions"),
             hjust = 0, vjust = 1, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
    scale_color_viridis_c(name = "log TF2G\nimportance") +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
    labs(title = paste0(tf, " regulon"), x = paste0("log2 fold change", if (kind == "r") " accessibility" else "", "\nRA vs C1H-1"),
         y = paste0("log2 fold change", if (kind == "r") " accessibility" else "", "\nPVALB-1 vs LAMP5")) +
    theme_fig() + theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
                        legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
                        plot.title = element_text(size = FIG_PT_AXIS_TITLE, family = FIG_FONT, hjust = 0))
  p <- egg::set_panel_size(p, width = unit(PANEL_IN, "in"), height = unit(PANEL_IN, "in"))
  stem <- file.path(D, if (kind == "g") "gene_expression" else "region_accessibility", tf)
  ggsave(paste0(stem, ".pdf"), p, width = W, height = H, device = cairo_pdf, family = FIG_FONT)
  ggsave(paste0(stem, ".svg"), p, width = W, height = H, bg = "transparent",
         device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
  tibble(TF = tf, kind = ifelse(kind == "g", "gene_expression", "region_accessibility"), n = nrow(d), pearson_r = r, spearman_rho = rho)
}

idx <- list()
for (tf in sort(unique(e$TF))) {
  et <- e %>% filter(TF == tf)
  if (n_distinct(et$Region) < MIN_REGIONS) next
  # genes: one row per target gene
  gi <- et %>% group_by(Gene) %>% summarise(imp = log(max(importance_TF2G)), .groups = "drop") %>% inner_join(g %>% select(gene, lfc_RA_vs_C1H1, lfc_PV1_vs_low), by = c("Gene" = "gene")) %>%
    transmute(label = nm(Gene), x = lfc_RA_vs_C1H1, y = lfc_PV1_vs_low, imp)
  idx[[length(idx) + 1]] <- plot_one(gi, tf, "g")
  # regions: one row per linked region; label = linked gene(s), color = highest importance among them
  ri <- et %>% group_by(Region) %>% summarise(label = paste(nm(unique(Gene)), collapse = "/"), imp = log(max(importance_TF2G)), .groups = "drop") %>%
    mutate(x = LFC_RA[Region], y = LFC_PV[Region])
  idx[[length(idx) + 1]] <- plot_one(ri, tf, "r")
}
IDX <- bind_rows(idx)
write_csv(IDX %>% mutate(across(where(is.double), ~ round(.x, 4))), file.path(D, "index.csv"))
cat("plots:", paste(names(table(IDX$kind)), table(IDX$kind), collapse = "; "), "->", D, "\n")
