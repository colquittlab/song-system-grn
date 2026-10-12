#!/usr/bin/env Rscript
# Scatter of REGION differential accessibility, in the style of the DEG scatter (scenicplus_mafb_onoff_scatter.R): one point per region behind the
# 27 genes that are up in both RA and PVALB-1, x = log2FC accessibility RA vs C1H-1, y = log2FC accessibility PVALB-1 vs LAMP5 (ArchR, all regions,
# whatever their significance; table from scenicplus_mafb_region_differential_accessibility_heatmap.R). Color is the log TF2G importance of the
# linked gene (MAFB -> gene, config37; the highest if a region links to several genes), viridis. Labels are gene names: the 10 regions with the
# largest smaller-of-the-two log2FC (both must exceed UP_LFC), plus the best region of KCNC1, ERBB4 and PVALB so each is labeled once. LOC genes
# are named by protein homology or ncRNA-N (loc_gene_labels_lonStrDom2.csv).
#
#   Rscript scenicplus_mafb_region_dar_scatter.R
suppressMessages({library(tidyverse); library(ggrepel); library(egg); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
PANEL_IN <- 1.6; UP_LFC <- 1; N_LABEL <- 10
LABEL_KEEP <- c("KCNC1", "ERBB4", "PVALB")
LABEL_MAP <- c(LOC110468166 = "NTM", LOC110470685 = "TLE4", LOC110474504 = "VIPR1", LOC116183613 = "ncRNA-1")

x <- read_csv(file.path(HERE, "scenicplus_mafb_up_both_regions_differential_accessibility_config37.csv"), show_col_types = FALSE)
imp <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c", importance_TF2G = "d"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name)) %>% group_by(Gene) %>% summarise(imp = max(importance_TF2G), .groups = "drop")
gene_imp <- setNames(log(imp$imp), imp$Gene)
nm <- function(g) ifelse(g %in% names(LABEL_MAP), LABEL_MAP[g], g)
x <- x %>% mutate(imp = sapply(strsplit(genes, ";"), function(g) max(gene_imp[g], na.rm = TRUE)),
                  label = sapply(strsplit(genes, ";"), function(g) paste(nm(g), collapse = "/")),
                  weak = pmin(lfc_RA, lfc_PV1))
r <- cor(x$lfc_RA, x$lfc_PV1); rho <- cor(x$lfc_RA, x$lfc_PV1, method = "spearman")
cat(sprintf("regions: %d | Pearson r = %.3f, Spearman rho = %.3f\n", nrow(x), r, rho))

top <- x %>% filter(lfc_RA > UP_LFC, lfc_PV1 > UP_LFC) %>% arrange(desc(weak)) %>% head(N_LABEL)
keep <- x %>% filter(sapply(strsplit(genes, ";"), function(g) any(g %in% LABEL_KEEP))) %>%
  mutate(fg = sapply(strsplit(genes, ";"), function(g) LABEL_KEEP[LABEL_KEEP %in% g][1])) %>% group_by(fg) %>% slice_max(weak, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(label = fg)
lab <- bind_rows(top, keep) %>% distinct(region, .keep_all = TRUE)
cat("labeled regions:", paste(lab$label, collapse = ", "), "\n")

p <- ggplot(x, aes(lfc_RA, lfc_PV1)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(aes(color = imp), size = 1, alpha = 0.9, stroke = 0) +
  geom_text_repel(data = lab, aes(label = label), xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                  segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, point.padding = 0.05, force = 3, max.time = 15,
                  max.iter = 300000, max.overlaps = Inf, seed = 1) +
  annotate("text", x = min(x$lfc_RA), y = max(x$lfc_PV1), label = paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(x), " regions"),
           hjust = 0, vjust = 1, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
  scale_color_viridis_c(name = "log TF2G\nimportance") +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
  labs(x = "log2FC accessibility\nRA vs C1H-1", y = "log2FC accessibility\nPVALB-1 vs LAMP5") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
p <- egg::set_panel_size(p, width = unit(PANEL_IN, "in"), height = unit(PANEL_IN, "in"))
W <- PANEL_IN + 1.75; H <- PANEL_IN + 0.85
stem <- file.path(OUT, MAIN, "mafb_up_both_regions_dar_scatter")
fig_save(p, stem, width = W, height = H)
ggsave(paste0(stem, ".png"), p, width = W, height = H, dpi = 600, bg = "transparent")
ggsave(paste0(stem, ".svg"), p, width = W, height = H, bg = "transparent",
       device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
cat("saved", stem, "\n")
