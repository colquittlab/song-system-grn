#!/usr/bin/env Rscript
# Scatter of MAFB +/+ target log2 fold changes in two MAFB on/off contrasts, in the style of the notebook's old RA vs C1H-1 / GABA-4 vs
# GABA-1 scatter (plot_deg_scatter in scenicplus_hybrid_all_template.qmd), with the contrasts from scenicplus_mafb_pvalb_vs_mafb_low_mge.R:
#   x = RA vs C1H-1     y = PVALB-1 (or PVALB-2) vs the MAFB-low MGE type (LAMP5)
# One point per target gene (the old plot used one row per region-gene link, which weights genes by their number of links). Color is
# log TF2G importance (log of the TF-to-gene importance of MAFB -> gene in config37, the notebook's importance_log_TF2G), viridis
# (the project's standard calls for single-hue sequential; viridis was approved for this figure). Of the genes with log2FC > 1 in BOTH contrasts (a directional rule, no padj), the 10 with the largest smaller-of-the-two log2FC are labeled, plus KCNC1, ERBB4 and PVALB, plus genes with literature support for fast-spiking / PV-class roles. Annotated with Pearson r and Spearman rho over the plotted genes.
#
#   Rscript scenicplus_mafb_onoff_scatter.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
PANEL_IN <- 1.6   # side of the square plotting region, inches
UP_LFC <- 1   # directional rule, no padj: a gene is a label candidate if log2FC > UP_LFC in BOTH contrasts
N_LABEL <- 10 # ... and only the N candidates with the largest smaller-of-the-two log2FC are labeled
LABEL_KEEP <- c("KCNC1", "ERBB4", "PVALB")   # always labeled (when plotted), whether or not they are in the top N
# Genes with literature support for a role in fast-spiking / PV-class interneurons that are MAFB +/+ targets in config37 (evidence and sources:
# pv_fs_literature_genes.csv). Labeled when plotted (same style as the rank-based labels).
LIT_GENES <- c("KCNC1", "KCNC2", "PVALB", "ERBB4", "MAF", "MAFB", "SOX6", "CNTNAP2")
# LOC genes: human counterpart by protein similarity where there is one, otherwise ncRNA-N (basis and evidence: loc_gene_labels_lonStrDom2.csv)
LABEL_MAP <- c(LOC110468166 = "NTM", LOC110470685 = "TLE4", LOC110474504 = "VIPR1", LOC116183613 = "ncRNA-1")

W <- read_csv(file.path(HERE, paste0("scenicplus_mafb_pvalb_vs_mafb_low_mge_", MAIN, ".csv")), show_col_types = FALSE)
imp <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c", importance_TF2G = "d"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name)) %>%
  group_by(gene = Gene) %>% summarise(importance_TF2G = max(importance_TF2G), .groups = "drop") %>%
  mutate(importance_log_TF2G = log(importance_TF2G))
D <- W %>% inner_join(imp, by = "gene")
cat("targets plotted:", nrow(D), "of", nrow(W), "testable\n")

make <- function(ycol, ylab, tag, min_abs = NULL) {
  d <- D %>% transmute(gene, x = lfc_RA_vs_C1H1, y = .data[[paste0("lfc_", ycol)]], importance_log_TF2G,
                       both = lfc_RA_vs_C1H1 > UP_LFC & .data[[paste0("lfc_", ycol)]] > UP_LFC) %>%
    filter(!is.na(x), !is.na(y))
  n_all <- nrow(d)
  if (!is.null(min_abs)) { d <- d %>% filter(abs(x) > min_abs | abs(y) > min_abs); tag <- paste0(tag, "_absLFC", min_abs) }   # drop genes with no |log2FC| > min_abs in either contrast
  r <- cor(d$x, d$y); rho <- cor(d$x, d$y, method = "spearman")
  cand <- d %>% filter(both) %>% mutate(weak = pmin(x, y)) %>% arrange(desc(weak))
  labeled <- union(union(head(cand$gene, N_LABEL), intersect(LABEL_KEEP, d$gene)), intersect(LIT_GENES, d$gene))
  d <- d %>% mutate(label = ifelse(gene %in% names(LABEL_MAP), LABEL_MAP[gene], gene))
  cat("  labeled:", paste(ifelse(labeled %in% names(LABEL_MAP), paste0(LABEL_MAP[labeled], " (", labeled, ")"), labeled), collapse = ", "), "\n")
  lab <- paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(d), " genes")
  cat(sprintf("%s: n=%d, Pearson r=%.3f, Spearman rho=%.3f\n", tag, nrow(d), r, rho))
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = importance_log_TF2G), size = 1, alpha = 0.9, stroke = 0) +
    geom_text_repel(data = d %>% filter(gene %in% labeled), aes(label = label), xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                    segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, point.padding = 0.05, force = 3, max.time = 15, max.iter = 300000, max.overlaps = Inf, seed = 1) +
    annotate("text", x = min(d$x), y = max(d$y), label = lab, hjust = 0, vjust = 1, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
    scale_color_viridis_c(name = "log TF2G\nimportance") +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
    labs(x = "log2 fold change\nRA vs C1H-1", y = ylab) + theme_fig() +
    theme(legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
  p <- egg::set_panel_size(p, width = unit(PANEL_IN, "in"), height = unit(PANEL_IN, "in"))   # the plotting region itself is PANEL_IN x PANEL_IN inches
  W <- PANEL_IN + 1.75; H <- PANEL_IN + 0.85
  fig_save(p, file.path(OUT, MAIN, paste0("mafb_onoff_scatter_", tag)), width = W, height = H)
  # Editable version: cairo_pdf puts all the data labels into ONE text object (an editor imports them as a single block), so also write an
  # SVG where every label is its own <text> element, in Arial (a live font, not outlines or embedded).
  ggsave(file.path(OUT, MAIN, paste0("mafb_onoff_scatter_", tag, ".svg")), p, width = W, height = H, bg = FIG_SURFACE,
         device = function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, system_fonts = list(sans = FIG_FONT), ...))
}
make("PV1_vs_low", "log2 fold change\nPVALB-1 vs LAMP5", "PVALB1_vs_LAMP5")
make("PV2_vs_low", "log2 fold change\nPVALB-2 vs LAMP5", "PVALB2_vs_LAMP5")
# same plots keeping only genes with |log2FC| > 0.5 in at least one of the two contrasts
make("PV1_vs_low", "log2 fold change\nPVALB-1 vs LAMP5", "PVALB1_vs_LAMP5", min_abs = 0.5)
make("PV2_vs_low", "log2 fold change\nPVALB-2 vs LAMP5", "PVALB2_vs_LAMP5", min_abs = 0.5)
# and at |log2FC| > 1
make("PV1_vs_low", "log2 fold change\nPVALB-1 vs LAMP5", "PVALB1_vs_LAMP5", min_abs = 1)
make("PV2_vs_low", "log2 fold change\nPVALB-2 vs LAMP5", "PVALB2_vs_LAMP5", min_abs = 1)
