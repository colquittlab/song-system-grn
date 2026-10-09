#!/usr/bin/env Rscript
# Scatter of MAFB +/+ target log2 fold changes in two MAFB on/off contrasts, in the style of the notebook's old RA vs C1H-1 / GABA-4 vs
# GABA-1 scatter (plot_deg_scatter in scenicplus_hybrid_all_template.qmd), with the contrasts from scenicplus_mafb_pvalb_vs_mafb_low_mge.R:
#   x = RA vs C1H-1     y = PVALB-1 (or PVALB-2) vs the MAFB-low MGE type (LAMP5)
# One point per target gene (the old plot used one row per region-gene link, which weights genes by their number of links). Color is
# log TF2G importance (log of the TF-to-gene importance of MAFB -> gene in config37, the notebook's importance_log_TF2G), single-hue
# sequential as the project's figure standard requires instead of viridis. Larger points are targets up in BOTH contrasts
# (padj < 0.05 and log2FC > 0.25 each). Annotated with Pearson r and Spearman rho over the plotted genes.
#
#   Rscript scenicplus_mafb_onoff_scatter.R
suppressMessages({library(tidyverse); library(ggrepel); library(here)})
source(here::here("config/figure_theme.R"))
fig_check_font()
MAIN <- "config37"
OUT <- path.expand("~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all/")
HERE <- here::here("grn/analysis")
LABEL <- c("PVALB", "KCNC1", "MAFB", "MAF", "RUNX2", "KIAA1217", "ERBB4", "GRIP1", "NXPH1", "RSPO2", "CNTNAP2", "COL19A1", "PCSK5", "CHN2", "ELFN1", "SGCZ")

W <- read_csv(file.path(HERE, paste0("scenicplus_mafb_pvalb_vs_mafb_low_mge_", MAIN, ".csv")), show_col_types = FALSE)
imp <- read_tsv(file.path(OUT, MAIN, "scenicplus_eRegulons.txt"), col_types = cols_only(TF = "c", Gene = "c", eRegulon_name = "c", importance_TF2G = "d"), progress = FALSE) %>%
  filter(TF == "MAFB", grepl("\\+/\\+$", eRegulon_name)) %>%
  group_by(gene = Gene) %>% summarise(importance_TF2G = max(importance_TF2G), .groups = "drop") %>%
  mutate(importance_log_TF2G = log(importance_TF2G))
D <- W %>% inner_join(imp, by = "gene")
cat("targets plotted:", nrow(D), "of", nrow(W), "testable\n")

make <- function(ycol, ylab, tag, min_abs = NULL) {
  d <- D %>% transmute(gene, x = lfc_RA_vs_C1H1, y = .data[[paste0("lfc_", ycol)]], importance_log_TF2G,
                       both = padj_RA_vs_C1H1 < 0.05 & lfc_RA_vs_C1H1 > 0.25 & .data[[paste0("padj_", ycol)]] < 0.05 & .data[[paste0("lfc_", ycol)]] > 0.25) %>%
    filter(!is.na(x), !is.na(y))
  n_all <- nrow(d)
  if (!is.null(min_abs)) { d <- d %>% filter(abs(x) > min_abs | abs(y) > min_abs); tag <- paste0(tag, "_absLFC", min_abs) }   # drop genes with no |log2FC| > min_abs in either contrast
  r <- cor(d$x, d$y); rho <- cor(d$x, d$y, method = "spearman")
  lab <- paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(d), " genes")
  cat(sprintf("%s: n=%d, Pearson r=%.3f, Spearman rho=%.3f, up in both=%d\n", tag, nrow(d), r, rho, sum(d$both)))
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = importance_log_TF2G, size = both), alpha = 0.9, stroke = 0) +
    geom_text_repel(data = d %>% filter(gene %in% LABEL), aes(label = gene), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, color = FIG_INK_PRIMARY,
                    segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, max.overlaps = Inf, seed = 1) +
    annotate("text", x = min(d$x), y = max(d$y), label = lab, hjust = 0, vjust = 1, size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
    scale_color_gradient(low = "#c6d8f0", high = "#0f3c78", name = "log TF2G\nimportance") +
    scale_size_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.6), labels = c(`FALSE` = "other", `TRUE` = "up in both"), name = NULL) +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
    labs(x = "log2 fold change, RA vs C1H-1", y = ylab) + theme_fig() +
    theme(aspect.ratio = 1, legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
  fig_save(p, file.path(OUT, MAIN, paste0("mafb_onoff_scatter_", tag)), width = 3.6, height = 2.8)
}
make("PV1_vs_low", "log2 fold change, PVALB-1 vs MAFB-low MGE (LAMP5)", "PVALB1_vs_LAMP5")
make("PV2_vs_low", "log2 fold change, PVALB-2 vs MAFB-low MGE (LAMP5)", "PVALB2_vs_LAMP5")
# same plots keeping only genes with |log2FC| > 0.5 in at least one of the two contrasts
make("PV1_vs_low", "log2 fold change, PVALB-1 vs MAFB-low MGE (LAMP5)", "PVALB1_vs_LAMP5", min_abs = 0.5)
make("PV2_vs_low", "log2 fold change, PVALB-2 vs MAFB-low MGE (LAMP5)", "PVALB2_vs_LAMP5", min_abs = 0.5)
# and at |log2FC| > 1
make("PV1_vs_low", "log2 fold change, PVALB-1 vs MAFB-low MGE (LAMP5)", "PVALB1_vs_LAMP5", min_abs = 1)
make("PV2_vs_low", "log2 fold change, PVALB-2 vs MAFB-low MGE (LAMP5)", "PVALB2_vs_LAMP5", min_abs = 1)
