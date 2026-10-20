#!/usr/bin/env Rscript
# Figure 4, panel F -- expression of predicted MAFB target genes in two MAFB on/off contrasts.
#
# x: log2 fold change RA vs C1H-1 (glutamatergic), y: log2 fold change PVALB-1 vs LGE-1 (GABAergic); pseudobulk DESeq2 on raw counts, one sample per library x
# bird x cell type with >= 10 cells, design ~ library + group when estimable. One point per MAFB +/+ target gene (config37) tested in both contrasts, keeping
# genes with |log2FC| > MIN_ABS in at least one contrast. Color: log TF-to-gene importance of MAFB -> gene (viridis). Labels: of the genes with log2FC >
# UP_LFC in both contrasts (a directional rule, no padj), the N_LABEL with the largest smaller-of-the-two log2FC, plus KCNC1, ERBB4 and PVALB and genes with
# literature support for fast-spiking / PV-class roles (grn/analysis/pv_fs_literature_genes.csv). Annotated with Pearson r, Spearman rho and n.
# Same computation as grn/analysis/scenicplus_pv1_vs_lge1_contrast.R + scenicplus_mafb_pvalb_vs_mafb_low_mge.R (RA vs C1H-1) +
# scenicplus_mafb_pvalb1_vs_lge1_scatters.R, which wrote
#   <scenicplus output>/config37/PVALB1_vs_LGE1/mafb_onoff_scatter_PVALB1_vs_LGE1_absLFC1.pdf
#
# Inputs: the hybrid-labeled Seurat object and the config37 SCENIC+ MuData (common.R).
# Output: FIG4_OUT/panel_F_mafb_onoff_scatter_PVALB1_vs_LGE1_absLFC1.{pdf,png,svg} and the plotted values (csv)
#
#   Rscript figures/figure4/panel_F_mafb_target_scatter.R
source(here::here("figures/figure4/common.R"))
suppressMessages({library(ggrepel); library(egg)})
PANEL_IN <- 1.6; UP_LFC <- 1; N_LABEL <- 10; MIN_ABS <- 1
LABEL_KEEP <- c("KCNC1", "ERBB4", "PVALB")
LIT_GENES <- c("KCNC1", "KCNC2", "PVALB", "ERBB4", "MAF", "MAFB", "SOX6", "CNTNAP2")
fig_check_font()

rna <- mafb_onoff_rna_contrasts()
imp <- read_eregulons() %>% filter(TF == "MAFB") %>% group_by(gene = Gene) %>% summarise(importance_TF2G = max(as.numeric(importance_TF2G)), .groups = "drop") %>%
  mutate(importance_log_TF2G = log(importance_TF2G))
d <- imp %>% inner_join(rna, by = "gene") %>% filter(!is.na(lfc_RA), !is.na(lfc_PV1)) %>% transmute(gene, x = lfc_RA, y = lfc_PV1, importance_log_TF2G)
cat("MAFB +/+ targets:", nrow(imp), "| tested in both contrasts:", nrow(d), "\n")
d <- d %>% filter(abs(x) > MIN_ABS | abs(y) > MIN_ABS)
r <- cor(d$x, d$y); rho <- cor(d$x, d$y, method = "spearman")
cand <- d %>% filter(x > UP_LFC, y > UP_LFC) %>% mutate(weak = pmin(x, y)) %>% arrange(desc(weak))
labeled <- union(union(head(cand$gene, N_LABEL), intersect(LABEL_KEEP, d$gene)), intersect(LIT_GENES, d$gene))
d <- d %>% mutate(label = gene_label(gene))
cat(sprintf("plotted: %d genes, Pearson r = %.3f, Spearman rho = %.3f\nlabeled: %s\n", nrow(d), r, rho, paste(gene_label(labeled), collapse = ", ")))
write_csv(d %>% mutate(across(where(is.double), ~ round(.x, 3))), file.path(FIG4_OUT, "panel_F_mafb_targets_log2fc.csv"))
p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(aes(color = importance_log_TF2G), size = 1, alpha = 0.9, stroke = 0) +
  geom_text_repel(data = d %>% filter(gene %in% labeled), aes(label = label), xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT,
                  color = FIG_INK_PRIMARY, segment.size = 0.2, segment.color = FIG_INK_MUTED, min.segment.length = 0.1, box.padding = 0.25, point.padding = 0.05, force = 3,
                  max.time = 15, max.iter = 300000, max.overlaps = Inf, seed = 1) +
  annotate("text", x = min(d$x), y = max(d$y), label = paste0("r = ", sprintf("%.2f", r), "\nρ = ", sprintf("%.2f", rho), "\nn = ", nrow(d), " genes"), hjust = 0, vjust = 1,
           size = fig_pt(FIG_PT_AXIS_TEXT), family = FIG_FONT, lineheight = 0.95) +
  scale_color_viridis_c(name = "log TF2G\nimportance") + scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) + coord_cartesian(clip = "off") +
  labs(x = "log2 fold change\nRA vs C1H-1", y = "log2 fold change\nPVALB-1 vs LGE-1") + theme_fig() +
  theme(plot.background = element_blank(), panel.background = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.25, "cm"), legend.title = element_text(size = FIG_PT_AXIS_TEXT), legend.text = element_text(size = FIG_PT_AXIS_TEXT))
p <- egg::set_panel_size(p, width = unit(PANEL_IN, "in"), height = unit(PANEL_IN, "in"))
save_gg(p, "panel_F_mafb_onoff_scatter_PVALB1_vs_LGE1_absLFC1", width = PANEL_IN + 1.75, height = PANEL_IN + 0.85)
