## Xenium spatial-localisation score aligned to the snRNA-seq positional
## specificity figure from snrna/trees/celltypes_hclust_all_hybrid.qmd.
##
## Rebuilds that notebook's vertical panel (cell types x dissection position,
## doubly-normalised, with the tissue-specificity bar) using the obj_subset.qs2
## it already saved, then appends the Xenium per-type score in the SAME row
## order so the two modalities can be read across.
##
## snRNA side  : tissue specificity of position proportions (calc_tissue_specificity)
## Xenium side : spatial localisation delta from spatial_variation_by_type.R
##               (JSD of a type's bin distribution from all cells, minus a
##               matched-n permutation null)
##
## Glut-DACH2-HVCra-Int is not scored on the Xenium side: after the
## confident + nCount>=400 filter it has only 24-74 cells per section, below
## the n=100 needed for a comparable index. It is drawn as an empty slot rather
## than given a value computed at a different n.

suppressMessages({
  library(Seurat); library(BPCells); library(tidyverse); library(qs2)
  library(ComplexHeatmap); library(circlize); library(here)
})
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "stats.R"))
select = dplyr::select

out_dir = here::here("xenium/label_transfer", "spatial_variation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj_fname = here::here("snrna/trees/celltypes_hclust_all_hybrid", "obj_subset.qs2")
res_to_use = "celltype_hybrid"

ct_order <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Int", "Glut-DACH2-HVCx",
             "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
             "Glut-DACH2-LMANsh", "Glut-DACH2-LMANco",
             "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8", "Glut-SATB2-1",
             "Glut-CACNA1H-RA", "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
             "Glut-GABA", "Glut-Im", "Glut-NB", "Glut-NSC",
             "GABA-LGE-1", "GABA-LGE-2",
             "GABA-MGE-SST-1", "GABA-MGE-SST-2", "GABA-MGE-SST-3",
             "GABA-MGE-PVALB-1", "GABA-MGE-PVALB-2", "GABA-MGE-PVALB-3",
             "GABA-MGE-LAMP5", "GABA-MGE-ST18", "GABA-MGE-LHX8",
             "GABA-CGE-1", "GABA-CGE-2", "GABA-Im",
             "Astro-1", "Astro-2", "Astro-3", "Epen",
             "Oligo-1", "Oligo-2", "Oligo-3", "OPC", "Micro", "Endo")

## ---- snRNA side: position x celltype, exactly as the tree notebook ---------
obj = qs_read(obj_fname, nthreads = 8)
stopifnot(res_to_use %in% colnames(obj@meta.data))
md = obj@meta.data %>% select(position, all_of(res_to_use)); rm(obj); gc()
stopifnot(setequal(ct_order, unique(md[[res_to_use]])))

tab = table(md$position, md[[res_to_use]])
tab = tab[c("hvc","nc","lman","nr","ra","arco"), ct_order]
tab = sweep(tab, 1, rowSums(tab), FUN = "/")
tab = sweep(tab, 2, colSums(tab), FUN = "/")
tab_spec = apply(tab, 2, calc_tissue_specificity)[ct_order]

## ---- Xenium side ----------------------------------------------------------
per_type = read_csv(file.path(out_dir,"spatial_localisation_by_type.csv"), show_col_types=FALSE)
per_sec  = read_csv(file.path(out_dir,"spatial_localisation_by_type_section.csv"), show_col_types=FALSE)
on = per_sec %>% group_by(first_type) %>%
  summarize(obs=mean(obs), null=mean(null), .groups="drop")
x = tibble(cluster = ct_order) %>%
  left_join(per_type %>% select(first_type, delta, delta_sd), by=c("cluster"="first_type")) %>%
  left_join(on, by=c("cluster"="first_type"))
scored = !is.na(x$delta)
cat("Xenium score available for ", sum(scored), "/", length(ct_order), " types\n", sep="")
cat("not scored: ", paste(x$cluster[!scored], collapse=", "), "\n", sep="")

## ---- shared colours (same scheme as the tree notebook) --------------------
glut_precursor = c("Glut-NSC","Glut-NB","Glut-Im"); gaba_precursor = c("GABA-Im")
bar_cols = data.frame(cluster_type=c("Glut","Glut-Pre","GABA","GABA-Pre","Non-neuron"),
                      col=c("#008000","#90ca90","#800000","#c37474","#808000"))
ct_df = data.frame(cluster = ct_order) %>%
  mutate(cluster_type = case_when(cluster %in% glut_precursor ~ "Glut-Pre",
                                  grepl("^Glut", cluster) ~ "Glut",
                                  cluster %in% gaba_precursor ~ "GABA-Pre",
                                  grepl("^GABA", cluster) ~ "GABA",
                                  TRUE ~ "Non-neuron")) %>%
  left_join(bar_cols, by="cluster_type")
bar_colors = ct_df$col

## Unscored types stay NA so their slot draws EMPTY. Filling them with 0 would
## render as a real measurement of "not localised", which is the opposite of
## "not measurable at comparable n".
delta_draw = x$delta; obs_draw = x$obs; null_draw = x$null

cols = colorRamp2(breaks = seq(0, 1, length.out = 9),
                  colors = scales::brewer_pal(palette="Greys")(9))
RH = 0.13; hm_h = length(ct_order)*RH; hm_w = nrow(tab)*RH

ann_snrna = HeatmapAnnotation(
  which = "row",
  "snRNA\nspecificity" = anno_barplot(tab_spec, 0, border=TRUE, gp=gpar(fill=bar_colors),
                                      width=unit(1.1,"in"), axis_param=list(gp=gpar(fontsize=5))),
  annotation_name_gp = gpar(fontsize=6), gap = unit(2,"mm"))
ann_xen = HeatmapAnnotation(
  which = "row",
  "Xenium\nlocalisation" = anno_barplot(delta_draw, 0, border=TRUE, gp=gpar(fill=bar_colors),
                                        width=unit(1.1,"in"), axis_param=list(gp=gpar(fontsize=5))),
  annotation_name_gp = gpar(fontsize=6), gap = unit(2,"mm"))
## observed vs permutation null, same rows
ann_pts = HeatmapAnnotation(
  which = "row",
  "obs vs null" = anno_points(cbind(obs_draw, null_draw), border=TRUE,
                              pch=c(16,1), size=unit(1.1,"mm"),
                              gp=gpar(col=c("black","grey55")),
                              width=unit(1.1,"in"), axis_param=list(gp=gpar(fontsize=5))),
  annotation_name_gp = gpar(fontsize=6), gap = unit(2,"mm"))

hm = Heatmap(t(tab), cluster_rows=FALSE, cluster_columns=FALSE,
             show_row_names=TRUE, row_names_side="left", col=cols,
             border=TRUE, border_gp=gpar(lwd=0.5),
             height=unit(hm_h,"in"), width=unit(hm_w,"in"),
             column_names_gp=gpar(fontsize=6), row_names_gp=gpar(fontsize=6),
             heatmap_legend_param=list(title="prop", labels_gp=gpar(fontsize=5),
                                       title_gp=gpar(fontsize=6)))

for (variant in c("bars","bars_points")) {
  ht = if (variant=="bars") hm + ann_snrna + ann_xen else hm + ann_snrna + ann_xen + ann_pts
  w = if (variant=="bars") 7.5 else 9
  pdf(file.path(out_dir, paste0("spatial_variation_aligned_",variant,".pdf")),
      height=hm_h+2, width=w); draw(ht, merge_legend=TRUE); dev.off()
  png(file.path(out_dir, paste0("spatial_variation_aligned_",variant,".png")),
      height=(hm_h+2)*300, width=w*300, res=300); draw(ht, merge_legend=TRUE); dev.off()
}
cat("wrote aligned figures\n")

## ---- do the two modalities agree? ----------------------------------------
cmp = x %>% mutate(snrna_spec = as.numeric(tab_spec[cluster]),
                   cluster_type = ct_df$cluster_type) %>%
  filter(!is.na(delta))
ct_test = cor.test(cmp$snrna_spec, cmp$delta, method="spearman")
cat(sprintf("\nSpearman snRNA specificity vs Xenium localisation: rho=%.3f, p=%.3g, n=%d\n",
            ct_test$estimate, ct_test$p.value, nrow(cmp)))
write_csv(cmp %>% select(cluster, cluster_type, snrna_spec, xenium_delta=delta,
                         xenium_obs=obs, xenium_null=null),
          file.path(out_dir,"snrna_vs_xenium_specificity.csv"))

p = ggplot(cmp, aes(snrna_spec, delta, colour=cluster_type)) +
  geom_point(size=2) +
  ggrepel::geom_text_repel(aes(label=cluster), size=2, max.overlaps=12, show.legend=FALSE) +
  scale_colour_manual(values=setNames(bar_cols$col, bar_cols$cluster_type), name=NULL) +
  labs(x="snRNA-seq positional specificity", y="Xenium spatial localisation (delta)",
       subtitle=sprintf("Spearman rho = %.2f, p = %.2g, n = %d",
                        ct_test$estimate, ct_test$p.value, nrow(cmp))) +
  cowplot::theme_cowplot(font_size=10)
ggsave(file.path(out_dir,"snrna_vs_xenium_specificity.png"), p, width=7.5, height=6, dpi=300, bg="white")
cat("DONE\n")
