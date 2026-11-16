## One clustering of the MAFB-domain accessible elements on four features: whole-tree and oscine-lineage
## (Anc222) mean phyloP and block % (% of bases in 11-bp blocks with mean phyloP > 1), each over the
## lobe's full extent, accessibility-weighted (element_profiles.py). Each column is z-scored; lobes
## are clustered by Ward's method on Euclidean distance and drawn with their dendrogram, unsplit.
## A strip marks the k groups (default 3: the Gaussian-mixture BIC minimum for these features);
## rows are labelled with the lobe's scored extent (1-based) and class (RA, G4, both).
##
##   Rscript plot_element_cluster.R [k=3]
## Writes results/mafb_domain/fig_element_cluster_meanblock_k<k>.{pdf,png} and element_cluster_meanblock_k<k>.tsv.
suppressMessages({library(ComplexHeatmap); library(circlize); library(data.table)})
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(if (length(.self)) normalizePath(file.path(dirname(.self), "../../../config/figure_theme.R")) else
  here::here("config/figure_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args)) as.integer(args[1]) else 3L
res <- "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
S <- fread(file.path(res, "element_rule_classes.tsv"))
feat <- c(whole_ext_mean = "whole-tree mean", whole_ext_block = "whole-tree block %",
          osc_ext_mean = "oscine mean", osc_ext_block = "oscine block %")
X <- scale(as.matrix(S[, names(feat), with = FALSE])); colnames(X) <- feat
hc <- hclust(dist(X), method = "ward.D2")
cl <- cutree(hc, k = k)
## name groups from their centroids: high (z > 0) on whole-tree and/or oscine features
cz <- apply(X, 2, function(v) tapply(v, cl, mean))
nm <- ifelse(rowMeans(cz[, 1:2, drop = FALSE]) > 0 & rowMeans(cz[, 3:4, drop = FALSE]) > 0, "deep+oscine",
       ifelse(rowMeans(cz[, 1:2, drop = FALSE]) > 0, "deep", ifelse(rowMeans(cz[, 3:4, drop = FALSE]) > 0, "oscine", "weak")))
nm <- make.unique(nm, sep = " "); S[, group := nm[cl]]
print(data.table(group = nm, n = as.vector(table(cl)), round(cz, 2)))
fwrite(S[, .(element, cls, group, chrom, ext_start, ext_end, n_lobes, lobe, whole_ext_mean, whole_ext_block,
             osc_ext_mean, osc_ext_block)], file.path(res, sprintf("element_cluster_meanblock_k%d.tsv", k)), sep = "\t")

lbl <- sprintf("%s:%s-%s%s  %s", S$chrom, format(S$ext_start + 1, big.mark = ",", trim = TRUE),
               format(S$ext_end, big.mark = ",", trim = TRUE),
               ifelse(S$n_lobes > 1, sprintf(" (lobe %d/%d)", S$lobe, S$n_lobes), ""), S$cls)
gcol <- setNames(c("#c8161a", "#1baf7a", "#9a9a9a", "#7a3ca0", "#eda100", "#2a78d6")[seq_along(nm)], nm)
gcol[grepl("^deep\\+oscine", nm)] <- "#7a3ca0"; gcol[grepl("^deep( |$)", nm)] <- "#c8161a"
gcol[grepl("^oscine", nm)] <- "#1baf7a"; gcol[grepl("^weak", nm)] <- "#9a9a9a"
ra <- rowAnnotation(group = S$group, class = S$cls,
                    col = list(group = gcol, class = c(RA = FIG_PAL[[2]], G4 = FIG_PAL[[1]], both = FIG_PAL[[3]])),
                    annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "bottom", annotation_name_rot = 45, simple_anno_size = unit(2, "mm"),
                    annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
ht <- Heatmap(X, name = "z", col = colorRamp2(c(-2, 0, 2), c(FIG_DIV_LOW, FIG_DIV_MID, FIG_DIV_HIGH)),
              cluster_rows = hc, row_dend_width = unit(18, "mm"), cluster_columns = FALSE,
              right_annotation = ra, row_labels = lbl, row_names_gp = gpar(fontsize = 5, fontfamily = "mono"),
              column_names_gp = gpar(fontsize = 6), column_names_rot = 45, width = unit(4 * 6, "mm"),
              column_title = sprintf("Ward, z-scored (k = %d)", k), column_title_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)))
h <- 1.4 + nrow(S) * 2 / 25.4
for (dev in c("pdf", "png")) {
  f <- file.path(res, sprintf("fig_element_cluster_meanblock_k%d.%s", k, dev))
  if (dev == "pdf") cairo_pdf(f, width = 5.2, height = h, family = FIG_FONT) else png(f, width = 5.2, height = h, units = "in", res = 600, family = FIG_FONT)
  draw(ht, merge_legend = TRUE); dev.off()
}
message("wrote ", file.path(res, sprintf("fig_element_cluster_meanblock_k%d", k)), ".pdf/.png")
