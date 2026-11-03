#!/usr/bin/env Rscript
# 02_process_oe.R -- chicken MAFB/EMX2 overexpression snRNA (PIPseq run 260312): cells, embryos, cell classes, plasmid.
#
# Libraries (one dissection each, several embryos pooled; embryos demultiplexed by souporcell):
#   MAFBn  dox-inducible MAFB, -dox   (souporcell k = 2)
#   MAFBp  dox-inducible MAFB, +dox   (k = 3)
#   EMX2p  dox-inducible EMX2, +dox   (k = 4)   -- no -dox EMX2 library
# Plasmid (pPIDNB2): mNeonGreen is constitutive (marks electroporated cells); the TF (MAFB-tg / EMX2-tg) and mScarlet are
# dox-inducible. The three libraries share no embryos: depth-filtered souporcell genotypes are equally discordant (0.32-0.45)
# for every cross-library cluster pair.
#
# Cell classes are assigned per cluster from canonical markers and printed for checking. Plasmid status per cell:
#   plasmid  = any plasmid read (mNeonGreen, mScarlet or the TF transgene) -- the cell was electroporated
#   inducible = any MAFB-tg / EMX2-tg / mScarlet read
# The original analysis (chick_mafb-emx2_run1.qmd) used MAFB-tg >= 1 / EMX2-tg >= 1 only, pooled across libraries, so -dox
# cells with leaky MAFB-tg were counted as MAFB-overexpressing.
#
#   Rscript overexpression/reanalysis/02_process_oe.R
# Writes output/oe_obj.qs2, output/oe_cells.csv.gz, tables/oe_cell_counts.csv, tables/oe_plasmid_by_library.csv.
suppressMessages({library(Seurat); library(qs2); library(rtracklayer)})
source(here::here("overexpression/reanalysis/oe_helpers.R"))
set.seed(10)

LIBS <- tribble(
  ~lib,    ~subdir,                 ~k, ~tf,    ~dox,
  "MAFBn", "bGalGal1.mat.mafb_v2",  2,  "MAFB", "minus",
  "MAFBp", "bGalGal1.mat.mafb_v2",  3,  "MAFB", "plus",
  "EMX2p", "bGalGal1.mat.emx2",     4,  "EMX2", "plus"
)
PLASMID_GENES <- c("MAFB-tg", "EMX2-tg", "mScarlet", "mNeonGreen")
QC <- list(max_mito = 2, min_count = 500, max_count = 20000, min_feature = 300)

## 1. counts, padded to the union of features (each library's reference only carries its own TF transgene)
mats <- pmap(LIBS, function(lib, subdir, ...) {
  f <- list.files(file.path(PIPSEQ_MAFB_EMX2_DIR, subdir, lib), "scRNA\\.filtered", full.names = TRUE)
  m <- ReadMtx(mtx = grep("matrix", f, value = TRUE), cells = grep("barcodes", f, value = TRUE), features = grep("features", f, value = TRUE))
  colnames(m) <- paste(lib, colnames(m), sep = "_")
  m
}) %>% set_names(LIBS$lib)
all_genes <- Reduce(union, map(mats, rownames))
mats <- map(mats, function(m) {
  miss <- setdiff(all_genes, rownames(m))
  if (length(miss)) m <- rbind(m, Matrix(0, length(miss), ncol(m), sparse = TRUE, dimnames = list(miss, colnames(m))))
  m[all_genes, ]
})
cnt <- do.call(cbind, mats)

## 2. souporcell: singlets only; embryo = library x souporcell cluster
sc <- pmap_dfr(LIBS, function(lib, subdir, k, ...) {
  read_tsv(file.path(PIPSEQ_MAFB_EMX2_DIR, subdir, lib, sprintf("souporcell_%d", k), "clusters.tsv"), show_col_types = FALSE) %>%
    transmute(cell = paste(lib, barcode, sep = "_"), status, assignment)
})
cells <- tibble(cell = colnames(cnt), lib = sub("_.*", "", colnames(cnt))) %>% left_join(sc, by = "cell") %>% left_join(LIBS %>% select(lib, tf, dox), by = "lib")

## 3. QC
annot <- import(GTF_CHICKEN_MAFB_EMX2)
is_mt <- as.character(seqnames(annot)) == "chrMT"
mt <- intersect(unique(c(annot$gene_name[is_mt], annot$gene[is_mt], annot$gene_id[is_mt], str_replace(annot$gene_id[is_mt], "_", "-"))), rownames(cnt))
cells <- cells %>% mutate(nCount = colSums(cnt), nFeature = colSums(cnt > 0), perc_mito = 100 * colSums(cnt[mt, , drop = FALSE]) / nCount,
                          pass_qc = perc_mito < QC$max_mito & nCount >= QC$min_count & nCount <= QC$max_count & nFeature >= QC$min_feature)
cat(sprintf("mitochondrial genes: %d\n", length(mt)))
print(cells %>% count(lib, status, pass_qc) %>% pivot_wider(names_from = pass_qc, values_from = n, names_prefix = "qc_"))
cells <- cells %>% filter(status == "singlet", pass_qc) %>% mutate(embryo = paste(lib, assignment, sep = "."))

## 4. plasmid status (before any normalization)
pl <- as.matrix(t(cnt[PLASMID_GENES, cells$cell]))
cells <- cells %>% mutate(mNG = pl[, "mNeonGreen"], mScarlet = pl[, "mScarlet"], tf_tg = ifelse(tf == "MAFB", pl[, "MAFB-tg"], pl[, "EMX2-tg"]),
                          plasmid = (mNG + mScarlet + tf_tg) > 0, inducible = (mScarlet + tf_tg) > 0)

## 5. Seurat: SCT, PCA, Harmony over library (clusters are for cell classes only; DE uses raw counts), clusters
obj <- CreateSeuratObject(cnt[, cells$cell], meta.data = cells %>% column_to_rownames("cell") %>% as.data.frame())
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- harmony::RunHarmony(obj, group.by.vars = "lib", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)

## 6. cell classes from markers: per cluster, the class whose markers have the highest mean detection fraction, centered
## per gene across clusters. (Scaling instead of centering inflates sparsely detected markers: a z-scored version called a
## CSF1R-negative cluster "Micro".)
MARKERS <- list(
  Glut = c("SLC17A6", "SLC17A7", "NEUROD6", "TBR1", "SATB2"),
  GABA = c("GAD1", "GAD2", "SLC32A1", "DLX5", "DLX6"),
  Progenitor = c("SOX2", "PAX6", "HES5", "FABP7", "MKI67", "TOP2A", "VIM"),
  Astro = c("AQP4", "SLC1A3", "GJA1", "ALDH1L1"),
  OPC_Oligo = c("PDGFRA", "OLIG2", "SOX10", "PLP1", "MBP"),
  Micro = c("CSF1R", "C1QB", "PTPRC", "CX3CR1"),
  Endo_Mural = c("CLDN5", "FLT1", "PECAM1", "PDGFRB", "RGS5"),
  Epen = c("FOXJ1", "CFAP43", "DNAH9")
)
MARKERS <- map(MARKERS, ~ intersect(.x, rownames(obj)))
mk <- unique(unlist(MARKERS))
det_cl <- sapply(split(seq_len(ncol(obj)), obj$seurat_clusters), function(i) Matrix::rowMeans(GetAssayData(obj, assay = "RNA", layer = "counts")[mk, i, drop = FALSE] > 0))
z <- det_cl - rowMeans(det_cl)
cls_score <- sapply(MARKERS, function(g) colMeans(z[g, , drop = FALSE]))
rownames(cls_score) <- colnames(z)
cl_margin <- apply(cls_score, 1, function(x) diff(sort(x, decreasing = TRUE)[2:1]))
cl_class <- setNames(ifelse(cl_margin < 0.03, "Unresolved", colnames(cls_score)[apply(cls_score, 1, which.max)]), rownames(cls_score))
obj$cell_class <- unname(cl_class[as.character(obj$seurat_clusters)])
cat("\ncluster -> class (mean centered marker detection; margin < 0.03 -> Unresolved):\n")
print(as.data.frame(cls_score) %>% rownames_to_column("cluster") %>% mutate(class = cl_class[cluster], n = as.integer(table(obj$seurat_clusters)[cluster]),
                                                                      margin = apply(cls_score, 1, function(x) diff(sort(x, decreasing = TRUE)[2:1]))) %>%
        mutate(across(where(is.double), ~ round(.x, 2))), row.names = FALSE)

## 7. tables
md <- obj@meta.data %>% rownames_to_column("cell") %>% as_tibble()
counts_tab <- md %>% group_by(lib, dox, tf, embryo, cell_class) %>%
  summarise(cells = n(), plasmid = sum(plasmid), mNG = sum(mNG > 0), inducible = sum(inducible), tf_tg = sum(tf_tg > 0), .groups = "drop")
write_csv(counts_tab, file.path(OE_TABLES, "oe_cell_counts.csv"))
cat("\nglutamatergic cells per embryo (plasmid = any plasmid read):\n")
print(as.data.frame(counts_tab %>% filter(cell_class == "Glut")), row.names = FALSE)

# Is the TF transgene dox-induced? Detection depends on depth, so compare per 10k UMI as well as per cell.
lib_tab <- md %>% group_by(lib, dox, tf) %>%
  summarise(cells = n(), median_umi = median(nCount_RNA),
            pct_mNG = 100 * mean(mNG > 0), pct_tf_tg = 100 * mean(tf_tg > 0), pct_mScarlet = 100 * mean(mScarlet > 0),
            tf_tg_umi_per_10k = 1e4 * sum(tf_tg) / sum(nCount_RNA), mScarlet_umi_per_10k = 1e4 * sum(mScarlet) / sum(nCount_RNA),
            mNG_umi_per_10k = 1e4 * sum(mNG) / sum(nCount_RNA),
            tf_tg_umi_per_pos_cell = sum(tf_tg) / max(1, sum(tf_tg > 0)),
            inducible_among_mNG = sum(inducible & mNG > 0) / max(1, sum(mNG > 0)), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(lib_tab, file.path(OE_TABLES, "oe_plasmid_by_library.csv"))
cat("\nplasmid transcripts by library:\n"); print(as.data.frame(lib_tab), row.names = FALSE)
# Depth bias: cells with more UMIs are more likely to show any plasmid read
cat("\nmedian UMI, plasmid+ vs plasmid- glutamatergic cells, by library:\n")
print(as.data.frame(md %>% filter(cell_class == "Glut") %>% group_by(lib, plasmid) %>% summarise(n = n(), median_umi = median(nCount_RNA), .groups = "drop")), row.names = FALSE)

qs_save(obj, file.path(OE_OUT, "oe_obj.qs2"), nthreads = 8)
write_csv(md %>% select(cell, lib, dox, tf, embryo, assignment, cell_class, seurat_clusters, nCount_RNA, nFeature_RNA, perc_mito, mNG, mScarlet, tf_tg, plasmid, inducible),
          file.path(OE_OUT, "oe_cells.csv.gz"))

## 8. QC figures: UMAP by class and by plasmid status
emb <- as_tibble(Embeddings(obj, "umap"), rownames = "cell") %>% setNames(c("cell", "u1", "u2")) %>% left_join(md, by = "cell")
class_lv <- names(sort(table(emb$cell_class), decreasing = TRUE))
p1 <- ggplot(emb, aes(u1, u2, color = factor(cell_class, levels = class_lv))) + geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(values = setNames(rep_len(FIG_PAL, length(class_lv)), class_lv), name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) + coord_equal() + labs(x = "UMAP 1", y = "UMAP 2") + theme_fig() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), legend.text = element_text(size = FIG_PT_AXIS_TEXT), legend.key.size = unit(6, "pt"))
fig_save(p1, file.path(OE_OUT, "umap_cell_class"), width = 3.2, height = 2.4)
p2 <- ggplot(emb %>% arrange(plasmid), aes(u1, u2, color = plasmid)) + geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(values = c(`FALSE` = FIG_NEUTRAL_FILL, `TRUE` = FIG_PAL[2]), labels = c("no plasmid read", "plasmid read"), name = NULL) +
  facet_wrap(~ lib) + guides(color = guide_legend(override.aes = list(size = 1.5))) + coord_equal() + labs(x = "UMAP 1", y = "UMAP 2") + theme_fig() +
  theme(axis.text = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        strip.text = element_text(size = FIG_PT_AXIS_TITLE))
fig_save(p2, file.path(OE_OUT, "umap_plasmid_by_library"), width = 4.5, height = 2)
cat("\ndone\n")
