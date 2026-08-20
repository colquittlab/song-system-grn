#library(SeuratDisk)
#library(SeuratWrappers)
library(Seurat)
#library(zellkonverter)
library(purrr)
library(qs2)
library(scCustomize)
library(reticulate)
library(tidyverse)
#library(sceasy)
#library(Rmagic)

#use_condaenv("/home/brad/micromamba/envs/SAMap")

in_dir = "/ssd/brad/rstudio/multiome/dev_full1/lonStrDom3/seurat/dev-full1_cellbender.0.1_multiome_seurat_preprocess/"
infile_obj = file.path(in_dir, "obj_clustered.qs2")

out_dir = file.path(in_dir)
h5ad_outfile = file.path(out_dir, "obj_clustered.h5ad")

obj = qs_read(infile_obj)


assay_to_export = "RNA"

DefaultAssay(obj) = assay_to_export
if (any(grepl("\\.[0-9]+", Layers(obj))))
  obj = JoinLayers(obj)

dimreducs = c("umap_weighted_rna40_atac30_mindist0.3nn30")
obj_slim = DietSeurat(
  obj,
  layers = c("counts", "data"),
  assays = assay_to_export,
  dimreducs = dimreducs
)


md = obj_slim@meta.data %>%
  rownames_to_column() %>%
  select(rowname, cluster, stage, ID, sex) %>%
  column_to_rownames()

obj_slim@meta.data = md


# Write h5ad directly with anndata — avoids sceasy/zellkonverter Seurat v5 incompatibilities
# and guarantees cellxgene-compatible obsm format
use_condaenv("/home/brad/micromamba/envs/magic")
ad <- reticulate::import("anndata")
np <- reticulate::import("numpy")

X_mat <- np$array(
  t(as.matrix(LayerData(obj_slim, assay = assay_to_export, layer = "counts"))),
  dtype = "float32"
)
umap_mat <- np$array(
  Embeddings(obj_slim, reduction = dimreducs),
  dtype = "float32"
)

adata <- ad$AnnData(
  X = X_mat,
  obs = reticulate::r_to_py(obj_slim@meta.data),
  var = reticulate::r_to_py(data.frame(row.names = rownames(obj_slim))),
  obsm = list(X_umap_seurat = umap_mat)
)
adata$write_h5ad(h5ad_outfile)


obj_h5ad = ad$read_h5ad(h5ad_outfile)
