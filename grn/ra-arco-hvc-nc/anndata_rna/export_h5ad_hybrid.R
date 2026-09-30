## Export the hybrid-labeled multiome object to the h5ad that anndata.ipynb reads.
##
## Same contract as the older convert_to_h5ad_adult-multiome.R (X = SCT `data`, cells x genes,
## the six obs columns downstream code expects), except that `cluster` is now `cluster_hybrid`
## from multiome/seurat/reduction_viz/combined_all_umap_hybrid.R. The object is read-only here;
## the h5ad is a separate copy.
##
## Usage: Rscript export_h5ad_hybrid.R [in_qs2] [out_h5ad]

suppressMessages({
  library(qs2)
  library(Seurat)
  library(reticulate)
})

args = commandArgs(trailingOnly = TRUE)
in_fname = if (length(args) >= 1) args[1] else
  "/ssd/brad/rstudio/multiome/song-system-grn/multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
out_fname = if (length(args) >= 2) args[2] else
  "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/data/adata.h5ad"
dir.create(dirname(out_fname), recursive = TRUE, showWarnings = FALSE)

use_python("/home/brad/micromamba/envs/scenicplus/bin/python", required = TRUE)
ad = import("anndata")

obj = qs_read(in_fname, nthreads = 8)
DefaultAssay(obj) = "SCT"
if (length(Layers(obj, assay = "SCT", search = "data")) > 1) obj = JoinLayers(obj)

md = obj@meta.data
stopifnot("cluster_hybrid has NA cells" = !anyNA(md$cluster_hybrid))
## The old script renamed sample_name -> sample_short; the object now carries both, so check they agree.
stopifnot("sample_name and sample_short disagree" = all(md$sample_name == md$sample_short))

obs = data.frame(
  cluster      = as.character(md$cluster_hybrid),
  position     = as.character(md$position),
  assignment   = as.character(md$assignment),
  sample_short = as.character(md$sample_short),
  run          = as.character(md$run),
  replicate    = as.character(md$replicate),
  row.names    = rownames(md),
  stringsAsFactors = FALSE
)

X = t(LayerData(obj, assay = "SCT", layer = "data"))  # cells x genes, sparse
stopifnot(identical(rownames(X), rownames(obs)))

## Both the hybrid-script UMAP (dims30...) and umap_rna_int are kept: anndata.ipynb aliases one of
## them to `archr_umap` for display only.
umap_names = c("dims30nn30mindist0.3", "umap_rna_int")
obsm = lapply(setNames(umap_names, paste0("X_", umap_names)), function(nm) {
  e = Embeddings(obj, reduction = nm)
  stopifnot(identical(rownames(e), rownames(obs)))
  unname(e)
})

adata = ad$AnnData(
  X = r_to_py(X),
  obs = r_to_py(obs),
  var = r_to_py(data.frame(name = colnames(X), row.names = colnames(X))),
  obsm = obsm
)
adata$write_h5ad(out_fname)

chk = ad$read_h5ad(out_fname)
print(chk)
print(table(chk$obs$cluster))
