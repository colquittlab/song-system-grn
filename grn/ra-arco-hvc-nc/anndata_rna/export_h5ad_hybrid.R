## Export the hybrid-labeled multiome object to the h5ad that anndata.ipynb reads.
##
## Same obs contract as the older convert_to_h5ad_adult-multiome.R (cells x genes, the six obs columns
## downstream code expects), except that `cluster` is now `cluster_hybrid` from
## multiome/seurat/reduction_viz/combined_all_umap_hybrid.R, and X is raw counts rather than SCT data
## (see below). The object is read-only here; the h5ad is a separate copy.
##
## Usage: Rscript export_h5ad_hybrid.R [in_qs2] [out_h5ad] [count_assay = RNA | CB | RAW] [x_layer = counts | sct_data]

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
## RNA and CB are the same CellBender-corrected counts (SCT was fit on these); RAW is the uncorrected
## CellRanger count.
count_assay = if (length(args) >= 3) args[3] else "RNA"
## "counts" (default, what SCENIC+ needs) or "sct_data" (legacy: SCT log data, only to rebuild the earlier input)
x_layer = if (length(args) >= 4) args[4] else "counts"
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

## Raw integer counts, not SCT output. anndata.ipynb does `adata.raw = adata` *before* normalizing, and
## SCENIC+ reads `.raw` (prepare_GEX_ACC, use_raw_for_GEX_anndata=True) and normalizes internally, so X
## must be counts. Feeding it SCT data got log-normalized a second time (an earlier export did exactly that).
## Genes are restricted to the SCT feature set so the gene universe matches the earlier run.
stopifnot(count_assay %in% c("RNA", "CB", "RAW"), x_layer %in% c("counts", "sct_data"))
genes = rownames(obj[["SCT"]])
if (x_layer == "counts") {
  counts = LayerData(obj, assay = count_assay, layer = "counts")[genes, rownames(obs)]
  X = t(as(counts, "dgCMatrix"))  # cells x genes, sparse
  stopifnot(identical(rownames(X), rownames(obs)), all(X@x == round(X@x)))
} else {
  ## LEGACY, for reproducing the earlier run's input only (make_adata_doublenorm.py): SCT `data`, which is
  ## already log-normalized. Never use this for SCENIC+; see above.
  X = t(as(LayerData(obj, assay = "SCT", layer = "data")[genes, rownames(obs)], "dgCMatrix"))
  stopifnot(identical(rownames(X), rownames(obs)))
}

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
