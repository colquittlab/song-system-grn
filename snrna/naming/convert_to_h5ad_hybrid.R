## Convert the hybrid-labelled object to h5ad for cellxgene / cellxgene VIP display.
##
## Follows the convention of snrna/clustering/convert_to_h5ad_snrna.R (write via anndata/reticulate
## rather than sceasy/zellkonverter, which are incompatible with Seurat v5 objects and produce obsm
## that cellxgene will not read), with these changes forced by this particular object and by the
## display requirement:
##
##  1. MAGIC runs on the SCT `data` layer (log1p of SCT-corrected counts), not on a manually
##     CP10K-normalised RNA layer, so the smoothed expression matches the MAGIC_SCT assay used
##     elsewhere for display (e.g. snrna/deg/song-surround_deg_glut.qmd plots MAGIC_SCT alongside
##     SCT). The object already carries a precomputed SCT assay with `counts`/`data`/`scale.data`,
##     unlike RNA (which has only `counts`), so no manual normalisation is needed here.
##  2. SCT is a strict gene subset of RNA (15,523 of RNA's 17,029 genes; SCTransform's default
##     min-cells filter drops the rest), so the exported var list is SCT's gene set. The raw-counts
##     layer in the h5ad is RNA counts subset to those same genes, kept for true, unsmoothed UMI
##     counts alongside the SCT-normalised and MAGIC-smoothed values -- SCT's own `counts` layer is
##     "corrected" counts, not raw UMIs, so it is not used for that layer.
##  3. The RNA counts layer is a BPCells on-disk IterableMatrix (class `RenameDims`), not an
##     in-memory matrix. as.matrix() on it is not the same operation as on a dgCMatrix, so it is
##     materialised explicitly first. (Same gotcha that broke the label-transfer notebook.) The SCT
##     layers are already in-memory dgCMatrix and need no such step.
##  4. 659 of 34,295 cells have celltype_hybrid = NA -- the clusters excluded by the hybrid scheme
##     (Glut-Arco-4, Glut-Nido-3, confirmed artifacts). They are KEPT here, so this file is a faithful
##     conversion of the object rather than a silently filtered subset; the NA is preserved as a
##     missing value in obs. Filter downstream (or set drop_unlabelled = TRUE below) if you want only
##     labelled cells.
##  5. The embedding is taken from the parameter-sweep object written by
##     snrna/reduction_viz/combined_all_umap_hybrid.R, reduction `dims40nn30mindist0.3`, NOT the
##     labels object's own `umap_rna`. That script drops the clusters the hybrid scheme excludes
##     before fitting the UMAP, so its embedding covers 33,636 of the 34,295 cells; the h5ad is
##     therefore restricted to those cells, because cellxgene needs coordinates for every cell it
##     displays. That makes point 4 moot in practice -- the excluded cells are absent rather than
##     NA-labelled -- but the flag is kept for the case of exporting against a different embedding.
##  6. Expression is MAGIC-smoothed for display. `X` holds the smoothed values, so colouring by gene
##     in cellxgene/VIP is denoised by default. The unsmoothed data is not thrown away: layers/counts
##     holds raw RNA counts (subset to SCT's genes) and layers/sct_data the SCT `data` values MAGIC
##     was run on, both sparse, so any quantitative use can switch layers in VIP.
library(Seurat)
library(purrr)
library(qs2)
library(Matrix)
library(reticulate)
library(tidyverse)

in_dir = "/ssd/brad/rstudio/multiome/song-system-grn/snrna/naming/hybrid_division_naming"
infile_obj = file.path(in_dir, "obj_hybrid_labels.qs2")

out_dir = in_dir
h5ad_outfile = file.path(out_dir, "obj_hybrid_labels.h5ad")

## Written to a temp path and renamed into place only after the round-trip checks pass: an h5ad
## write that fails partway (e.g. HDF5 refusing to lock the file while a cellxgene server has it
## open) otherwise truncates the previous good output to zero bytes.
h5ad_tmpfile = paste0(h5ad_outfile, ".tmp")

## Embedding source: the UMAP parameter sweep, whose cell set defines the cells exported here.
umap_obj_file = "/ssd/brad/rstudio/multiome/song-system-grn/snrna/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2"
umap_reduction = "dims40nn30mindist0.3"

drop_unlabelled = FALSE

## MAGIC parameters. Defaults of the magic-impute package (knn = 5, decay = 1, t = "auto"); only the
## thread count and solver are pinned. "exact" rather than "approximate": at 34k cells the exact
## diffusion operator is affordable here (~950 GB RAM free) and avoids the landmark approximation.
magic_n_jobs = 16L
magic_solver = "exact"

## Read the sweep object first and keep only its embedding, so both 3 GB objects are never resident
## at once.
umap_obj = qs_read(umap_obj_file, nthreads = 8)
stopifnot(
  "requested UMAP reduction not in the sweep object" =
    umap_reduction %in% Reductions(umap_obj)
)
umap_emb = Embeddings(umap_obj, reduction = umap_reduction)
rm(umap_obj); gc(verbose = FALSE)
message("Embedding ", umap_reduction, ": ", nrow(umap_emb), " cells")
stopifnot("embedding contains NA coordinates" = !anyNA(umap_emb))

obj = qs_read(infile_obj, nthreads = 8)

assay_to_export = "SCT"
DefaultAssay(obj) = assay_to_export
if (any(grepl("\\.[0-9]+", Layers(obj))))
  obj = JoinLayers(obj)

if (drop_unlabelled) {
  keep = colnames(obj)[!is.na(obj$celltype_hybrid)]
  message("Dropping ", ncol(obj) - length(keep), " cells with no hybrid label")
  obj = subset(obj, cells = keep)
}

## Restrict to the cells the embedding covers, in the object's own order.
stopifnot(
  "embedding has cells absent from the labels object" =
    all(rownames(umap_emb) %in% colnames(obj))
)
keep_umap = intersect(colnames(obj), rownames(umap_emb))
message("Restricting to the ", length(keep_umap), " cells in ", umap_reduction,
        " (dropping ", ncol(obj) - length(keep_umap), ")")
obj = subset(obj, cells = keep_umap)
umap_emb = umap_emb[colnames(obj), , drop = FALSE]

## Raw RNA counts, subset to the SCT gene set, kept as the h5ad's unsmoothed counts layer. Pulled
## before DietSeurat drops the RNA assay. RNA counts is a BPCells on-disk IterableMatrix (class
## `RenameDims`), not an in-memory matrix, so it is materialised explicitly first.
raw_counts = as(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"), "dgCMatrix")
stopifnot(
  "SCT genes are not a subset of RNA genes" = all(rownames(obj[["SCT"]]) %in% rownames(raw_counts))
)
raw_counts = raw_counts[rownames(obj[["SCT"]]), , drop = FALSE]

## No reductions carried over: the embedding comes from the sweep object, not this one.
obj_slim = DietSeurat(
  obj,
  layers = c("counts", "data"),
  assays = assay_to_export,
  dimreducs = NULL
)
rm(obj); gc(verbose = FALSE)
raw_counts = raw_counts[, colnames(obj_slim), drop = FALSE]

## Metadata carried over: the new labels, the original cluster they came from (so the relabelling is
## traceable in the h5ad itself), dissection position, and the library/individual + QC columns.
md = obj_slim@meta.data %>%
  rownames_to_column() %>%
  select(rowname, celltype_hybrid, cluster, position, sample_name, assignment,
         nCount_RNA, nFeature_RNA, perc_mito) %>%
  column_to_rownames()
obj_slim@meta.data = md

stopifnot(
  "metadata rows and object cells are out of sync" =
    identical(rownames(obj_slim@meta.data), colnames(obj_slim)),
  "embedding rows and object cells are out of sync" =
    identical(rownames(umap_emb), colnames(obj_slim))
)

# Write h5ad directly with anndata — avoids sceasy/zellkonverter Seurat v5 incompatibilities
# and guarantees cellxgene-compatible obsm format
use_condaenv("/home/brad/micromamba/envs/magic")
## convert = FALSE throughout: the MAGIC output is a 34k x 17k dense array and reticulate's default
## auto-conversion would materialise it as an R double matrix (~4.7 GB, and then float64 in the
## h5ad). Keeping it as a python object hands it straight to anndata as float32.
ad <- reticulate::import("anndata", convert = FALSE)
np <- reticulate::import("numpy", convert = FALSE)
sp <- reticulate::import("scipy.sparse", convert = FALSE)
mg <- reticulate::import("magic", convert = FALSE)

sct_data = SeuratObject::LayerData(obj_slim, assay = assay_to_export, layer = "data")
stopifnot(
  "raw_counts and sct_data are out of sync" =
    identical(dim(raw_counts), dim(sct_data)) && identical(colnames(raw_counts), colnames(sct_data))
)

## Sparse CSR rather than the template's dense np$array: this matrix is 15,523 x 34,295, which is
## ~2 GB dense in float32, versus a small fraction of that sparse. A genes x cells dgCMatrix is
## already CSC, whose slots are exactly the CSR slots of its cells x genes transpose, so this is a
## relabelling rather than a transpose.
as_csr_cells_by_genes = function(m) {
  sp$csr_matrix(
    reticulate::tuple(
      np$array(m@x, dtype = "float32"),
      np$array(m@i, dtype = "int32"),
      np$array(m@p, dtype = "int32")
    ),
    shape = reticulate::tuple(ncol(m), nrow(m))
  )
}

X_counts = as_csr_cells_by_genes(raw_counts)
X_sct_data = as_csr_cells_by_genes(sct_data)
rm(sct_data); gc(verbose = FALSE)

## ---- MAGIC ----------------------------------------------------------------------------------
## genes = "all_genes": every gene is imputed, so any gene is colourable in cellxgene rather than
## only a preselected panel. Output is dense cells x genes (~2 GB float32). Input is the SCT `data`
## layer (log1p of SCT-corrected counts), matching the MAGIC_SCT assay used elsewhere for display.
message("Running MAGIC on ", ncol(raw_counts), " cells x ", nrow(raw_counts), " genes ...")
magic_op = mg$MAGIC(n_jobs = magic_n_jobs, solver = magic_solver, verbose = 1L)
X_magic = magic_op$fit_transform(X_sct_data, genes = "all_genes")
X_magic = np$asarray(X_magic, dtype = "float32")
message("MAGIC done; smoothed matrix shape ",
        paste(unlist(reticulate::py_to_r(X_magic$shape)), collapse = " x "))

umap_mat <- np$array(umap_emb, dtype = "float32")

adata <- ad$AnnData(
  X = X_magic,
  obs = reticulate::r_to_py(obj_slim@meta.data),
  var = reticulate::r_to_py(data.frame(row.names = rownames(obj_slim))),
  obsm = list(X_umap_seurat = umap_mat),
  layers = list(counts = X_counts, sct_data = X_sct_data)
)
adata$uns = reticulate::r_to_py(list(
  X_description = "MAGIC-smoothed SCT data (log1p of SCT-corrected counts) expression",
  umap_reduction = umap_reduction,
  umap_source = umap_obj_file,
  magic_params = list(
    input_layer = "sct_data",
    normalisation = "SCTransform",
    genes = "all_genes",
    solver = magic_solver,
    n_jobs = magic_n_jobs,
    magic_version = reticulate::py_to_r(mg$`__version__`)
  )
))
if (file.exists(h5ad_tmpfile)) unlink(h5ad_tmpfile)
adata$write_h5ad(h5ad_tmpfile)

## Read back and check the round trip actually preserved shape, labels, embedding and the smoothed
## and unsmoothed matrices, rather than assuming the write succeeded because it did not error.
## The checks run on the python side: X is a 34k x 17k dense float32 array and pulling it into R
## would double it to ~4.7 GB for no reason.
obj_h5ad = ad$read_h5ad(h5ad_tmpfile)
py$adata_check = obj_h5ad
st = reticulate::py_run_string('
import numpy as np
_a = adata_check
_check_genes = [g for g in ["GAD1", "SLC17A6", "SLC17A7", "AQP4", "PVALB"] if g in _a.var_names]
_cors = {}
for _g in _check_genes:
    _j = _a.var_names.get_loc(_g)
    _ln = np.asarray(_a.layers["sct_data"][:, _j].todense()).ravel()
    _mg = np.asarray(_a.X[:, _j]).ravel()
    _cors[_g] = float(np.corrcoef(_ln, _mg)[0, 1]) if _ln.std() > 0 else float("nan")
_stats = {
    "n_obs": int(_a.n_obs),
    "n_vars": int(_a.n_vars),
    "n_nan": int(np.isnan(_a.X).sum()),
    "x_min": float(np.nanmin(_a.X)),
    "x_max": float(_a.X.max()),
    "frac_nz_magic": float((_a.X > 0).mean()),
    "frac_nz_sct_data": float(_a.layers["sct_data"].nnz) / (_a.n_obs * _a.n_vars),
    "layers": list(_a.layers.keys()),
    "obs_cols": list(map(str, _a.obs.columns)),
    "obsm_keys": list(_a.obsm.keys()),
    "gene_cors": _cors,
}
')$`_stats`

cat("wrote:", h5ad_outfile, "\n")
cat("shape:", st$n_obs, "x", st$n_vars,
    " (expected", ncol(obj_slim), "x", nrow(obj_slim), ")\n")
cat("obs cols:", paste(st$obs_cols, collapse = ", "), "\n")
cat("obsm:", paste(unlist(st$obsm_keys), collapse = ", "), "\n")
cat("layers:", paste(unlist(st$layers), collapse = ", "), "\n")
cat("NaNs in X:", st$n_nan, "\n")
cat("X (MAGIC) min/max:", round(st$x_min, 4), "/", round(st$x_max, 4), "\n")
cat("fraction nonzero -- X (MAGIC):", round(st$frac_nz_magic, 3),
    " sct_data:", round(st$frac_nz_sct_data, 3), "\n")
if (length(st$gene_cors) > 0)
  cat("per-gene cor(sct_data, MAGIC):",
      paste(sprintf("%s=%.2f", names(st$gene_cors), unlist(st$gene_cors)), collapse = ", "), "\n")

stopifnot(
  "cell count changed in the round trip" = st$n_obs == ncol(obj_slim),
  "gene count changed in the round trip" = st$n_vars == nrow(obj_slim),
  "X is not MAGIC-smoothed: no denser than the sct_data input" =
    st$frac_nz_magic > st$frac_nz_sct_data,
  "MAGIC produced negative expression" = st$x_min >= 0,
  "MAGIC produced NaN expression values" = st$n_nan == 0,
  "counts/sct_data layers missing" = all(c("counts", "sct_data") %in% unlist(st$layers)),
  "umap embedding missing from obsm" = "X_umap_seurat" %in% unlist(st$obsm_keys),
  "celltype_hybrid missing from obs" = "celltype_hybrid" %in% st$obs_cols
)
## Only now is the previous output replaced.
try(obj_h5ad$file$close(), silent = TRUE)  # in-memory read, but release any handle
rm(obj_h5ad); py$adata_check = NULL
invisible(gc(verbose = FALSE))
stopifnot("could not move the validated h5ad into place" =
            file.rename(h5ad_tmpfile, h5ad_outfile))
cat("OK: MAGIC-smoothed h5ad ready for cellxgene VIP ->", h5ad_outfile, "\n")
