## Run CCA integration via IntegrateLayers for a single parameter set.
## Usage: Rscript slurm/integrate.R --norm [lognorm|sct] --script-name [name]
##        Task index taken from SLURM_ARRAY_TASK_ID (1-based).
## Saves UMAP embeddings + metadata as a compact data frame per task.

args = commandArgs(trailingOnly = TRUE)
get_arg = function(flag, default = NULL) {
  i = which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
norm        = get_arg("--norm",        "lognorm")
script_name = get_arg("--script-name", "snrna_zaremba-chk_integration-cca_glut")
task_id     = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))

library(Seurat)
library(tidyverse)
library(qs)
library(here)

source(here::here("snrna/integration/slurm/config.R"))

n_threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
RhpcBLASctl::blas_set_num_threads(n_threads)
RhpcBLASctl::omp_set_num_threads(n_threads)
RcppParallel::setThreadOptions(numThreads = n_threads)

## IntegrateLayers uses future internally; raise the global size limit so large
## Seurat objects are not rejected by the default 500 MiB cap.
options(future.globals.maxSize = Inf)

norm_dir = here::here("snrna/integration", script_name, norm)
int_dir  = file.path(norm_dir, "integrated")
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)

i           = task_id
dims_to_use = transfer_params$dims_to_use[[i]]
ndims       = max(dims_to_use)
reduction   = transfer_params$reduction[i]
int_method  = if (reduction == "cca") CCAIntegration else RPCAIntegration
new_reduc   = paste0("integrated.", reduction)

cat(sprintf("Task %d: k.anchor=%s k.filter=%s k.score=%s max.features=%s dims=1:%s reduction=%s\n",
            i, transfer_params$k.anchor[i], transfer_params$k.filter[i],
            transfer_params$k.score[i], transfer_params$max.features[i], ndims, reduction))

out_fname = file.path(int_dir,
                      sprintf("umap_%04d_k%s_f%s_s%s_m%s_d%s_%s.rds",
                              i,
                              transfer_params$k.anchor[i],   transfer_params$k.filter[i],
                              transfer_params$k.score[i],    transfer_params$max.features[i],
                              ndims, reduction))

if (file.exists(out_fname)) {
  cat("Already exists, skipping.\n")
  quit(save = "no")
}

obj_bf  = qread(file.path(norm_dir, "obj_bf_common.qs"))
obj_chk = qread(file.path(norm_dir, "obj_chk_common.qs"))

obj = merge(obj_bf, obj_chk)
obj[["RNA"]] = JoinLayers(obj[["RNA"]])

if (norm == "lognorm") {
  obj[["RNA"]] = split(obj[["RNA"]], f = obj$dataset, layers = c("counts", "data"))
  obj = obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = transfer_params$max.features[i], verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = ndims, verbose = FALSE)

  obj = IntegrateLayers(obj,
                        method         = int_method,
                        orig.reduction = "pca",
                        new.reduction  = new_reduc,
                        dims           = dims_to_use,
                        k.anchor       = transfer_params$k.anchor[i],
                        k.filter       = transfer_params$k.filter[i],
                        k.score        = transfer_params$k.score[i],
                        verbose        = FALSE)
} else {
  ## Use RNA counts for SCT integration — the pre-applied SCT assay from prep
  ## is v3 format and loses variable features on merge. Run SCTransform fresh
  ## on the merged split object as required by Seurat v5 IntegrateLayers.
  DefaultAssay(obj) = "RNA"
  obj[["SCT"]]      = NULL
  obj[["RNA"]]      = split(obj[["RNA"]], f = obj$dataset, layers = c("counts", "data"))
  obj = obj %>%
    SCTransform(verbose = FALSE) %>%
    RunPCA(npcs = ndims, verbose = FALSE)

  obj = IntegrateLayers(obj,
                        method               = int_method,
                        orig.reduction       = "pca",
                        new.reduction        = new_reduc,
                        normalization.method = "SCT",
                        dims                 = dims_to_use,
                        k.anchor             = transfer_params$k.anchor[i],
                        k.filter             = transfer_params$k.filter[i],
                        k.score              = transfer_params$k.score[i],
                        verbose              = FALSE)
}

obj = RunUMAP(obj, reduction = new_reduc, dims = dims_to_use, verbose = FALSE)

umap_df = as.data.frame(Embeddings(obj, "umap")) %>%
  rownames_to_column("cell") %>%
  mutate(dataset      = obj$dataset,
         cluster      = obj$cluster,
         anno_level_3 = obj$anno_level_3,
         k.anchor     = transfer_params$k.anchor[i],
         k.filter     = transfer_params$k.filter[i],
         k.score      = transfer_params$k.score[i],
         max.features = transfer_params$max.features[i],
         ndims        = ndims,
         reduction    = reduction)

saveRDS(umap_df, out_fname)
cat("Done task", i, "\n")
