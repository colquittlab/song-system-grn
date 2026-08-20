## Run label transfer for a single parameter set (one Slurm array task).
## Usage: Rscript slurm/transfer.R --norm [lognorm|sct] --script-name [name]
##        Task index taken from SLURM_ARRAY_TASK_ID (1-based).

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

out_dir     = here::here("snrna/integration", script_name)
norm_dir    = file.path(out_dir, norm)
anchors_dir = file.path(norm_dir, "anchors")
dir.create(anchors_dir, recursive = TRUE, showWarnings = FALSE)

obj_bf_common  = qread(file.path(norm_dir, "obj_bf_common.qs"))
obj_chk_common = qread(file.path(norm_dir, "obj_chk_common.qs"))

if (norm == "sct") {
  sct_features = readRDS(file.path(norm_dir, "sct_features.rds"))
}

i           = task_id
dims_to_use = transfer_params$dims_to_use[[i]]
ndims       = max(dims_to_use)
reduction   = transfer_params$reduction[i]

out_fname = file.path(norm_dir,
                      sprintf("pred_long_k%s_f%s_s%s_m%s_d%s_%s.rds",
                              transfer_params$k.anchor[i],   transfer_params$k.filter[i],
                              transfer_params$k.score[i],    transfer_params$max.features[i],
                              ndims, reduction))
if (file.exists(out_fname)) {
  cat("Task", i, "already complete, skipping.\n")
  quit(save = "no")
}

cat(sprintf("Task %d: k.anchor=%s k.filter=%s k.score=%s max.features=%s dims=1:%s reduction=%s\n",
            i, transfer_params$k.anchor[i], transfer_params$k.filter[i],
            transfer_params$k.score[i], transfer_params$max.features[i], ndims, reduction))

anchor_fname = file.path(anchors_dir,
                         sprintf("anchor_k%s_f%s_s%s_m%s_d%s_%s.qs",
                                 transfer_params$k.anchor[i], transfer_params$k.filter[i],
                                 transfer_params$k.score[i],  transfer_params$max.features[i],
                                 ndims, reduction))

## weight.reduction for TransferData: CCA uses the CCA embedding;
## RPCA uses the projected PCA of the query ("pcaproject")
weight_reduction = if (reduction == "cca") "cca" else "pcaproject"

if (file.exists(anchor_fname)) {
  anchors = qread(anchor_fname)
} else if (norm == "lognorm") {
  anchors = FindTransferAnchors(reference = obj_bf_common,
                                normalization.method = "LogNormalize",
                                k.anchor     = transfer_params$k.anchor[i],
                                k.filter     = transfer_params$k.filter[i],
                                k.score      = transfer_params$k.score[i],
                                max.features = transfer_params$max.features[i],
                                reduction    = reduction,
                                query        = obj_chk_common,
                                dims         = dims_to_use)
  qsave(anchors, anchor_fname)
} else {
  anchors = FindTransferAnchors(reference = obj_bf_common,
                                normalization.method = "SCT",
                                features  = sct_features,
                                k.anchor     = transfer_params$k.anchor[i],
                                k.filter     = transfer_params$k.filter[i],
                                k.score      = transfer_params$k.score[i],
                                max.features = transfer_params$max.features[i],
                                reduction    = reduction,
                                query        = obj_chk_common,
                                dims         = dims_to_use)
  qsave(anchors, anchor_fname)
}

predictions_cluster = TransferData(anchorset        = anchors,
                                   weight.reduction = weight_reduction,
                                   refdata          = obj_bf_common$cluster,
                                   dims             = dims_to_use) %>%
  mutate(predicted.id_cluster = predicted.id)

predictions_cluster$anno_level_3 = obj_chk_common$anno_level_3

pred_long_i = predictions_cluster %>%
  dplyr::select(contains("prediction.score"), "anno_level_3") %>%
  dplyr::select(-prediction.score.max) %>%
  pivot_longer(-anno_level_3, names_to = "col_tmp") %>%
  mutate(col1 = str_extract(col_tmp, "Glut.*$")) %>%
  mutate(cluster = gsub("\\.", "-", col1)) %>%
  mutate(k.anchor     = transfer_params$k.anchor[i],
         k.filter     = transfer_params$k.filter[i],
         k.score      = transfer_params$k.score[i],
         max.features = transfer_params$max.features[i],
         ndims        = ndims,
         reduction    = reduction)

saveRDS(pred_long_i, out_fname)
cat("Done task", i, "\n")
