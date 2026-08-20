## Prepare normalized objects for label transfer.
## Usage: Rscript slurm/prep.R --norm [lognorm|sct] --script-name [name]
## Run from the project root.

args = commandArgs(trailingOnly = TRUE)
get_arg = function(flag, default = NULL) {
  i = which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
norm        = get_arg("--norm",        "lognorm")
script_name = get_arg("--script-name", "snrna_zaremba-chk_integration-cca_glut")

library(Seurat)
library(qs)
library(here)
library(Matrix.utils)
library(Matrix)

source(here::here("config/paths.R"))

n_threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
RhpcBLASctl::blas_set_num_threads(n_threads)
RhpcBLASctl::omp_set_num_threads(n_threads)

out_dir         = here::here("snrna/integration", script_name)
norm_dir        = file.path(out_dir, norm)
bf_out_fname    = file.path(out_dir, "obj_bf.qs")
chk_out_fname   = file.path(out_dir, "obj_chk.qs")

obj_bf  = qread(bf_out_fname)
obj_chk = qread(chk_out_fname)

common_features = intersect(Features(obj_chk), Features(obj_bf))
obj_chk_common  = obj_chk[common_features,]
counts_bf       = as(GetAssayData(obj_bf, assay = "RNA", layer = "counts"), "dgCMatrix")
obj_bf_common   = CreateSeuratObject(counts = counts_bf[common_features, ],
                                     meta.data = obj_bf@meta.data)

if (norm == "lognorm") {

  DefaultAssay(obj_bf_common)  = "RNA"
  DefaultAssay(obj_chk_common) = "RNA"

  obj_bf_common  = obj_bf_common  %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  obj_chk_common = obj_chk_common %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  obj_chk_common = AddMetaData(obj_chk_common, FetchData(obj_chk, "anno_level_3"))

  qsave(obj_bf_common,  file.path(norm_dir, "obj_bf_common.qs"))
  qsave(obj_chk_common, file.path(norm_dir, "obj_chk_common.qs"))

} else if (norm == "sct") {

  obj_bf_common  = obj_bf_common  %>% SCTransform(verbose = FALSE) %>% RunPCA(verbose = FALSE)
  obj_chk_common = obj_chk_common %>% SCTransform(verbose = FALSE) %>% RunPCA(verbose = FALSE)

  sct_features   = SelectIntegrationFeatures(list(obj_bf_common, obj_chk_common), nfeatures = 3000)
  obj_bf_common  = PrepSCTIntegration(list(obj_bf_common),  anchor.features = sct_features)[[1]]
  obj_chk_common = PrepSCTIntegration(list(obj_chk_common), anchor.features = sct_features)[[1]]
  obj_chk_common = AddMetaData(obj_chk_common, FetchData(obj_chk, "anno_level_3"))

  qsave(obj_bf_common,  file.path(norm_dir, "obj_bf_common.qs"))
  qsave(obj_chk_common, file.path(norm_dir, "obj_chk_common.qs"))
  saveRDS(sct_features, file.path(norm_dir, "sct_features.rds"))

} else {
  stop("--norm must be 'lognorm' or 'sct'")
}

cat("prep.R done:", norm, script_name, "\n")
