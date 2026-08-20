## Aggregate per-task UMAP data frames from integrate.R into umap_long.rds.
## Usage: Rscript slurm/aggregate_integrate.R --norm [lognorm|sct] --script-name [name]

args = commandArgs(trailingOnly = TRUE)
get_arg = function(flag, default = NULL) {
  i = which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
norm        = get_arg("--norm",        "lognorm")
script_name = get_arg("--script-name", "snrna_zaremba-chk_integration-cca_glut")

library(tidyverse)
library(here)

source(here::here("snrna/integration/slurm/config.R"))

int_dir  = here::here("snrna/integration", script_name, norm, "integrated")
nparams  = nrow(transfer_params)

fnames = file.path(int_dir,
                   sprintf("umap_%04d_k%s_f%s_s%s_m%s_d%s_%s.rds",
                           seq_len(nparams),
                           transfer_params$k.anchor,
                           transfer_params$k.filter,
                           transfer_params$k.score,
                           transfer_params$max.features,
                           map_int(transfer_params$dims_to_use, max),
                           transfer_params$reduction))

missing = fnames[!file.exists(fnames)]
if (length(missing) > 0) {
  stop("Missing result files:\n", paste(missing, collapse = "\n"))
}

umap_long = map(fnames, readRDS) %>% bind_rows()
out_fname = file.path(here::here("snrna/integration", script_name, norm), "umap_long.rds")
saveRDS(umap_long, out_fname)
cat("Aggregated", nparams, "parameter sets ->", out_fname, "\n")
