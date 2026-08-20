## Aggregate per-task pred_long files into a single pred_long.rds.
## Usage: Rscript slurm/aggregate.R --norm [lognorm|sct] --script-name [name]

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

norm_dir = here::here("snrna/integration", script_name, norm)
nparams  = nrow(transfer_params)

fnames = file.path(norm_dir,
                   sprintf("pred_long_k%s_f%s_s%s_m%s_d%s_%s.rds",
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

pred_long = map(fnames, readRDS) %>% bind_rows()
saveRDS(pred_long, file.path(norm_dir, "pred_long.rds"))
cat("Aggregated", nparams, "parameter sets ->", file.path(norm_dir, "pred_long.rds"), "\n")
