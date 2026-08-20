#!/usr/bin/env Rscript
## Prints a SLURM array spec for chunks that have no result yet, so a partially
## completed run can be resubmitted without redoing finished work.
## Usage: sbatch --array=$(Rscript 04_missing.R <hpc_rctd_proseg_dir>) scripts/rctd_array.sbatch <hpc_rctd_proseg_dir>
args <- commandArgs(trailingOnly = TRUE)
BASE <- args[1]
m <- readRDS(file.path(BASE, "manifest.rds"))
f <- list.files(file.path(BASE, "results"), pattern = "^result_.*qs2$")
have <- as.integer(sub("result_(\\d+).*", "\\1", f))
miss <- setdiff(seq_len(m$n_chunks), have)
if (!length(miss)) { message("nothing missing"); cat("\n"); quit(save = "no", status = 1) }
cat(paste(miss, collapse = ","), "\n")
