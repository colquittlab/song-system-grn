#!/usr/bin/env Rscript
## Step 3 — combine array results into one table. Run on the cluster after the
## array finishes, or locally after pulling `results/` back.
## Reports any missing chunks rather than silently returning a short table.
##
## Usage: Rscript 03_gather.R <hpc_rctd_dir> [out.csv.gz]

suppressMessages({ library(tidyverse); library(qs2) })

args <- commandArgs(trailingOnly = TRUE)
BASE <- args[1]
OUT  <- if (length(args) > 1) args[2] else file.path(BASE, "rctd_all.csv.gz")

manifest <- readRDS(file.path(BASE, "manifest.rds"))
files <- list.files(file.path(BASE, "results"), pattern = "^result_.*qs2$", full.names = TRUE)
have <- as.integer(str_match(basename(files), "result_(\\d+)")[, 2])
missing <- setdiff(seq_len(manifest$n_chunks), have)

cat("chunks expected:", manifest$n_chunks, " present:", length(have), "\n")
if (length(missing)) {
  cat("MISSING:", paste(missing, collapse = ","), "\n")
  cat("resubmit with: sbatch --array=", paste(missing, collapse = ","),
      " rctd_array.sbatch ", BASE, "\n", sep = "")
} else cat("all chunks present\n")

res <- map_dfr(files, qs_read)
cat("cells:", nrow(res), "\n\n")
print(table(res$spot_class))
cat("\ntop first_type:\n")
print(head(sort(table(res$first_type), decreasing = TRUE), 15))

write_csv(res, OUT)
cat("\nwritten:", OUT, "\n")
