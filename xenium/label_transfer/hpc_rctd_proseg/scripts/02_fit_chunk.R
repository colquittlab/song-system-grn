#!/usr/bin/env Rscript
## Step 2 (RUNS ON HPC, one SLURM array task per chunk) — fit one pre-sliced
## RCTD chunk. Idempotent: an existing output is left alone, so a partially
## completed array can be resubmitted with --array over the missing indices.
##
## Usage: Rscript 02_fit_chunk.R <chunk_dir> <out_dir> <index> [cores]

suppressMessages({
  library(tidyverse); library(qs2); library(Matrix); library(spacexr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: 02_fit_chunk.R <chunk_dir> <out_dir> <index> [cores]")
CHUNK_DIR <- args[1]; OUT_DIR <- args[2]
IDX <- as.integer(args[3])
CORES <- if (length(args) > 3) as.integer(args[4]) else 4

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
in_f  <- file.path(CHUNK_DIR, sprintf("chunk_%03d.qs2", IDX))
out_f <- file.path(OUT_DIR,  sprintf("result_%03d.qs2", IDX))

manifest <- readRDS(file.path(dirname(CHUNK_DIR), "manifest.rds"))

## A stale --array upper bound (hardcoded in the sbatch script, only correct
## if manually reconciled against manifest.rds$n_chunks after 01_prepare_chunks.R
## runs) previously caused a guaranteed one-task failure per submission on any
## dataset whose chunk count differs from the number carried over from a prior
## run. An index beyond n_chunks means the array range is just wider than the
## data — exit cleanly rather than erroring, so a stale range reads as "no-op",
## not "failure".
if (IDX > manifest$n_chunks) {
  message("chunk ", IDX, " exceeds n_chunks (", manifest$n_chunks,
          ") — array range wider than the data, exiting cleanly")
  quit(save = "no", status = 0)
}

if (file.exists(out_f)) { message("chunk ", IDX, " already done — exiting"); quit(save = "no") }
if (!file.exists(in_f)) {
  stop("chunk ", IDX, " is within range (n_chunks=", manifest$n_chunks,
       ") but missing input: ", in_f)
}

t0 <- Sys.time()
r <- qs_read(in_f)
## spacexr scales poorly within a process (measured ~2.4 of 22 cores busy), so
## keep per-task cores modest and get parallelism from array width instead.
r@config$max_cores <- CORES
r <- fitPixels(r, doublet_mode = "doublet")

res <- r@results$results_df %>%
  rownames_to_column("cell") %>%
  as_tibble() %>%
  mutate(across(c(first_type, second_type),
                ~ unname(manifest$name_map[as.character(.x)])),
         chunk = IDX)

qs_save(res, out_f)
message("chunk ", IDX, ": ", nrow(res), " cells in ", format(Sys.time() - t0))
print(table(res$spot_class))
