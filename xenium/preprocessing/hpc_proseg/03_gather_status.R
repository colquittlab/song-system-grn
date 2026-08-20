#!/usr/bin/env Rscript
## Run after pulling results/ back (or directly on the cluster). Reports which
## sections finished and which are missing, rather than silently proceeding
## with a partial set. Doesn't merge anything across sections — proseg output
## is per-section by design (see README) — this just checks completeness.
##
## Usage: Rscript 03_gather_status.R <hpc_proseg_dir>

suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
BASE <- args[1]

manifest <- read_tsv(file.path(BASE, "manifest.tsv"), show_col_types = FALSE)

status <- manifest %>%
  mutate(result_dir = file.path(BASE, "results", section_id),
         done = file.exists(file.path(result_dir, "cell-metadata.parquet")),
         n_transcript_meta = file.exists(file.path(result_dir, "transcript-metadata.parquet")),
         n_counts = file.exists(file.path(result_dir, "expected-counts.csv.gz")))

print(as.data.frame(status %>% select(idx, section_id, done, n_transcript_meta, n_counts)),
      row.names = FALSE)

missing <- status %>% filter(!done) %>% pull(idx)
if (length(missing)) {
  cat("\nMISSING:", paste(missing, collapse = ","), "\n")
  cat("resubmit with: sbatch --array=", paste(missing, collapse = ","),
      " scripts/proseg_array.sbatch ", BASE, "\n", sep = "")
} else {
  cat("\nall", nrow(status), "sections present\n")
}
