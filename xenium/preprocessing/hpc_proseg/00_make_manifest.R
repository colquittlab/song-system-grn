#!/usr/bin/env Rscript
## Run LOCALLY. Enumerates the ten Xenium regions from experiment.xenium (the
## authoritative source for section identity — see xenium/README.md on why the
## delivered object/file names cannot be trusted), and stages a SELF-CONTAINED
## directory: manifest + transcripts + scripts all under one tree, so the
## whole thing can be moved with a single `rclone sync` where a direct
## rsync/ssh between machines isn't available (e.g. VPN restrictions blocking
## a direct HPC <-> workstation path).
##
## transcripts.parquet files are hardlinked into src/, not copied — same
## filesystem, so this costs no extra disk space or time, and a change to
## either path is invisible to the other (they share an inode).
##
## Usage: Rscript 00_make_manifest.R [outdir]

suppressMessages({ library(tidyverse); library(jsonlite); library(here) })
source(here::here("config/paths.R"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) > 0) args[1] else
  file.path(path.expand(XENIUM_PROC_DIR), "hpc_proseg")
dir.create(file.path(OUT, "src"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
## SLURM opens --output/--error under logs/ before the array's own script body
## runs, so that directory must already exist when sbatch is invoked — an
## in-script mkdir is too late for the very first task. rclone/rsync also
## commonly drop empty directories in transfer (rclone needs
## --create-empty-src-dirs to keep them), so an empty logs/ created here is not
## guaranteed to survive the sync. A placeholder file forces it to transfer as
## a non-empty directory regardless of rclone flags.
file.create(file.path(OUT, "logs", ".keep"))

region_dirs <- list.dirs(XENIUM_RAW_DIR, recursive = FALSE)
region_dirs <- region_dirs[grepl("^output-XETG", basename(region_dirs))]

sections <- map_dfr(region_dirs, function(d) {
  e <- jsonlite::fromJSON(file.path(d, "experiment.xenium"))
  tibble(dir = d, bird = e$cassette_name, section = e$region_name)
}) %>%
  mutate(section_id = paste(bird, section, sep = "_")) %>%
  arrange(bird, section) %>%
  mutate(idx = row_number(),
         src_transcripts = file.path(dir, "transcripts.parquet"),
         ## relative to OUT, so the manifest is portable across machines
         transcripts_path = file.path("src", section_id, "transcripts.parquet"))

stopifnot(all(file.exists(sections$src_transcripts)))

for (i in seq_len(nrow(sections))) {
  dest_dir <- file.path(OUT, "src", sections$section_id[i])
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(OUT, sections$transcripts_path[i])
  if (!file.exists(dest)) {
    ok <- file.link(sections$src_transcripts[i], dest)
    if (!ok) file.copy(sections$src_transcripts[i], dest)  # fallback across filesystems
  }
}

manifest <- sections %>% select(idx, section_id, transcripts_path)
write_tsv(manifest, file.path(OUT, "manifest.tsv"))

## copy the SLURM/array scripts + conda env spec into the same tree — see
## step 2 in README
script_dir <- here::here("xenium/preprocessing/hpc_proseg")
file.copy(file.path(script_dir, c("proseg_array.sbatch", "03_gather_status.R",
                                  "environment.yml")),
          file.path(OUT, "scripts"), overwrite = TRUE)

cat(nrow(manifest), "sections staged under", OUT, "\n")
print(manifest %>% select(idx, section_id, transcripts_path))

cat("\nTotal transcripts.parquet size:",
    round(sum(file.size(file.path(OUT, manifest$transcripts_path))) / 1e9, 1), "GB\n")
cat("(hardlinked, not duplicated on disk — `du -sh", OUT, "` will not show this twice)\n")
cat("\nEverything needed is now under:", OUT, "\n")
cat("rclone it as one tree, e.g.:\n")
cat("  rclone sync ", OUT, " remote:hpc_proseg\n", sep = "")
