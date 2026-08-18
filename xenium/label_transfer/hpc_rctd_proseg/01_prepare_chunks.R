#!/usr/bin/env Rscript
## Step 1 (RUN LOCALLY) — RCTD label transfer over the proseg-resegmented
## Xenium data (all ten sections), mirroring xenium/label_transfer/hpc/ for the
## delivered 5um-expansion segmentation. Builds the setup once, pre-slices into
## per-chunk objects for a SLURM array. See that directory's README for the
## sizing rationale (choose_sigma_c is a fixed cost best run on few cores;
## fitPixels scales and should be chunked at ~20k cells).
##
## Difference from the delivered-segmentation run: the query is built from
## proseg's per-section `expected-counts.csv.gz` + `cell-metadata.parquet`
## rather than the merged `obj_xenium.qs2`, and proseg's cell IDs restart at 0
## in every section, so cells are renamed `<section_id>_proseg_<cell>` to stay
## globally unique.
##
## A real gotcha found while building this: proseg writes its gene panel to
## each section's zarr store in a DIFFERENT order per section (confirmed —
## OR52YW26_1_4 starts GAD1/SLC1A3/..., OR69PU4_2_7 starts SLC1A3/SNAP25/...).
## Row position is NOT safe to assume matches across sections; every section's
## matrix is re-indexed by gene NAME before combining. Gene order per section
## was extracted once (via a python env with zarr — proseg has no plain-text
## gene list output) into gene_lists/<section_id>.txt, read below.
##
## Usage:  Rscript 01_prepare_chunks.R [outdir]

suppressMessages({
  library(tidyverse); library(Seurat); library(BPCells); library(qs2)
  library(Matrix); library(spacexr); library(here)
})
source(here::here("config/paths.R"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) > 0) args[1] else
  file.path(path.expand(XENIUM_PROC_DIR), "hpc_rctd_proseg")
CHUNK <- 2e4
UMI_MIN <- 20
SETUP_CORES <- 4
MAX_REF_PER_LABEL <- 1000
EXCLUDE <- c("Glut-Nido-3")
ASSAY_REF <- "RAW"

dir.create(file.path(OUT, "chunks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
file.create(file.path(OUT, "logs", ".keep"))   # see hpc_proseg/README on why
t0 <- Sys.time()

sections <- c("OR52YW26_1_4", "OR52YW26_1_7", "OR52YW26_2_2", "OR52YW26_2_4", "OR52YW26_2_7",
             "OR69PU4_1_4", "OR69PU4_1_7", "OR69PU4_2_2", "OR69PU4_2_4", "OR69PU4_2_7")

## ---- reference (identical config to the delivered-segmentation run) -----
label_set <- read_csv(here::here("xenium/label_transfer/panel_resolution_benchmark",
                                 "recommended_label_set.csv"), show_col_types = FALSE)
panel <- read_csv(XENIUM_PANEL_CSV, show_col_types = FALSE)$Gene

o <- qs_read(file.path(SNRNA_SEURAT_DIR, "obj_clustered.qs2"))
o <- subset(o, cells = Cells(o)[!grepl("/", o$assignment)])
ca <- LayerData(o, assay = ASSAY_REF, layer = "counts")
ps <- intersect(panel, rownames(ca))
md <- o@meta.data %>% rownames_to_column("cell") %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(label_set %>% select(cluster, recommended_label), by = "cluster") %>%
  filter(!is.na(recommended_label), !recommended_label %in% EXCLUDE)
set.seed(10)
md <- md %>% group_by(recommended_label) %>% slice_sample(n = MAX_REF_PER_LABEL) %>% ungroup()

counts_ref <- as(ca[ps, md$cell], "dgCMatrix")
labels_rctd <- setNames(factor(md$recommended_label), md$cell)
levels(labels_rctd) <- str_replace_all(levels(labels_rctd), "/", "_")
name_map <- setNames(str_replace_all(levels(labels_rctd), "_", "/"), levels(labels_rctd))
rm(ca, o); gc()
message("reference: ", ncol(counts_ref), " cells, ", nlevels(labels_rctd), " labels, ",
        length(ps), " genes")

## ---- query: load each proseg section, re-index genes by NAME ------------
gene_list_dir <- here::here("xenium/label_transfer/hpc_rctd_proseg/gene_lists")

load_section <- function(sid) {
  d <- file.path(path.expand(XENIUM_PROSEG_DIR), sid)
  genes_sid <- readLines(file.path(gene_list_dir, paste0(sid, ".txt")))
  stopifnot(all(ps %in% genes_sid))   # every panel/reference gene must be present

  m <- Matrix::readMM(gzfile(file.path(d, "expected-counts.csv.gz")))
  m <- as(t(m), "CsparseMatrix")      # -> genes x cells
  rownames(m) <- genes_sid
  m <- m[ps, , drop = FALSE]          # re-index to the shared gene set, BY NAME
  m@x <- round(m@x); m <- Matrix::drop0(m)

  cmeta <- nanoparquet::read_parquet(file.path(d, "cell-metadata.parquet"))
  cell_names <- paste0(sid, "_proseg_", cmeta$cell)
  colnames(m) <- cell_names

  list(counts = m,
       meta = tibble(cell = cell_names, section_id = sid,
                     original_cell_id = cmeta$original_cell_id,
                     x_centroid = cmeta$centroid_x, y_centroid = cmeta$centroid_y,
                     volume = cmeta$volume))
}

message("loading ", length(sections), " proseg sections...")
per_section <- map(sections, load_section)
names(per_section) <- sections

counts_q_all <- do.call(cbind, map(per_section, "counts"))
mdx <- map_dfr(per_section, "meta")
rm(per_section); gc()
message("query: ", ncol(counts_q_all), " cells across ", length(sections), " sections; ",
        format(Sys.time() - t0))

write_csv(mdx, file.path(OUT, "proseg_cell_metadata.csv.gz"))

## ---- offset sections so they don't overlap in coordinate space ----------
## (doublet-mode assignment does not use coordinates; this only keeps the
## stored geometry sane, matching the delivered-segmentation run)
offs <- mdx %>% group_by(section_id) %>%
  summarize(w = max(x_centroid, na.rm = TRUE) - min(x_centroid, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_offset = lag(cumsum(w * 1.1), default = 0))
coords <- mdx %>% left_join(offs, by = "section_id") %>%
  transmute(cell, x = x_centroid + x_offset, y = y_centroid) %>% as.data.frame()
rownames(coords) <- coords$cell; coords$cell <- NULL
coords <- coords[colnames(counts_q_all), ]

## ---- setup: paid once ----------------------------------------------------
ref <- Reference(counts_ref, labels_rctd, colSums(counts_ref),
                 min_UMI = UMI_MIN, n_max_cells = MAX_REF_PER_LABEL)
sp_all <- SpatialRNA(coords, counts_q_all, colSums(counts_q_all))
r <- create.RCTD(sp_all, ref, max_cores = SETUP_CORES,
                 UMI_min = UMI_MIN, counts_MIN = UMI_MIN)
r@config$RCTDmode <- "doublet"
message("create.RCTD done: ", format(Sys.time() - t0))
r <- fitBulk(r)
r <- choose_sigma_c(r)
message("setup complete: ", format(Sys.time() - t0))

## ---- pre-slice ------------------------------------------------------------
fit_cells <- colnames(r@spatialRNA@counts)
chunks <- split(fit_cells, ceiling(seq_along(fit_cells) / CHUNK))
message("slicing ", length(chunks), " chunks of ", CHUNK)

for (i in seq_along(chunks)) {
  rc <- r
  rc@spatialRNA <- restrict_puck(r@spatialRNA, chunks[[i]])
  if (.hasSlot(rc, "originalSpatialRNA")) rc@originalSpatialRNA <- rc@spatialRNA
  qs_save(rc, file.path(OUT, "chunks", sprintf("chunk_%03d.qs2", i)))
  if (i %% 10 == 0) message("  sliced ", i, "/", length(chunks))
}

saveRDS(list(name_map = name_map, n_chunks = length(chunks),
             chunk_size = CHUNK, genes = ps,
             cells_per_chunk = lengths(chunks)),
        file.path(OUT, "manifest.rds"))

## copy the SLURM/array scripts + env spec into the same self-contained tree
script_dir <- here::here("xenium/label_transfer/hpc_rctd_proseg/scripts")
file.copy(list.files(script_dir, full.names = TRUE),
          file.path(OUT, "scripts"), overwrite = TRUE)

## Patch the STAGED sbatch's array upper bound to the actual chunk count.
## A manual "update --array=1-N if it differs from 80" step was documented
## here but not enforced, and the reconciliation was in fact skipped once
## (79 actual chunks vs. 80 carried over from the delivered-segmentation run)
## — a guaranteed one-task failure per submission, masked at the time by an
## unrelated fatal error (stale BASH_SOURCE path) that failed every task
## first regardless. Patching the copy under OUT means what actually gets
## rclone'd to the cluster is always correct, independent of what value is
## checked into the repo template.
staged_sbatch <- file.path(OUT, "scripts", "rctd_array.sbatch")
sbatch_lines <- readLines(staged_sbatch)
array_idx <- grep("^#SBATCH --array=", sbatch_lines)
throttle <- str_extract(sbatch_lines[array_idx], "%[0-9]+")
new_array_line <- paste0("#SBATCH --array=1-", length(chunks), if (!is.na(throttle)) throttle else "")
if (sbatch_lines[array_idx] != new_array_line) {
  message("updating ", staged_sbatch, ":\n  ", sbatch_lines[array_idx], "  ->  ", new_array_line)
  sbatch_lines[array_idx] <- new_array_line
  writeLines(sbatch_lines, staged_sbatch)
}

message("DONE in ", format(Sys.time() - t0),
        " — ", length(chunks), " chunks in ", file.path(OUT, "chunks"))
message("Total staged size: ",
        round(sum(file.size(list.files(OUT, recursive = TRUE, full.names = TRUE))) / 1e9, 2), " GB")
