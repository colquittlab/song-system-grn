#!/usr/bin/env Rscript
## Step 1 (RUN LOCALLY) — build the RCTD setup once and pre-slice it into
## per-chunk objects for a SLURM array.
##
## Why pre-slice: run.RCTD's fixed setup (gene selection, fitBulk,
## choose_sigma_c) does not scale with cell count, so it is paid once. But the
## fitted object carries the full 1.59M-cell spatialRNA (~19 GB resident). If
## every array task loaded it, 80 tasks would need ~1.5 TB. Slicing here makes
## each task load only its own ~20k cells (~1-2 GB).
##
## Why chunk at 20k: gather_results assigns rows into a sparse Matrix, which is
## O(N^2) in chunk size. Total gather cost across the dataset is therefore
## proportional to chunk size, while per-chunk scheduler/startup overhead grows
## as 1/chunk size. Measured optimum is 10-25k.
##
## Usage:  Rscript 01_prepare_chunks.R [outdir]

suppressMessages({
  library(tidyverse); library(Seurat); library(BPCells); library(qs2)
  library(Matrix); library(spacexr); library(here)
})
source(here::here("config/paths.R"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) > 0) args[1] else
  file.path(path.expand(XENIUM_PROC_DIR), "hpc_rctd")
CHUNK <- 2e4
UMI_MIN <- 20
SETUP_CORES <- 4          # few: this stage is dominated by worker startup
MAX_REF_PER_LABEL <- 1000
EXCLUDE <- c("Glut-Nido-3")
ASSAY_REF <- "RAW"        # NOT the CellBender-corrected assay — see qmd

dir.create(file.path(OUT, "chunks"), recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()

## ---- reference ----------------------------------------------------------
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
## RCTD rejects "/" in cell type names
levels(labels_rctd) <- str_replace_all(levels(labels_rctd), "/", "_")
name_map <- setNames(str_replace_all(levels(labels_rctd), "_", "/"), levels(labels_rctd))
rm(ca, o); gc()
message("reference: ", ncol(counts_ref), " cells, ", nlevels(labels_rctd), " labels")

## ---- query --------------------------------------------------------------
x <- qs_read(file.path(path.expand(XENIUM_PROC_DIR), "xenium_preprocess", "obj_xenium.qs2"))
cq_all <- LayerData(x, assay = "Xenium", layer = "counts")
cells_all <- colnames(x)
mdx <- x@meta.data %>% rownames_to_column("cell")
rm(x); gc()

## Offset sections along x so they do not overlap in coordinate space.
## Doublet-mode assignment does not use coordinates; this only keeps the stored
## geometry sane.
offs <- mdx %>% group_by(section_id) %>%
  summarize(w = max(x_centroid, na.rm = TRUE) - min(x_centroid, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_offset = lag(cumsum(w * 1.1), default = 0))
coords <- mdx %>% left_join(offs, by = "section_id") %>%
  transmute(cell, x = x_centroid + x_offset, y = y_centroid) %>% as.data.frame()
rownames(coords) <- coords$cell; coords$cell <- NULL
coords <- coords[cells_all, ]

## ---- setup: paid once ---------------------------------------------------
ref <- Reference(counts_ref, labels_rctd, colSums(counts_ref),
                 min_UMI = UMI_MIN, n_max_cells = MAX_REF_PER_LABEL)
sp_all <- SpatialRNA(coords, cq_all[ps, cells_all], colSums(cq_all[ps, cells_all]))
r <- create.RCTD(sp_all, ref, max_cores = SETUP_CORES,
                 UMI_min = UMI_MIN, counts_MIN = UMI_MIN)
r@config$RCTDmode <- "doublet"
message("create.RCTD done: ", format(Sys.time() - t0))
r <- fitBulk(r)
r <- choose_sigma_c(r)
message("setup complete: ", format(Sys.time() - t0))

## ---- pre-slice ----------------------------------------------------------
fit_cells <- colnames(r@spatialRNA@counts)
chunks <- split(fit_cells, ceiling(seq_along(fit_cells) / CHUNK))
message("slicing ", length(chunks), " chunks of ", CHUNK)

for (i in seq_along(chunks)) {
  rc <- r
  rc@spatialRNA <- restrict_puck(r@spatialRNA, chunks[[i]])
  ## originalSpatialRNA also carries all cells; drop it to keep slices small
  if (.hasSlot(rc, "originalSpatialRNA")) rc@originalSpatialRNA <- rc@spatialRNA
  qs_save(rc, file.path(OUT, "chunks", sprintf("chunk_%03d.qs2", i)))
  if (i %% 10 == 0) message("  sliced ", i, "/", length(chunks))
}

saveRDS(list(name_map = name_map, n_chunks = length(chunks),
             chunk_size = CHUNK, genes = ps,
             cells_per_chunk = lengths(chunks)),
        file.path(OUT, "manifest.rds"))

## Keep the array upper bound in rctd_array.sbatch in sync with the actual
## chunk count automatically. A manual "update --array=1-N if it differs"
## step was documented but not enforced, and was in fact forgotten once for
## the proseg variant of this pipeline (79 actual chunks vs. 80 hardcoded,
## carried over from this run) — a guaranteed one-task failure per submission,
## masked at the time by an unrelated fatal error that failed every task
## first. This script is what gets rsynced to the cluster (see README step 2),
## so patching it here keeps the deployed copy correct without relying on a
## human to reconcile it after every run.
sbatch_path <- here::here("xenium/label_transfer/hpc/rctd_array.sbatch")
sbatch_lines <- readLines(sbatch_path)
array_idx <- grep("^#SBATCH --array=", sbatch_lines)
throttle <- str_extract(sbatch_lines[array_idx], "%[0-9]+")
new_array_line <- paste0("#SBATCH --array=1-", length(chunks), if (!is.na(throttle)) throttle else "")
if (sbatch_lines[array_idx] != new_array_line) {
  message("updating ", sbatch_path, ":\n  ", sbatch_lines[array_idx], "  ->  ", new_array_line)
  sbatch_lines[array_idx] <- new_array_line
  writeLines(sbatch_lines, sbatch_path)
}

message("DONE in ", format(Sys.time() - t0),
        " — ", length(chunks), " chunks in ", file.path(OUT, "chunks"))
message("Total slice size: ",
        round(sum(file.size(list.files(file.path(OUT, "chunks"), full.names = TRUE))) / 1e9, 2), " GB")
