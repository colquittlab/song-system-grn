#!/usr/bin/env Rscript
## Step 1 (RUN LOCALLY) — RCTD label transfer over the proseg-resegmented
## Xenium data (all ten sections), using the HYBRID DIVISION LABELS from
## snrna/naming/hybrid_division_naming.qmd rather than the older
## `cluster` + recommended_label_set.csv scheme.
##
## Same pipeline as xenium/label_transfer/hpc_rctd_proseg/ — see that
## directory's 01_prepare_chunks.R and README for the sizing rationale, the
## fixed-vs-scaling cost split, the 20k chunking, and the proseg gene-order
## gotcha. All unchanged here. What differs is only the REFERENCE:
##
##   reference object : snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2
##   label column     : celltype_hybrid   (47 labels)
##
## Notes on the hybrid reference, all verified against the object:
##
##  * Cluster exclusions are ALREADY BAKED IN as NA in `celltype_hybrid`
##    (hybrid_division_naming.qmd's `excluded_clusters` = Glut-Nido-3,
##    Glut-Arco-4 — 659 cells). Cells are NOT removed from the object, so the
##    reference build must drop NA labels. That means this run is exclusion-
##    equivalent to hpc_rctd_proseg_no_arco4 (which excluded both) rather than
##    to hpc_rctd_proseg (which excluded only Glut-Nido-3). No separate
##    EXCLUDE list is needed, but one is kept below for future use.
##
##  * The hybrid scheme is NOT a pure rename of the old labels. Beyond the
##    division-based renaming (Glut-NC-*/Glut-NR-* -> Glut-DACH2-*,
##    Glut-Arco-* -> Glut-CACNA1H-*, Glut-NR-5 -> Glut-SATB2-1), it MERGES
##    Glut-NC-1 + Glut-Nido-1 into a single `Glut-DACH2-1`, and introduces
##    developmental/mixed labels (Glut-NSC, Glut-NB, Glut-Im, GABA-Im,
##    Glut-GABA). So per-label results are not directly comparable to the
##    earlier runs by name alone — use hybrid_label_lookup.csv to map.
##
##  * Song-region labels are anatomical rather than positional, which is the
##    main reason to prefer this reference for interpretation:
##      Glut-HVC-1  -> Glut-DACH2-HVCra    (HVC -> RA projecting)
##      Glut-HVC-1a -> Glut-DACH2-HVCra-Int
##      Glut-HVC-2  -> Glut-DACH2-HVCx     (HVC -> Area X)
##      Glut-RA     -> Glut-CACNA1H-RA
##      Glut-LMAN-1 -> Glut-DACH2-LMANsh   (shell)
##      Glut-LMAN-2 -> Glut-DACH2-LMANco   (core)
##    The shell/core assignment independently matches the spatial core/surround
##    structure seen in the LMAN zoom figures.
##
##  * The RAW assay is present in this object (checked) — required, since
##    CellBender-corrected counts break the transfer (see
##    xenium/label_transfer/xenium_label_transfer.qmd).
##
## Usage:  Rscript 01_prepare_chunks.R [outdir]

suppressMessages({
  library(tidyverse); library(Seurat); library(BPCells); library(qs2)
  library(Matrix); library(spacexr); library(here)
})
source(here::here("config/paths.R"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) > 0) args[1] else
  file.path(path.expand(XENIUM_PROC_DIR), "hpc_rctd_proseg_hybrid")
CHUNK <- 2e4
UMI_MIN <- 20
SETUP_CORES <- 4
MAX_REF_PER_LABEL <- 1000
LABEL_COL <- "celltype_hybrid"
## Extra exclusions on top of the NA labels already in celltype_hybrid.
## Empty: Glut-Nido-3 and Glut-Arco-4 are already NA there.
EXCLUDE <- character(0)
ASSAY_REF <- "RAW"

REF_OBJ <- here::here("snrna/naming/hybrid_division_naming", "obj_hybrid_labels.qs2")

dir.create(file.path(OUT, "chunks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
file.create(file.path(OUT, "logs", ".keep"))   # see hpc_proseg/README on why
t0 <- Sys.time()

sections <- c("OR52YW26_1_4", "OR52YW26_1_7", "OR52YW26_2_2", "OR52YW26_2_4", "OR52YW26_2_7",
             "OR69PU4_1_4", "OR69PU4_1_7", "OR69PU4_2_2", "OR69PU4_2_4", "OR69PU4_2_7")

## ---- reference: hybrid division labels ----------------------------------
panel <- read_csv(XENIUM_PANEL_CSV, show_col_types = FALSE)$Gene

o <- qs_read(REF_OBJ)
stopifnot(LABEL_COL %in% colnames(o@meta.data))
stopifnot(ASSAY_REF %in% names(o@assays))
o <- subset(o, cells = Cells(o)[!grepl("/", o$assignment)])
ca <- LayerData(o, assay = ASSAY_REF, layer = "counts")
ps <- intersect(panel, rownames(ca))

md <- o@meta.data %>% rownames_to_column("cell") %>%
  mutate(hybrid_label = as.character(.data[[LABEL_COL]])) %>%
  filter(!is.na(hybrid_label), !hybrid_label %in% EXCLUDE)
set.seed(10)
md <- md %>% group_by(hybrid_label) %>% slice_sample(n = MAX_REF_PER_LABEL) %>% ungroup()

counts_ref <- as(ca[ps, md$cell], "dgCMatrix")
labels_rctd <- setNames(factor(md$hybrid_label), md$cell)
## RCTD rejects "/" in cell type names — none of the hybrid labels contain one,
## but substitute defensively and map back, same as the other variants.
levels(labels_rctd) <- str_replace_all(levels(labels_rctd), "/", "_")
name_map <- setNames(str_replace_all(levels(labels_rctd), "_", "/"), levels(labels_rctd))
## NOTE: hybrid labels legitimately contain "_"?  They do not (all use "-"),
## so the reverse map above is identity in practice. Guard against a future
## label with "_" silently becoming "/":
if (any(grepl("_", levels(labels_rctd)))) {
  name_map <- setNames(levels(labels_rctd), levels(labels_rctd))
  message("NOTE: labels contain '_', name_map left as identity")
}
rm(ca, o); gc()
message("reference: ", ncol(counts_ref), " cells, ", nlevels(labels_rctd), " labels, ",
        length(ps), " genes  [", LABEL_COL, "]")

## ---- query: load each proseg section, re-index genes by NAME ------------
gene_list_dir <- here::here("xenium/label_transfer/hpc_rctd_proseg_hybrid/gene_lists")

load_section <- function(sid) {
  d <- file.path(path.expand(XENIUM_PROSEG_DIR), sid)
  genes_sid <- readLines(file.path(gene_list_dir, paste0(sid, ".txt")))
  genes_sid <- genes_sid[nzchar(genes_sid)]
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
             label_col = LABEL_COL, ref_obj = REF_OBJ,
             cells_per_chunk = lengths(chunks)),
        file.path(OUT, "manifest.rds"))

## copy the SLURM/array scripts + env spec into the same self-contained tree
script_dir <- here::here("xenium/label_transfer/hpc_rctd_proseg_hybrid/scripts")
file.copy(list.files(script_dir, full.names = TRUE),
          file.path(OUT, "scripts"), overwrite = TRUE)

## Patch the STAGED sbatch's array upper bound to the actual chunk count -- see
## hpc_rctd_proseg/01_prepare_chunks.R for why this is automated rather than a
## documented manual step.
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
