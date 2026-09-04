# Explore + export the ADULT chicken Seurat object (Zaremba et al. 2025) to HDF5 pieces
# Python can assemble into an h5ad. Same dgCMatrix-via-hdf5r approach as the dev chicken
# export (SeuratDisk/zellkonverter not installed in r44).
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(hdf5r)
})

RDS <- "/private/home/bcolquit/group/sc_datasets/zaremba_2025/Gg_adult_snRNA_seq_srt.rds/Gg_adult_snRNA_seq_srt.rds"
OUT <- "/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-gg-adult/data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("=== loading", RDS, "===\n"); flush.console()
t0 <- Sys.time()
srt <- readRDS(RDS)
cat("loaded in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n\n")

cat("=== object ===\n"); print(srt)
cat("\nassays:", paste(Assays(srt), collapse = ", "), "\n")
cat("default assay:", DefaultAssay(srt), "\n")
cat("reductions:", paste(Reductions(srt), collapse = ", "), "\n\n")

cat("=== meta.data columns ===\n")
md <- srt@meta.data
print(colnames(md))
cat("\ndims:", nrow(md), "cells\n\n")

cat("=== candidate label columns ===\n")
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (is.factor(v) || is.character(v)) {
    u <- length(unique(v))
    if (u > 1 && u <= 800) cat(sprintf("  %-30s %4d levels\n", cn, u))
  }
}
cat("\n")
for (a in Assays(srt)) {
  cat("assay", a, "layers:", paste(SeuratObject::Layers(srt[[a]]), collapse = ", "), "\n")
}
cat("\n=== gene name sample (rownames) ===\n")
print(head(rownames(srt), 20))
cat("\n=== a candidate cluster column, if present ===\n")
for (cn in c("cluster","seurat_clusters","cell_type","celltype","annotation","region","position","dissection")) {
  if (cn %in% colnames(md)) { cat(cn, ":\n"); print(table(md[[cn]])[1:min(20,length(table(md[[cn]])))]) }
}

# --- export raw counts -------------------------------------------------------
assay <- if ("RNA" %in% Assays(srt)) "RNA" else DefaultAssay(srt)
cat("\nexporting counts from assay:", assay, "\n")
m <- SeuratObject::LayerData(srt[[assay]], layer = "counts")
if (is.null(m) || nrow(m) == 0) m <- SeuratObject::LayerData(srt[[assay]], layer = "data")
m <- as(m, "CsparseMatrix")

cat("matrix:", nrow(m), "genes x", ncol(m), "cells; nnz =", length(m@x), "\n")
cat("integral counts?", all(m@x == round(m@x)), " max =", max(m@x), "\n")

h5 <- H5File$new(file.path(OUT, "gg_adult_counts.h5"), mode = "w")
h5[["data"]]    <- m@x
h5[["indices"]] <- m@i
h5[["indptr"]]  <- m@p
h5[["shape"]]   <- c(nrow(m), ncol(m))
h5[["genes"]]   <- rownames(m)
h5[["cells"]]   <- colnames(m)
h5$close_all()

md$CellID <- rownames(md)
write.csv(md, file.path(OUT, "gg_adult_metadata.csv"), row.names = FALSE)

cat("\nwrote", file.path(OUT, "gg_adult_counts.h5"), "and gg_adult_metadata.csv\n")
cat("done at", format(Sys.time()), "\n")
