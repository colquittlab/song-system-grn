## Inputs for the RPCA re-run on OrthoFinder orthologs (plusSATB2_noPre finch/chicken subset).
##
## The original RPCA run (rpca/) used reciprocal-best-BLASTP-hit orthologs (9,492 genes);
## GSI used OrthoFinder 2.5.4 one-to-one orthologs. This re-run keeps everything else
## fixed -- the same cells (the median-subsampled lists of the original run) and the same
## labels -- and changes only the gene set, to GSI's OrthoFinder table
## (gsi/data/orthologs_used_plusSATB2_noPre.csv; columns finch / mouse = finch / chicken
## symbol, already restricted to genes present in both datasets).
##
## Counts:
##   finch    RNA assay counts (CellBender-corrected UMIs) of the hybrid-labelled object.
##            This is what the cluster's bf_adult_hybrid.h5ad held: on GSI's subsampled
##            plusSATB2_noPre cells, RNA counts reproduce GSI's SCTransform gene list
##            exactly (8,219/8,219 genes detected in >= 5 cells); RAW counts do not.
##   chicken  RNA counts of Zaremba et al. 2025's adult object, rounded to integers and
##            translated from Ensembl 87 gene IDs to gene names exactly as
##            datasets/snrna-bf-adult_snrna-gg-adult/build_gg_adult_h5ad.py does (IDs
##            without a name keep the ID; duplicated names made unique). The ID -> name
##            map (data/ens87_gg_geneid_to_name.tsv) was rebuilt from Ensembl's release-87
##            Gallus_gallus-5.0 GTF, because the cluster's copy is not available here.
##
## Output: data/seurat/<TAG>_{query,ref}_{matrix.mtx,genes.tsv,cells.tsv,meta.csv}, the
## layout export_for_seurat.py writes and finch-integration-toolkit/R/cca_transfer.R reads,
## with both species' genes named by the finch symbol of each ortholog pair.
##
##   Rscript --no-init-file snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/prepare_inputs.R
renv_lib <- Sys.getenv("RENV_LIB", "/ssd/brad/rstudio/multiome/song-system-grn/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))
suppressMessages({library(Seurat); library(qs2); library(Matrix)})

here_dir <- "snrna/integration/gg_glutonly_hybrid/rpca_orthofinder"
gg_dir <- "snrna/integration/gg_glutonly_hybrid"
old <- file.path(gg_dir, "rpca/data/seurat/gg_glutonly_hybrid_plusSATB2_noPre")
tag <- "gg_glutonly_hybrid_plusSATB2_noPre_orthofinder"
out <- file.path(here_dir, "data/seurat", tag)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

FINCH_QS2 <- Sys.getenv("FINCH_QS2", "/ssd/brad/rstudio/multiome/song-system-grn/snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2")
CHICKEN_RDS <- Sys.getenv("CHICKEN_RDS", "/hdd/sc_datasets/zaremba_gg_dev_2024/Gg_adult_snRNA_seq_srt.rds")

orth <- read.csv(file.path(gg_dir, "gsi/data/orthologs_used_plusSATB2_noPre.csv"))
stopifnot("ortholog table is not 1:1" = !anyDuplicated(orth$finch) && !anyDuplicated(orth$mouse))

write_side <- function(m, labels, side) {
  stopifnot(identical(colnames(m), names(labels)))
  writeMM(m, paste0(out, "_", side, "_matrix.mtx"))
  writeLines(rownames(m), paste0(out, "_", side, "_genes.tsv"))
  writeLines(colnames(m), paste0(out, "_", side, "_cells.tsv"))
  write.csv(data.frame(label = unname(labels), row.names = colnames(m)), paste0(out, "_", side, "_meta.csv"))
  cat(sprintf("%s: %d genes x %d cells, %d labels\n", side, nrow(m), ncol(m), length(unique(labels))))
}

## SCTransform keeps genes detected in >= 5 cells; GSI's SCT gene list (same cells, before UMI
## thinning only removes counts) therefore checks both the counts and the gene mapping.
check_against_gsi_sct <- function(m, gsi_prefix, genes_in_gsi_order) {
  gsi_cells <- readLines(paste0(gsi_prefix, "_cells.tsv"))
  sct_genes <- readLines(paste0(gsi_prefix, "_sct_genes.tsv"))
  stopifnot("GSI cells differ from the RPCA cells" = setequal(gsi_cells, colnames(m)))
  detected <- genes_in_gsi_order[Matrix::rowSums(m > 0) >= 5]
  cat(sprintf("  check vs GSI SCT genes: %d detected in >= 5 cells, %d kept by SCT, %d kept-but-undetected\n",
              length(detected), length(sct_genes), sum(!sct_genes %in% detected)))
  stopifnot("counts/gene mapping disagree with GSI's SCT gene list" = all(sct_genes %in% detected))
}

## ---- finch (query) ---------------------------------------------------------------------
q_cells <- readLines(paste0(old, "_query_cells.tsv"))
q_meta <- read.csv(paste0(old, "_query_meta.csv"), row.names = 1)
stopifnot(identical(rownames(q_meta), q_cells))
o <- qs_read(FINCH_QS2, nthreads = 8)
stopifnot(all(q_cells %in% colnames(o)), all(orth$finch %in% rownames(o[["RNA"]])))
lab_obj <- sub("^Glut-DACH2-HVCra-Pre$", "Glut-DACH2-HVCra-Int", as.character(o$celltype_hybrid[q_cells]))
stopifnot("labels differ from the object" = identical(lab_obj, q_meta$label))
fm <- as(LayerData(o, assay = "RNA", layer = "counts")[orth$finch, q_cells], "dgCMatrix")
stopifnot(all(fm@x == round(fm@x)))
rm(o); invisible(gc())
check_against_gsi_sct(fm, file.path(gg_dir, "gsi/data/finch_plusSATB2_noPre"), orth$finch)
write_side(fm, setNames(q_meta$label, q_cells), "query")

## ---- chicken (reference) -------------------------------------------------------------------
r_cells <- readLines(paste0(old, "_ref_cells.tsv"))
r_meta <- read.csv(paste0(old, "_ref_meta.csv"), row.names = 1)
stopifnot(identical(rownames(r_meta), r_cells))
srt <- readRDS(CHICKEN_RDS)
cm <- LayerData(srt[["RNA"]], layer = "counts")
stopifnot(all(r_cells %in% colnames(cm)))
cm <- as(cm[, r_cells], "CsparseMatrix")
rm(srt); invisible(gc())
map <- read.delim(file.path(here_dir, "data/ens87_gg_geneid_to_name.tsv"), header = FALSE, col.names = c("id", "name"))
keep <- rownames(cm) %in% map$id
cm <- cm[keep, ]
nm <- map$name[match(rownames(cm), map$id)]
## anndata var_names_make_unique: 2nd, 3rd ... occurrence get -1, -2 ...
dup_n <- ave(seq_along(nm), nm, FUN = seq_along) - 1
nm <- ifelse(dup_n > 0, paste0(nm, "-", dup_n), nm)
rownames(cm) <- nm
miss <- setdiff(orth$mouse, rownames(cm))
stopifnot("OrthoFinder chicken genes missing after ID translation" = length(miss) == 0)
cm <- cm[orth$mouse, ]
cm@x <- round(cm@x)
cm <- drop0(cm)
check_against_gsi_sct(cm, file.path(gg_dir, "gsi/data/chicken_plusSATB2_noPre"), orth$mouse)
rownames(cm) <- orth$finch   # common (finch) symbol space, as export_for_seurat.py does
write_side(cm, setNames(r_meta$label, r_cells), "ref")
