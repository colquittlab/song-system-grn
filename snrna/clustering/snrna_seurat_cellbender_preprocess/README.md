# snRNA preprocessing outputs — the main object

This is the **main** snRNA object directory. Everything downstream reads from here.

`snrna_seurat_cellbender_preprocess.qmd` is at its published state (tag `v1.0`) and writes to
`SNRNA_SEURAT_DIR` (`~/hdd/rstudio/snRNA/snrna_cellranger/snrna_seurat_cellbender.0.01_preprocess/`,
see `config/paths.R`), in legacy `qs` format. This directory is the repo-local view of that:

| entry | what it is |
|---|---|
| `obj_clustered.qs` | symlink → `SNRNA_SEURAT_DIR/obj_clustered.qs` (3.26 GB, 2025-12-11) |
| `obj_clustered.qs2` | symlink → `SNRNA_SEURAT_DIR/obj_clustered.qs2` (3.12 GB, 2026-05-11) — the object the naming and analysis notebooks read |
| `obj.qs`, `obj_processed.qs` | local copies of the merged / doublet-filtered intermediates (`qs::qread`, **not** `qs2::qs_read`) |
| `md_doubletfinder.rds` | DoubletFinder classifications |
| `bpcells/` | on-disk BPCells matrices referenced by the `obj*.qs` objects — the paths are baked in, so this directory must keep this name |
| `obj_clustered.h5ad` | AnnData export, written by `../convert_to_h5ad_snrna.R` |
| `umap_*.png` | QC/marker UMAPs from the notebook |

## Format gotcha

`obj.qs` / `obj_processed.qs` / `obj_clustered.qs` are legacy `qs`. Read them with `qs::qread()`;
`qs2::qs_read()` fails on them with "qs-legacy format detected". Only `obj_clustered.qs2` is qs2.

## History

The 2026-08-12 souporcell/genotype re-run of the preprocessing notebook (branch
`snrna-reclustering`) previously occupied this path. It has been demoted to
`../snrna_seurat_cellbender_preprocess_rerun_20260812/`, which also holds the version of the
notebook that produced it. This directory is exactly the pre-re-run state that was preserved as
`..._backup_20260812`, restored to its proper name.
