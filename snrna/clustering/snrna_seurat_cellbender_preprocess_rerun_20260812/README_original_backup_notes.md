# Backup of snrna_seurat_cellbender_preprocess outputs, 2026-08-12

Taken before re-running `snrna_seurat_cellbender_preprocess.qmd` with genotype-derived
individual assignments (branch `snrna-reclustering`). Every `redo` flag in that notebook is
`TRUE`, so it overwrites its outputs in place.

## What is here

Copies of the legacy-`qs` outputs. **These are `qs` format, not `qs2` -- read them with
`qs::qread()`, not `qs2::qs_read()`**, which errors with "qs-legacy format detected":

| file | size | note |
|---|---|---|
| `obj.qs` | 428 MB | merged, pre-QC-filter object |
| `obj_processed.qs` | 1.3 GB | after doublet removal, SCT, clustering |
| `md_doubletfinder.rds` | 200 KB | DoubletFinder classifications (plain RDS, unaffected) |
| `bpcells/` | 2.5 GB | BPCells matrix dirs, rewritten with `overwrite = TRUE` |
| `umap_*.png` | 9 MB | |

Symlinks, copied as symlinks (`cp -P`) and still pointing at the untouched originals:

- `obj_clustered.qs` -> `/hdd/.../snrna_seurat_cellbender.0.01_preprocess/obj_clustered.qs`
  (3.26 GB, 2025-12-11)
- `obj_clustered.qs2` -> `/hdd/.../obj_clustered.qs2` (3.12 GB, 2026-05-11)

The re-run no longer writes any of the `.qs` names: the notebook was migrated to qs2 and its
outputs are now `obj.qs2`, `obj_processed.qs2`, `obj_clustered.qs2` and the `obj_*.qs2`
subset objects. So the local `.qs` copies above are belt-and-braces rather than strictly
necessary -- except `bpcells/`, which is still rewritten in place.

## The change made outside this directory

`obj_clustered.qs2` is now the notebook's write target (`qs_save(obj_filt, file = out_fname)`),
and it was a symlink to a 3.12 GB file on `/hdd` -- writing through it would have truncated
that file, which copying the symlink alone would not have prevented. So it was re-pointed at a
new, dated target:

    snrna_seurat_cellbender_preprocess/obj_clustered.qs2
      -> /hdd/.../snrna_seurat_cellbender.0.01_preprocess/obj_clustered_20260812.qs2

The new run therefore lands on `/hdd` as before (keeping ~3 GB off `/ssd`, which is at 93%),
and the May object is untouched. `obj_clustered.qs` was left pointing at its original December
target: the notebook no longer writes that path, and several scripts still read it
(`snrna/reduction_viz/combined_all_umap*.R`, `snrna/deg/all_deg.qmd`,
`snrna/trees/celltypes_hclust_all.qmd`, `snrna/integration/...`, `grn/analysis/*.qmd`).

**Consequence worth knowing:** after the re-run, `obj_clustered.qs` is the *old* December
object and `obj_clustered.qs2` is the new one. Scripts reading the `.qs` name will silently get
stale data. The `.qs2` readers (`snrna/trees/celltypes_hclust_glut*.qmd`, the xenium
label-transfer notebooks) pick up the new object automatically.

## Restoring

    cd snrna/clustering/snrna_seurat_cellbender_preprocess
    BK=../snrna_seurat_cellbender_preprocess_backup_20260812
    cp -a $BK/{obj.qs,obj_processed.qs,md_doubletfinder.rds} .   # qs format; qs::qread to load
    rm -rf bpcells && cp -a $BK/bpcells .
    ln -sfn /hdd/brad/rstudio/snRNA/snrna_cellranger/snrna_seurat_cellbender.0.01_preprocess/obj_clustered.qs2 obj_clustered.qs2
