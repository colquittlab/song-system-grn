# Demoted: the 2026-08-12 re-run outputs

**Not the main object.** The main object is `../snrna_seurat_cellbender_preprocess/`.

These are the outputs of the 2026-08-12 re-run of `snrna_seurat_cellbender_preprocess.qmd` on
branch `snrna-reclustering` — the version migrated to `qs2`, writing into the repo tree rather than
to `SNRNA_SEURAT_DIR`, and using souporcell genotype-derived `individual`/`batch` assignments in
place of the per-library `assignment`. That version of the notebook is preserved here as
`snrna_seurat_cellbender_preprocess_rerun_20260812.qmd`; the notebook in the repo has been returned
to its published `v1.0` state.

Nothing in the analysis chain reads this directory. Kept so the re-run is not lost.

## Caveats if you ever come back to these

- The `obj*.qs2` objects reference BPCells matrices by **absolute path**, baked in as
  `.../snrna_seurat_cellbender_preprocess/bpcells/...`. That path now resolves to the *main*
  object's matrices, not the `bpcells/` copy in this directory. Loading these objects and touching
  the counts will read the wrong matrices — rename this directory back to
  `snrna_seurat_cellbender_preprocess` (moving the main one aside) before using them.
- `obj_clustered.qs2` here is a symlink to
  `SNRNA_SEURAT_DIR/obj_clustered_20260812.qs2` (~3 GB, the re-run's clustered object, on `/hdd`).
- `obj_clustered.qs` is a symlink to the *December 2025* object — i.e. it is **not** from this
  re-run; it was never rewritten because the re-run only wrote `qs2` names.
- `README_original_backup_notes.md` is the note written when this state was first backed up on
  2026-08-12, kept verbatim for the record. Its "Restoring" section is already carried out.
