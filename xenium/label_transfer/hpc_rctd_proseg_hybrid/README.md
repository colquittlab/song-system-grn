# RCTD over proseg data — hybrid division labels

Same pipeline as `xenium/label_transfer/hpc_rctd_proseg/` (proseg-resegmented
Xenium, all ten sections, RAW-assay reference, 20k chunks, `min_UMI = 20`,
1000 reference cells/label). See that README for the sizing rationale and the
segmentation/reference findings this run inherits.

**What differs: the reference labels.**

| | `hpc_rctd_proseg` | this run |
|---|---|---|
| reference object | `snrna/clustering/.../obj_clustered.qs2` | `snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2` |
| label source | `cluster` + `panel_resolution_benchmark/recommended_label_set.csv` | `celltype_hybrid` |
| labels | 48 | 47 |

`celltype_hybrid` comes from `snrna/naming/hybrid_division_naming.qmd`, the
production naming notebook. `hybrid_label_lookup.csv` in that output directory
maps old cluster names to hybrid labels.

## Three things about this reference worth knowing

**1. Exclusions are already baked in.** `hybrid_division_naming.qmd`'s
`excluded_clusters` (`Glut-Nido-3`, `Glut-Arco-4`) are set to `NA` in
`celltype_hybrid` — 659 cells. The cells are *not* removed from the object, so
`01_prepare_chunks.R` drops `NA` labels rather than naming exclusions itself.
This makes the run **exclusion-equivalent to `hpc_rctd_proseg_no_arco4`**, not
to `hpc_rctd_proseg`. Both of those exclusions were independently supported by
the earlier RCTD work (Glut-Nido-3 as a sink; Glut-Arco-4 redistributing
cleanly to neighbouring Arco subtypes).

**2. It is not a pure rename.** Beyond the division-based renaming, it merges
`Glut-NC-1` + `Glut-Nido-1` into one `Glut-DACH2-1` (2,561 cells), and adds
developmental/mixed labels (`Glut-NSC`, `Glut-NB`, `Glut-Im`, `GABA-Im`,
`Glut-GABA`). Per-label comparison against the earlier runs therefore needs
`hybrid_label_lookup.csv` — matching on label name alone will mislead.

Note the merge differs from the one the panel benchmark chose: the benchmark
merged `Glut-Nido-1` with `Glut-NC-2` (having been asked to protect
`Glut-NC-1`), whereas the hybrid scheme merges it with `Glut-NC-1`.

**3. Song-region labels are anatomical, not positional** — the main reason to
prefer this reference for interpretation:

| old | hybrid |
|---|---|
| `Glut-HVC-1` | `Glut-DACH2-HVCra` (HVC → RA projecting) |
| `Glut-HVC-1a` | `Glut-DACH2-HVCra-Pre` |
| `Glut-HVC-2` | `Glut-DACH2-HVCx` (HVC → Area X) |
| `Glut-RA` | `Glut-CACNA1H-RA` |
| `Glut-LMAN-1` | `Glut-DACH2-LMANsh` (shell) |
| `Glut-LMAN-2` | `Glut-DACH2-LMANco` (core) |

The shell/core assignment independently agrees with the spatial structure in
the LMAN zoom figures (LMAN-2 dense inner core, LMAN-1 broader surround).

## Steps

Identical to the other variants, pointed here:

```bash
cd <repo root>
Rscript xenium/label_transfer/hpc_rctd_proseg_hybrid/01_prepare_chunks.R
# stages $XENIUM_PROC_DIR/hpc_rctd_proseg_hybrid/

D=~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid
rclone sync "$D" remote:hpc_rctd_proseg_hybrid
# on the cluster:
rclone copy remote:hpc_rctd_proseg_hybrid /scratch/$USER/hpc_rctd_proseg_hybrid
cd /scratch/$USER/hpc_rctd_proseg_hybrid
micromamba activate xenium-rctd   # reuse the existing env
mkdir -p logs
sbatch scripts/rctd_array.sbatch .
```

`--partition=medium` and no `--account` are already set, and `--array` is
patched to the real chunk count by step 1 — nothing to edit before submitting
on phoenix.

Gather and pull back:

```bash
Rscript scripts/03_gather.R .
rclone sync /scratch/$USER/hpc_rctd_proseg_hybrid remote:hpc_rctd_proseg_hybrid
# locally:
rclone copy remote:hpc_rctd_proseg_hybrid/rctd_all.csv.gz \
  ~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid/
```

## After it lands

- `proseg_ncount.csv.gz` is **not** produced by the array — extract it from the
  chunks the same way as the other runs before any depth-filtered
  (`nCount >= 400`) plotting, or the spatial-map scripts will fail on a missing
  file.
- Spatial maps: copy `proseg_spatial_maps.qmd` to a `_hybrid` variant, point
  `hpc_dir`/`out_dir` at this run, and **update `song_types`** to the hybrid
  names (the old `Glut-HVC-1` etc. no longer exist, so the existing
  `stopifnot(all(song_types %in% ...))` will fail loudly rather than silently
  plot nothing).
- Cell naming (`<section_id>_proseg_<cell>`) is unchanged, so ROI polygons,
  the depth filter, and the comparison scripts all apply as-is.
