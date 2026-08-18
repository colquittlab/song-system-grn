# RCTD over proseg-resegmented Xenium data — SLURM array

Label transfer for all ten Xenium sections **after** proseg resegmentation
(`xenium/preprocessing/hpc_proseg/`), following the same design as
`xenium/label_transfer/hpc/` (which ran this over the delivered 5 µm
nucleus-expansion segmentation). Read that directory's README first — the
sizing rationale, the `choose_sigma_c`/`fitPixels` cost split, and the
20k-cell chunking are unchanged and not re-derived here.

## What's different from the delivered-segmentation run

**The query is built from proseg's per-section output**
(`expected-counts.csv.gz` + `cell-metadata.parquet` in `XENIUM_PROSEG_DIR`)
instead of the merged `obj_xenium.qs2`, and combined across all ten sections in
`01_prepare_chunks.R` rather than already being one object.

**proseg's gene order is not consistent across sections — confirmed, not
assumed.** Checked directly:

```
OR52YW26_1_4 starts: GAD1, SLC1A3, SNAP25, PLP1, B2M
OR69PU4_2_7  starts: SLC1A3, SNAP25, RGS10, GAD1, PVALB
```

Row position in `expected-counts.csv.gz` (a plain MatrixMarket matrix — no row
names) is only meaningful together with that section's own gene list, which
proseg only writes into its zarr store (`tables/table/var/gene`), not as plain
text. `gene_lists/<section_id>.txt` was extracted once (via a Python
environment with `zarr` — `~/micromamba/envs/regvelo` locally had it; any env
with `zarr` works) and is read by `01_prepare_chunks.R` to re-index every
section's matrix **by gene name** before combining. Skipping this would
silently scramble expression profiles in nine of the ten sections without
erroring — regenerate `gene_lists/` before rerunning `01_prepare_chunks.R` if
proseg is ever re-run.

**Cell IDs restart at 0 in every proseg section.** Renamed to
`<section_id>_proseg_<cell>` to stay globally unique in the combined object.
`cell-metadata.parquet`'s `original_cell_id` still maps back to the delivered
segmentation's IDs if that's ever needed.

**Reference config is unchanged**: `RAW` assay (not CellBender-corrected —
see `xenium/label_transfer/xenium_label_transfer.qmd` for why), Glut-Nido-3
excluded (confirmed sink, see `panel_resolution_benchmark`'s merge log and the
label-transfer script's exclusion rationale), 1000 cells/label max, `min_UMI = 20`.

## Steps

### 0. Environment (once, on the cluster)

Same env as the delivered-segmentation run — reuse it if that cluster account
already has it from the earlier run, no need to recreate.

```bash
micromamba create -f scripts/environment.yml
micromamba activate xenium-rctd
```

(No proseg/rust toolchain needed here — this is pure R/spacexr over already-
segmented data.)

### 1. Locally — build and stage

```bash
cd <repo root>
Rscript xenium/label_transfer/hpc_rctd_proseg/01_prepare_chunks.R
# stages $XENIUM_PROC_DIR/hpc_rctd_proseg/ — manifest.rds, chunks/, scripts/, logs/
```

Takes ~15-25 min: loading + gene re-indexing for ten sections, then the same
`create.RCTD`/`fitBulk`/`choose_sigma_c` setup as before (paid once, ~9-17 min
at 4 cores per the delivered-segmentation run's measurement).

The reported chunk count (proseg's ~1.59M cells at 20k/chunk is ~80, close to
but not guaranteed identical to the delivered-segmentation run) no longer
needs manual reconciliation: this script patches the staged
`scripts/rctd_array.sbatch`'s `--array` upper bound to the actual chunk count
automatically after slicing. (This step was previously manual and was in fact
missed once — 79 actual chunks vs. 80 hardcoded — a guaranteed one-task
failure per submission; see `02_fit_chunk.R`'s own guard for the same failure
mode's second line of defense, in case a stale range ever gets shipped by some
other path.)

### 2. Move it with rclone

Same self-contained-directory reasoning as `hpc_proseg/` (VPN blocks a direct
rsync/ssh path to this cluster): everything needed —
`manifest.rds`, `chunks/`, `scripts/`, `logs/` — lives under one tree, moved as
a whole.

```bash
D=~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg
rclone sync "$D" remote:hpc_rctd_proseg
# on the cluster:
rclone copy remote:hpc_rctd_proseg /scratch/$USER/hpc_rctd_proseg
```

### 3. Submit

`--partition=medium` and no `--account` are already set (confirmed working on
phoenix/prism), and the `--array` count is already correct — patched
automatically by step 1. Nothing to edit before submitting on this cluster.

```bash
cd /scratch/$USER/hpc_rctd_proseg
mkdir -p logs   # defensive — see the note in rctd_array.sbatch on why
sbatch scripts/rctd_array.sbatch .
```

Idempotent per chunk. Resume a partial run:

```bash
sbatch --array=$(Rscript scripts/04_missing.R .) scripts/rctd_array.sbatch .
```

### 4. Gather and pull back

```bash
Rscript scripts/03_gather.R .   # writes rctd_all.csv.gz, reports any missing chunks
rclone sync /scratch/$USER/hpc_rctd_proseg remote:hpc_rctd_proseg
# locally:
rclone copy remote:hpc_rctd_proseg/rctd_all.csv.gz ~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg/
```

## After it lands

- **Cell metadata for ROI/depth analysis** is in `proseg_cell_metadata.csv.gz`
  (written by `01_prepare_chunks.R`, includes `x_centroid`/`y_centroid` in
  original per-section coordinates, not the offset combined-object coordinates
  used only for the RCTD run itself) — join on `cell` against `rctd_all.csv.gz`.
- **Apply the ROI polygons**, not cell-ID selections — proseg cells have new
  IDs (see `xenium/preprocessing/roi_polygons.qmd` and
  `xenium/preprocessing/hpc_proseg/README.md`'s "Before using the output").
- **Composition still needs a depth filter** — proseg's per-cell counts run
  lower than the delivered segmentation (measured: median 375 vs 437 on
  OR52YW26_2_2). Localization (ROI enrichment) held up across the depth range
  already tested there; composition claims should use
  `composition_min_counts` from `xenium_label_transfer.qmd`, likely re-tuned
  for proseg's lower depth rather than reused as-is.
- **The real comparison**: re-run the ROI enrichment test from
  `xenium_label_transfer.qmd` (Glut-HVC-1 in HVC vs. surround, all five paired
  sections) on this output and set it beside the delivered-segmentation result
  already in `hpc_roi_enrichment.csv` / `hpc_roi_enrichment_pooled.csv` — that
  single-section test showed proseg improving singlet rate and composition at
  the cost of some raw depth; this is the full-dataset version of that check.
