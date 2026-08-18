# Proseg resegmentation on HPC — SLURM array

Reseguments all ten Xenium sections with Proseg (transcript-level, density-aware
segmentation), replacing the delivered 5 µm nucleus-expansion segmentation.

## Why: what the single-section test showed

Measured on OR52YW26_2_2 against the delivered segmentation, using RCTD (the
production label-transfer method) restricted to the paired HVC/surround ROI:

| min transcript count | | singlet frac | Glut frac | Astro frac | Glut-HVC-1 in HVC |
|---|---|---|---|---|---|
| 0   | 5 µm expansion / proseg | 41.7% / **61.0%** | 20.6% / **32.2%** | 42.5% / **28.0%** | 10.4% / **15.1%** |
| 400 | 5 µm expansion / proseg | 28.8% / **53.1%** | 31.6% / **51.6%** | 38.7% / **21.4%** | 18.5% / **29.3%** |

Proseg wins on every metric not distorted by a small-N artifact: far more clean
singlets, composition that looks like real telencephalon instead of an
astrocyte-dominated blend, and more Glut-HVC-1 cells actually detected inside
HVC. The cost is real and worth stating: proseg's median transcripts/cell is
lower than the delivered segmentation (375 vs 437 on this section), because it
correctly discards transcripts that don't belong to the seeded cell rather than
sweeping them in with a fixed-radius expansion. See
`xenium/label_transfer/xenium_label_transfer/proseg_vs_expansion_rctd_depth.csv`
for the full comparison.

## Parallelism: two levels, only one of which is ours to design

**Level 1 — across sections.** Ten independent tissue pieces with no spatial
relationship to each other. This is the embarrassingly-parallel cut the SLURM
array uses: one task per section.

**Level 2 — within a section, proseg already does this itself.** It is not a
per-FOV or whole-tissue-serial tool. It partitions space into quads
(`--quad-size`, default 150 µm) and updates non-adjacent quads concurrently — a
checkerboard/graph-coloring scheme — targeting `--cells-per-chunk` (default 100)
cells per chunk. On the 165,094-cell test section that is **~1,650 independent
chunks**, far more than the 22 cores used locally could saturate. `-t` is what
parallelizes over this internal scheme.

The consequence: **do not hand-tile a section spatially.** Splitting by FOV or
by a manual bounding-box grid would cut across the sampler's spatial
dependencies at the tile boundaries and is not how the tool is built to be
split — proseg's own chunking already does this correctly, with graph coloring
guaranteeing chunks updated in the same round don't share a boundary. The
correct lever is simply **cores per task**: give each section-task a full
node's cores (`-t` / `--cpus-per-task`) rather than inventing a coarser
decomposition, since chunk count is nowhere near saturating even a large node.

## Sizing

Local reference point: OR52YW26_2_2, 96.7M transcripts, 165,094 cells, 22
threads → ~45 min wall, ~2 GB transcript-metadata output.

`proseg_array.sbatch` requests 32 CPUs / 64 GB / 2 h per task as a starting
point — **edit `--cpus-per-task` to match your cluster's actual node size**,
since more cores should help further given the chunk-count headroom above.
Check `sacct -j <jobid> --format=Elapsed,MaxRSS,CPUTime` after the array
completes and tighten `--time`/`--mem` for a resubmission if they're generous.

## Everything lives in one self-contained directory

No direct rsync/ssh path between this workstation and the cluster (VPN
restrictions), so the whole job is staged as a single directory tree and moved
with `rclone` via whatever remote you have configured, rather than assembled
from pieces on each side:

```
hpc_proseg/
├── manifest.tsv                              # idx, section_id, transcripts_path (relative)
├── src/<section_id>/transcripts.parquet      # hardlinked, not copied — see below
├── scripts/
│   ├── proseg_array.sbatch
│   ├── 03_gather_status.R
│   └── environment.yml                       # conda/micromamba env for the cluster
├── logs/                                     # SLURM stdout/stderr land here
└── results/<section_id>/...                  # created by the array as it runs
```

`transcripts_path` in the manifest is relative to this directory, so the tree
is portable — it doesn't matter what absolute path it lands at on the cluster.
`00_make_manifest.R` hardlinks each `transcripts.parquet` into `src/` rather
than copying it (same filesystem locally, so this costs no extra disk space or
time — `du -sh` on the staged directory will still show ~17.5 GB despite the
originals also existing elsewhere, since a hardlink is a second name for the
same inode, not a second copy). `rclone` follows hardlinks as regular files, so
the copy that reaches the cluster is a real, independent file there.

## Steps

### 0. Environment (once, on the cluster)

No admin rights needed. One conda/micromamba env (`scripts/environment.yml`,
travels with the tree) covers both things the cluster side needs: a rust
toolchain to build the `proseg` binary (needs rustc >= 1.88; conda-forge's
`rust` package satisfies this), and a minimal R (readr/dplyr/stringr only) for
`scripts/03_gather_status.R`, which just checks which sections finished — not
required to run proseg itself.

```bash
micromamba create -f scripts/environment.yml
micromamba activate xenium-proseg
cargo install proseg   # fetches from crates.io — confirm the cluster can reach it
```

`cargo install` puts the `proseg` binary at `$CARGO_HOME/bin` (default
`~/.cargo/bin`) regardless of whether rustc came from conda or elsewhere;
`proseg_array.sbatch` already puts that on `PATH` and activates this env before
running, in case the binary needs any of the toolchain's shared libraries at
runtime.

If the cluster's network policy blocks crates.io specifically (unlike the
rclone remote), the alternative is building the binary locally and bundling it
into `src/` or a `bin/` directory in this same tree before syncing — ask and
I'll set that up instead.

### 1. Locally — stage the directory

```bash
cd <repo root>
Rscript xenium/preprocessing/hpc_proseg/00_make_manifest.R
# stages $XENIUM_PROC_DIR/hpc_proseg/ — manifest, hardlinked transcripts, scripts
```

### 2. Move it with rclone

```bash
D=~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_proseg
rclone sync "$D" remote:hpc_proseg          # remote: = whatever you've configured
# on the cluster:
rclone copy remote:hpc_proseg /scratch/$USER/hpc_proseg
```

### 3. Submit

Edit `--account` and `--partition` in `scripts/proseg_array.sbatch` first, and
`--cpus-per-task` to match a real node.

```bash
cd /scratch/$USER/hpc_proseg
mkdir -p logs                        # see note below — cheap and idempotent
sbatch scripts/proseg_array.sbatch .
```

`logs/` ships pre-created in the staged tree (with a `.keep` placeholder so it
survives transfer even if a sync tool drops empty directories), and
`results/<section_id>/` is created per-section by the script as it runs. The
`mkdir -p logs` immediately before `sbatch` is defense in depth: SLURM opens
`--output`/`--error` before the array's own script body executes, so nothing
inside the job can fix a missing `logs/` for its own first task — if the
directory didn't survive the rclone hop, the very first submission would fail
immediately with "couldn't open output file" and no other error. The command
above is a no-op if the directory is already there.

Idempotent: a task exits immediately if its section's
`results/<section_id>/cell-metadata.parquet` already exists, so a failed array
can be resubmitted as-is without redoing finished sections. Check progress
with:

```bash
Rscript scripts/03_gather_status.R /scratch/$USER/hpc_proseg
```

### 4. Pull back

Destination is `XENIUM_PROSEG_DIR` (`config/paths.R`) — the same place the
single-section local test (OR52YW26_2_2) already lives, so the full run joins
it rather than needing anything reorganized, and downstream R code can resolve
the path from config rather than a hardcoded string.

```bash
rclone sync /scratch/$USER/hpc_proseg remote:hpc_proseg
# locally:
rclone copy remote:hpc_proseg/results ~/hdd/rstudio/xenium/260811_brainard_adult-425g/proseg/
```

Each section directory then matches the layout already used for OR52YW26_2_2:
`cell-metadata.parquet`, `expected-counts.csv.gz` (gzipped MatrixMarket —
`Matrix::readMM(gzfile(...))`, genes in row order matching the panel; see
`xenium_label_transfer.qmd`'s proseg comparison code for the loader),
`transcript-metadata.parquet`, `cell-polygons.geojson.gz`.

## Before using the output

- **Apply the ROI polygons, not the collaborator's cell-ID selections.**
  Resegmentation assigns new cell IDs, so the ID-based selections (and the
  `1um_*` files from a prior segmentation run) will not match. Use
  `xenium/preprocessing/roi_polygons/roi_polygons.csv` and the
  `cells_in_polygon()` function from `roi_polygons.qmd` against the new
  `centroid_x`/`centroid_y`, which is exactly what the single-section test did.
- **`original_cell_id`** in `cell-metadata.parquet` maps back to the delivered
  segmentation's cell IDs when that's more convenient than polygons.
- **Composition still needs the depth filter**, more so than for the delivered
  segmentation — proseg's per-cell counts run lower (see the trade-off above).
  Apply `composition_min_counts` from `xenium_label_transfer.qmd` before
  reporting cell-type fractions; localization claims (enrichment in an ROI) are
  robust across the depth range already tested.
