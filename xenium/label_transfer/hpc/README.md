# RCTD on HPC — SLURM array

Runs RCTD label transfer for the ten Xenium sections as a SLURM job array.

## Why it is split this way

RCTD's cost splits into two very different parts, measured on this dataset:

| Stage | Scales with cell count? | Measured cost | Cores that help |
|---|---|---|---|
| `create.RCTD` (gene selection) | no | 1.1 min | few |
| `fitBulk` | no | 0.7 s | few |
| `choose_sigma_c` | **no** — fits `N_fit = 1000` cells for ≤8 epochs | 15.5 min at 4 cores | **fewer is faster** |
| `fitPixels` | **yes** | ~17 min per 20k cells | array width |
| `gather_results` | O(N²) *in chunk size* | ~4 min per 20k chunk | none (serial) |

Two counter-intuitive measurements drive the design:

1. **`choose_sigma_c` is ~6× faster on 4 cores than on 44.** It decomposes only
   1,000 cells per epoch, but spacexr spawns a fresh PSOCK cluster each epoch and
   every worker re-loads the R library stack. Startup swamps the work.
2. **spacexr only keeps ~2.4 cores busy** inside one process, whatever
   `max_cores` says (measured: 242% CPU with `max_cores = 22`). So parallelism
   must come from **array width**, not cores per task.

Hence: pay the fixed setup once locally, slice the query, and fan the scaling
stage out over independent array tasks with 4 cores each.

**Chunk size 20k** is a measured optimum. `gather_results` writes rows into a
sparse `Matrix`, which is O(N²) in chunk size, so total gather cost across the
dataset is proportional to chunk size — while per-task overhead grows as
1/chunk size. At 200k chunks the run projected to ~28 h; at 20k it is ~21 min
per chunk.

## Expected wall time

~80 chunks × ~21 min each. Serial that is ~28 h; the array collapses it:

| Concurrent tasks | Wall time |
|---|---|
| 20 | ~1.5 h |
| 40 (default `%40`) | ~45 min |
| 80 | ~25 min |

Each task needs 4 CPUs and peaks near 6 GB, so 80 concurrent tasks is ~320 CPUs
and ~500 GB across the cluster — modest for an array.

## Steps

### 1. Locally — build and slice (~25 min)

```bash
cd <repo root>
Rscript xenium/label_transfer/hpc/01_prepare_chunks.R
# writes $XENIUM_PROC_DIR/hpc_rctd/{chunks/,manifest.rds}
```

This runs the fixed setup once and writes ~80 pre-sliced RCTD objects
(~1.2 GB total). Slicing matters: the fitted setup object carries the full
1.59M-cell query and is ~19 GB resident. Shipping *that* to 80 tasks would need
~1.5 TB of RAM; a slice needs ~1–2 GB.

### 2. Transfer

```bash
rsync -avP ~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd/ \
    user@cluster:/scratch/$USER/hpc_rctd/
rsync -avP xenium/label_transfer/hpc/ user@cluster:/scratch/$USER/hpc_rctd/scripts/
```

### 3. On the cluster — environment (once)

```bash
micromamba create -f scripts/environment.yml
micromamba activate xenium-rctd
Rscript -e 'remotes::install_github("dmcable/spacexr", upgrade = "never",
                                    lib = "PATH_TO_ENV/lib/R/library")'
```

**The `lib =` argument above is not optional — confirmed as a real failure
mode, not a hypothetical one.** On any account where `R_LIBS_USER` is set
(commonly via `~/.Renviron`), `remotes::install_github()` installs into that
user library instead of the active conda/micromamba env, silently — the
install appears to succeed, but `library(spacexr)` inside the env then fails
with "there is no package called 'spacexr'" the first time an array task
actually needs it. `~/.Renviron` re-asserts `R_LIBS_USER` even when the
environment variable is overridden elsewhere in the shell, so activating the
env is not sufficient protection on its own. Get `PATH_TO_ENV` from
`micromamba env list` (or `echo $CONDA_PREFIX` right after activating) — e.g.
`~/micromamba/envs/xenium-rctd/lib/R/library`. Confirm with
`Rscript -e 'library(spacexr); packageVersion("spacexr")'` after installing,
inside the activated env, before trusting the env is ready — don't assume a
clean `install_github()` exit means the package landed where it needs to.

A related latent hazard, not yet triggered but worth knowing: that user
library has no R-version subdirectory and may be shared across every R
version on the account (e.g. a system R 4.3.2 alongside this env's R 4.4),
which is a cross-version ABI risk if it's ever relied on for anything else.

### 4. Submit

`--partition=medium` and no `--account` are already set (confirmed working on
phoenix/prism via srun + sacct) — edit `rctd_array.sbatch` only if deploying
to a different cluster.

```bash
cd /scratch/$USER/hpc_rctd
mkdir -p logs
sbatch scripts/rctd_array.sbatch /scratch/$USER/hpc_rctd
```

### 5. Gather

```bash
Rscript scripts/03_gather.R /scratch/$USER/hpc_rctd
```

Reports any missing chunks rather than silently returning a short table. To
resubmit only what failed:

```bash
sbatch --array=$(Rscript scripts/04_missing.R /scratch/$USER/hpc_rctd) \
       scripts/rctd_array.sbatch /scratch/$USER/hpc_rctd
```

Then pull `rctd_all.csv.gz` back and feed it to `xenium_label_transfer.qmd`.

## Reference configuration baked into step 1

Both settings were established empirically and matter more than any tuning
parameter here:

- **`ASSAY_REF = "RAW"`**, not the CellBender-corrected assay. CellBender
  sharpens reference profiles by removing ambient RNA, while Xenium cells are
  blurred by segmentation spillover; a sharp reference against a blurred query
  sends everything to the blurriest reference type. Measured on OR52YW26_2_2:
  CellBender used 17 of 49 labels and called **zero** Glut-HVC-1 cells; RAW used
  41 labels and put Glut-HVC-1 at 30× enrichment in the HVC ROI.
- **`EXCLUDE = "Glut-Nido-3"`**, which acts as a sink. Excluding it took HVC
  enrichment from 30× to 90×.

## Known open issue

~14% of cells still land on whichever reference type has the most mixed profile
(`Endo`, or `Glut-Pre-1` once `Endo` is removed). Deleting the sink label only
relocates the sink — verified. These are spillover-blended query cells that
match no pure reference type, which is precisely what RCTD's doublet mode and
`reject` class are designed to handle, and the reason for running RCTD rather
than transfer anchors. Check the `spot_class` distribution in the gathered
output: cells landing in `doublet_certain`/`reject` rather than on a sink label
is the signal that this worked.
