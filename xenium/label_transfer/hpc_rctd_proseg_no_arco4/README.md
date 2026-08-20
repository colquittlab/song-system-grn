# RCTD over proseg data — Glut-Arco-4 excluded from the reference

Identical to `xenium/label_transfer/hpc_rctd_proseg/` (proseg-resegmented
Xenium data, all ten sections, RAW-assay reference, Glut-Nido-3 already
excluded there) with one change: **`Glut-Arco-4` is also excluded** from the
reference. See that directory's README for the full rationale on segmentation,
sizing, and the two hard-won reference settings (`RAW` assay, Glut-Nido-3
exclusion) this run inherits unchanged.

Everything else — chunk size (20k), `UMI_min`/`counts_MIN` (20), max
1000 reference cells/label, section list, gene re-indexing procedure, sbatch
sizing, self-contained/rclone workflow — is unchanged. `gene_lists/` is a
plain copy (proseg wasn't re-run, so the per-section gene orders are
identical).

## Steps

Same as `hpc_rctd_proseg/`, pointed at this directory instead:

```bash
cd <repo root>
Rscript xenium/label_transfer/hpc_rctd_proseg_no_arco4/01_prepare_chunks.R
# stages $XENIUM_PROC_DIR/hpc_rctd_proseg_no_arco4/

D=~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_no_arco4
rclone sync "$D" remote:hpc_rctd_proseg_no_arco4
# on the cluster:
rclone copy remote:hpc_rctd_proseg_no_arco4 /scratch/$USER/hpc_rctd_proseg_no_arco4
cd /scratch/$USER/hpc_rctd_proseg_no_arco4
micromamba activate xenium-rctd   # reuse the existing env, nothing new to install
mkdir -p logs
sbatch scripts/rctd_array.sbatch .
```

The `--array` count no longer needs manual reconciliation: `01_prepare_chunks.R`
patches the staged `scripts/rctd_array.sbatch`'s upper bound to the actual
chunk count automatically after slicing (this run: 79, matching the proseg
sibling — one fewer reference label changes the fit-cell count only
trivially). `02_fit_chunk.R` also exits cleanly rather than erroring if the
array range is ever wider than the data, as a second line of defense.

Gather and pull back exactly as before:

```bash
Rscript scripts/03_gather.R .
rclone sync /scratch/$USER/hpc_rctd_proseg_no_arco4 remote:hpc_rctd_proseg_no_arco4
# locally:
rclone copy remote:hpc_rctd_proseg_no_arco4/rctd_all.csv.gz ~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_no_arco4/
```

## Comparing against the original run

Both `rctd_all.csv.gz` files share the same cell-naming convention
(`<section_id>_proseg_<cell>`), so composition, ROI enrichment, and the
spatial-map scripts can be pointed at either output directory with no other
change. The interesting comparison is what happens to the cells that were
previously called `Glut-Arco-4` — same check as the Glut-Nido-3 removal test:
do they land on a real neighboring type, or reveal a new sink.
