# multiome/archr

ArchR analyses of the motor-pathway multiome, run against the **hybrid** cluster
labels.

## Where the labels come from

`multiome/seurat/label_transfer_hybrid.qmd` transferred the hybrid snRNA labels
(`snrna/naming/hybrid_division_naming.qmd`, `celltype_hybrid`) onto the multiome
object and curated them to one label per multiome cluster by majority vote. Its
`cluster_hybrid_majority_vote.csv` is the whole of that curation, and it is tracked.

`hybrid_labels.R` applies it here:

```r
source(here::here("config/paths.R"))
source(here::here("multiome/archr/hybrid_labels.R"))

proj = loadArchRProject(ARCHR_PROJ_DIR)
proj = add_cluster_hybrid(proj)     # adds proj$cluster_hybrid, in memory only
```

Reading the per-cluster CSV rather than the 2 GB `obj_clustered_hybrid.qs2` is a
shortcut, and it is a checked one: `verify_cluster_hybrid_map()` confirms the map
reproduces every per-cell label in that object (16,588/16,588 at the time of writing).
`chromvar.qmd` and `pos_regulators.qmd` do read the object itself, for expression.

The helper also carries the things the old label names used to encode implicitly —
`hybrid_ct_order`, `hybrid_pairs_glut`, `hybrid_pairs_gaba`, `hybrid_groups`,
`hybrid_position`. Use those instead of re-deriving groups with `grepl()`: the
pre-hybrid scripts built their families out of the labels themselves (`"Arco|RA"`,
`"HVC|NC|Nido"`, `"Pre"`), and the hybrid names do not support that.

## What the rename did

Mostly one-to-one. The substantive changes:

| pre-hybrid | hybrid | note |
|---|---|---|
| `Glut-RA` | `Glut-CACNA1H-RA` | |
| `Glut-Arco-1` | `Glut-CACNA1H-1` | |
| `Glut-Arco-2`, `-6`, `-7` | `Glut-CACNA1H-2` | **three clusters merged into one** |
| `Glut-HVC-1` / `-2` | `Glut-DACH2-HVCra` / `-HVCx` | |
| `Glut-NC-1`…`-4` | `Glut-DACH2-1`…`-4` | |
| `Glut-Pre-1` / `-2` / `-3` | `Glut-NSC` / `Glut-NB` / `Glut-Im` | |
| `GABA-1-1`, `-1-2` | `GABA-LGE-1`, `-2` | |
| `GABA-2-1`, `-3`, `-4-1` | `GABA-MGE-SST-1`, `-PVALB-1`, `-PVALB-2` | |
| `GABA-5`, `-6`, `-8`, `-Pre` | `GABA-CGE`, `-MGE-LAMP5`, `-MGE-LHX8`, `GABA-Im` | |
| `Astro-1`, `Astro-2` | `Astro` | **two clusters merged into one** |
| `Glut-Nido-3` | *(none)* | artefactual; its 745 cells are dropped |

27 hybrid labels over 16,316 of the project's 17,061 cells. The Astro merge collapsed
`peaks.qmd`'s two astro-vs-oligo contrasts into one.

## Conventions

- **Nothing writes to the ArchR project.** `add_cluster_hybrid()` modifies only the
  in-memory `cellColData`. The peak set is the project's, called once in
  `archr_processing.qmd` grouped by `cluster`, and reused here as a fixed feature
  space — relabeling cells does not change which regions are open. Re-calling peaks
  by `cluster_hybrid` would write a new `PeakCalls/` and `PeakMatrix` into the shared
  84 GB project; if that is ever wanted it belongs in `archr_processing.qmd`, beside
  the call that made the current set. BED exports that used to land in
  `<proj_dir>/beds` now go to each script's own `out_dir`.
- **Outputs go to `<script>_hybrid/`**, so the pre-hybrid results beside them are
  left intact for comparison.
- `archr_processing.qmd` builds the project, so its import, pseudobulk-coverage and
  peak-calling steps stay on `cluster` — that column is the provenance the hybrid
  labels were transferred onto. Only its analysis half (from "Differential peaks")
  regroups.

## Run them one at a time

These scripts share the ArchR project's HDF5 Arrow files, and `chromvar.qmd`'s
`addDeviationsMatrix` *writes* MotifHSMatrix into them. Rendering two at once fails:
`differential_accessibility.qmd` died 51 minutes in with `H5Fopen(): Unable to open
file` while chromvar was writing, and `pos_regulators.qmd` failed reading the matrix
chromvar was mid-write on. Concurrent renders also collide in `multiome/archr/.quarto`,
which fails the HTML write *after* the R code has already finished. Render serially.

`differential_accessibility.qmd` takes ~90 minutes uncontended (351 pairwise
contrasts) and its `peaks_glut.rds` is **26 GB** -- check free space first.
`peaks_full.qmd` is byte-identical to it apart from `script_name` and one output
filename, so running both writes that 26 GB twice for the same result.

## Note on the notebooks

These were written as interactive scratchpads and several were not runnable top to
bottom: `peaks_glut.qmd` and `peaks_gaba.qmd` used `mat_avg` and `fdr_thresh` ~120
lines above where they were defined and saved a data frame before building it, and
`peaks2geneslinks.qmd` set `names(proj_regions)` from an unnamed vector so its
per-region BED loop wrote nothing. Those were reordered/fixed so the files run as
written; no computation changed. `differential_accessibility.qmd` and
`peaks_full.qmd` are still the same script under two names, differing only in
`script_name` and one output filename.
