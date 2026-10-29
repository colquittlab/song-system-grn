# GRN material in `main`: inventory and a plan to move development material to `grn-dev`

`main` is to become the publication-associated repository. This note inventories the GRN (SCENIC+) material that is tracked in `main` and proposes how to move
the development part to `grn-dev`. **Nothing here has been applied to `main`.** The classification is provisional (rules in `build_inventory.py`); it is a
starting point for you to edit.

**Settled (2026-10-20):** (1) the paper's GRN uses the hybrid label set, so the older-label subset runs are classified as development (superseded); (2) the astro-oligo analysis is excluded from the revised manuscript; (3) `overexpression/chick_mafb-emx2_run1.qmd` will be replaced by new data/analysis, so it no longer constrains this cleanup. The remaining decisions are listed in section 4.

- Snapshot: `origin/main` at `8a15dc8` (2026-10-20). Per-file table: [`inventory.csv`](inventory.csv) (path, size, category, group, rationale, first added, last modified, last commit subject).
- Rebuild it: `python build_inventory.py [ref]`. Build the removal branch locally: `./make_cleanup_branch.sh` (dry run by default; never pushes).

## 1. At a glance

`grn/` holds **294 tracked files, 39.5 MB**. For scale, all of `main` is about 533 MB: `snrna/` 389 MB (871 files), `multiome/` 93 MB, `grn/` 39.5 MB, `data/` 8 MB, `xenium/` 3 MB.
So the GRN material is a small share of the repository's size; the motive for moving it is clarity, not size.

| category | files | MB | meaning |
|---|---:|---:|---|
| PUBLICATION | 77 | 7.9 | candidate to stay in `main` (config37 pipeline, MAFB/regulon figure code, cisTarget builds) |
| DEVELOPMENT | 162 | 29.0 | to move to `grn-dev`: the four older-label runs (27.6 MB, mostly notebook outputs), the sweep/control configs, diagnostics |
| DECIDE | 55 | 2.6 | needs your call (section 4): hybrid subset runs, cross-config MAFB analyses, configs 34-39 |

### Publication candidates
| group | files | MB |
|---|---:|---:|
| hybrid all-cell: pycisTopic | 4 | 6.2 |
| hybrid all-cell: scenicplus | 3 | 1.1 |
| analysis: MAFB / regulon figure code | 52 | 0.5 |
| hybrid all-cell: config37 | 3 | 0.0 |
| hybrid all-cell: inputs | 2 | 0.0 |
| hybrid all-cell: cisTarget | 8 | 0.0 |
| cistarget build scripts | 2 | 0.0 |

### Development (move to `grn-dev`)
| group | files | MB |
|---|---:|---:|
| older-label subset run: ra-arco-hvc-nc_glut | 9 | 8.4 |
| older-label subset run: ra-arco-hvc-nc_gaba | 9 | 8.0 |
| older-label subset run: ra-arco-hvc-nc_astro-oligo | 9 | 6.5 |
| older-label subset run: glut-ra_glut-arco-1 | 9 | 4.7 |
| hybrid all-cell: other configs | 99 | 1.1 |
| analysis: older-label report wrappers | 4 | 0.1 |
| analysis: sweep / control diagnostics | 10 | 0.1 |
| hybrid all-cell: pycisTopic | 4 | 0.0 |
| hybrid all-cell: scenicplus | 1 | 0.0 |
| hybrid all-cell: inputs | 2 | 0.0 |

### To decide
| group | files | MB |
|---|---:|---:|
| glut hybrid reference run | 28 | 1.2 |
| analysis: coherence tiers (DEVELOPMENT) | 3 | 1.2 |
| hybrid all-cell: other configs | 15 | 0.2 |
| analysis: single-config report (PUBLICATION: template, `assemble_slim_*`) | 3 | 0.05 |
| analysis: config-sweep reports (DEVELOPMENT) | 3 | 0.05 |
| analysis: cross-config MAFB analysis | 5 | 0.1 |
| hybrid all-cell: scenicplus | 5 | 0.0 |
| analysis: local helpers | 2 | 0.0 |
| hybrid all-cell: pycisTopic | 2 | 0.0 |

## 2. What the evidence says (why the categories look like this)

- **`config37` is the main run.** Commit `7a134c2` "Make config37 the main comparison config". It uses the hybrid cisTarget database and the ArchR song-pair DARs, and all the MAFB
  shared-program / regulon-quadrant / region-heatmap scripts are built around it (their outputs are named `*_config37*`). The other 38 configs are the sweep that led there:
  1-11 one-at-a-time thresholds, 12-13 normalization control (the old double-normalised `.raw`), 14-15 topic number, 16-28 a Latin-hypercube threshold sweep,
  29-33 GABA/LGE region-set removal test, 34-35 the 40-topic object, 36-39 the hybrid-database / ArchR-pair comparisons (from the comments in `ra-arco-hvc-nc/scenicplus/make_configs.py`).
- **The MAFB analyses are manuscript figure code, not scratch.** The unmerged branch `worktree-grn-k40-nosongpairs` is 21 commits ahead of `main` and its subjects read "Panel A expansion", "Fig. 4D track rows",
  "KCNC1 expression bars", UR strip plots. Those files live in `grn/analysis/` and `scripts/`, so they belong in `main` (see Order of work).
- **The README must be rewritten around the hybrid set.** It still lists `ra-arco-hvc-nc` as "config3", `ra-arco-hvc-nc_glut` "config15", `_gaba` "config13", `glut-ra_glut-arco-1` "config18"
  (the older labels), while `config/paths.R` names config14 / config13 / config11 for the same runs, and the main result is config37 in `ra-arco-hvc-nc` (hybrid).
- **The older-label runs are superseded and nothing blocks removing them.** `overexpression/chick_mafb-emx2_run1.qmd` reads their SCENIC+ *outputs* (hdd run directories named in the `SCENICPLUS_*_CONFIG_DIR`
  constants in `config/paths.R`: `glut-ra_glut-arco-1`, `ra-arco_glut` and an older all-cell `ra-arco-hvc-nc_seurat-...` run), not the code in these directories, so removing the code cannot break it; and that analysis is being replaced.
  The astro-oligo run and its report are out of the revised manuscript. Two things to tidy when the new analysis lands: the `SCENICPLUS_*_CONFIG_DIR` constants and their config14/13/11 comments in `config/paths.R`.
  (`ra-arco_glut` never had a directory under `grn/`, and `grn/ra-arco-hvc-nc` in `main` is the hybrid pipeline, not the older all-cell run that `SCENICPLUS_EMX2_CONFIG_DIR` points at.)

## 3. Findings that affect the move

- **Notebooks carry their outputs.** Eight notebooks are 0.6-6.2 MB each (29 MB in total, about 74% of `grn/`): the five `pycisTopic.ipynb` (6.2, 6.1, 6.0, 4.8, 3.2 MB) and three `anndata.ipynb`.
  Measured: they are 29.5 MB with outputs and **0.3 MB without** (`nbstripout` saves 29.2 MB if all eight stay; any that move to `grn-dev` need no stripping). Git history keeps the old blobs.
- **Duplicates.** `genome_annotation.tsv` (0.98 MB) exists in 6 directories (4.9 MB redundant), `chromsizes.tsv` in 6, and the Snakefile is copied into 44 config folders (1.1 MB redundant; the configs are *generated* by `make_configs.py`).
  Keeping only config37 and the generator removes most of that.
- **Dependencies of PUBLICATION files on DEVELOPMENT files** (`build_inventory.py` prints them): only 3 of 74 publication files mention a development file or config.
  Two earlier cases were resolved by reclassifying: `make_interneuron_dars_archr.R` is called by three analysis scripts, and `scenicplus_mafb_targets.py` plus the unsuffixed `scenicplus_mafb_shared_program*.csv`
  compare MAFB across configs 1-35 (so they are DECIDE, not publication). The three remaining mentions are benign: a comment in `export_h5ad_hybrid.R` naming `make_adata_doublenorm.py`; `cisTarget/prism/README.md` naming `make_configs.py` (recommendation 7 keeps a trimmed version); and config37's YAML naming `chromsizes.tsv` / `genome_annotation.tsv` (a false match: it points at the `ra-arco-hvc-nc` copies, which are publication files; the development copies are the duplicates under `ra-arco-hvc-nc_glut_hybrid`). The string `config35` also appears as a comparison row inside the `*_config37*` CSVs, which is why config35 is DECIDE.
- **Reproducing config37 needs more than config37.** Its YAML points at `shared_k40/`, which `config34` builds (accessibility+expression object and region-to-gene fit); `prepare_shared_upstream.sh` makes it. Config34 and that script are therefore marked DECIDE/keep-candidates.
- **GRN-adjacent files outside `grn/`** (not in the inventory table; check paths when moving): README rows for `grn/`; `config/paths.R` and `paths.sh` (the `SCENICPLUS_*` constants, machine-specific `/hdd` paths); `envs/scenicplus.yaml`;
  `scripts/scenicplus.R` (sourced by the report templates and `overexpression/`); `data/lambert_tfs/tfs.csv`; `overexpression/chick_mafb-emx2_run1.qmd` (being replaced); `multiome/archr/` (the ArchR DARs that config37 uses).
- **Branches.** Unmerged GRN-related work: `worktree-grn-k40-nosongpairs` (21 commits, figure code), `worktree-grn-glut-hybrid` (PR #11 -> `grn-dev`, 15 commits, development material). `worktree-hybrid-glut-labels` and `worktree-phylop-perbase` are not GRN.
  Merged `worktree-*` branches (archr-hybrid, grn-hybrid-rerun, hybrid-magic-sct, rename-hvcra-int, add-adult-finch-integration) can be deleted on the remote once you no longer need them.

## 4. Decisions needed (with a recommendation)

Resolved: **label set = hybrid** (2026-10-20).

| # | decision | recommendation |
|---|---|---|
| 1 | The older-label directories (36 files, 27.6 MB) and their report wrappers | Move to `grn-dev` now (nothing depends on their code; astro-oligo is out of the paper, the chick analysis is being replaced). Tidy the `SCENICPLUS_*_CONFIG_DIR` constants in `config/paths.R` when the new analysis lands. |
| 2 | Which subsets does the paper need in hybrid form? The README listed glut, gaba, astro-oligo and RA vs Arco-1; astro-oligo is now out. Hybrid versions: glut (in `main`); gaba, RA vs C1H-1 (with and without run3) and all-CACNA1H (on `grn-dev` via PR #11) | Promote only what a published figure or table uses: `git checkout grn-dev -- <dir>` into a normal PR against `main`. |
| 3 | The hybrid glut run (`ra-arco-hvc-nc_glut_hybrid`, 28 files in `main`) | Keep if the paper uses a glut-only SCENIC+ run, otherwise `grn-dev`. |
| 4 | Hybrid report templates | **Decided 2026-10-20:** the sweep reports (`all.qmd`, `glut.qmd`, `scenicplus_hybrid_glut_template.qmd`) go to `grn-dev`. `main` keeps the single-config report: `scenicplus_hybrid_all_template.qmd` (it takes `params$configs`) with the new wrapper `config37.qmd` (`configs: ["config37"]`, not in the inventory because it did not exist at `8a15dc8`). |
| 5 | Local helpers (`assemble_slim_*`) | **Kept** (the single-config report reads the slim h5mu they build). |
| 6 | Coherence tiers (3 files, 1.2 MB) | **Decided 2026-10-20:** development, move to `grn-dev`. |
| 7 | Cross-config MAFB analyses (5 files) | Keep only if the config-robustness analysis is in the supplement; they need `config_parameters.tsv` and the sweep results, which would otherwise be gone. |
| 8 | Configs 34, 35, 36, 38, 39 | Keep config34 (needed to rebuild config37's upstream) or describe it in the README; 35, 36, 38, 39 to `grn-dev` unless cited. |
| 9 | `make_configs.py` | Keep a trimmed version that defines config37 only; the full 39-config generator goes to `grn-dev`. |
| 10 | Cluster helpers (`submit_all.sh`, `link_results_to_hdd.sh`, `prepare_shared_upstream.sh`) | Keep `prepare_shared_upstream.sh` (config37 needs it); the other two to `grn-dev`. |

Resolve a decision by editing the row's `category` in `inventory.csv` to `PUBLICATION` or `DEVELOPMENT` (edit the rules in `build_inventory.py` if you prefer to regenerate).

## 5. Strategy

Principle: **`grn-dev` keeps everything; `main` loses the development part through one reviewable PR; nothing is rewritten or force-pushed.**

Order of work:
1. **Land PR #11 on `grn-dev`** (it already targets it). `grn-dev` then holds main plus all development work.
2. **Settle the figure branch first.** `worktree-grn-k40-nosongpairs` (figure code, touches `grn/analysis/`) should go to `main` *before* the cleanup, so the cleanup PR does not conflict with it.
3. **Tag the pre-cleanup state:** `git tag grn-pre-cleanup origin/main && git push origin grn-pre-cleanup` (a permanent pointer to everything).
4. **Resolve the decisions** in section 4 by editing `inventory.csv`.
5. **Build the cleanup branch locally**, review, push, open a PR against `main`:
   ```bash
   cd grn/cleanup
   ./make_cleanup_branch.sh                 # dry run: counts and directories
   APPLY=1 ./make_cleanup_branch.sh         # local branch main-cleanup from origin/main, git rm, one commit
   git push -u origin main-cleanup          # then open the PR with base main
   ```
6. **In the same or a follow-up PR on `main`:** rewrite the `grn/` rows of `README.md` around the published set; strip notebook outputs (`nbstripout <notebooks that stay>`); drop duplicate `genome_annotation.tsv` / `chromsizes.tsv` copies (keep one, referenced by path); reconcile the `SCENICPLUS_*_CONFIG_DIR` comments in `config/paths.R` (config14/13/11 vs the README).
7. **Keep `grn-dev` intact after the PR merges (important, see below).**
8. **At publication:** tag `main` (e.g. `v1.0-manuscript`) and archive that tag (Zenodo gives a DOI).

### The merge trap, and the recipe
After the cleanup PR merges, `main` has *deleted* the files that `grn-dev` still has. A plain `git merge main` into `grn-dev` would apply those deletions and remove the development material from `grn-dev`.
Use the revert recipe once, right after the PR merges:
```bash
git switch grn-dev
git merge main                 # brings in the cleanup commit C (and anything else new in main)
git revert --no-edit <C>       # puts the files back on grn-dev; later merges of main will not delete them again
```
Fixes made in `main` still flow into `grn-dev` this way. Never merge `grn-dev` into `main` wholesale after the cleanup; take individual publication-worthy files with `git checkout grn-dev -- <path>` in a normal PR.

### What this does not do
- It does not shrink or clean the **history**: removed files stay in `main`'s past (the whole pack is only 64 MiB; the largest blobs in history are not GRN: `snrna/.../chicken_obs.csv` 14 MB, two ~11 MB peak BEDs under `multiome/archr`).
  If you ever want a history without development material or large blobs, that is a separate and destructive step (`git filter-repo` plus a force-push, or publishing a fresh snapshot repo). It is not needed for a normal cleanup and I would not do it without a backup and your explicit go-ahead.
- It does not touch the non-GRN parts of the repository (`snrna/` is 73% of its size).

## 6. Files here

| file | what |
|---|---|
| `inventory.csv` | one row per tracked file under `grn/` in `origin/main`, with category, group and rationale |
| `build_inventory.py` | regenerates the inventory and runs the dependency scan; the classification rules live in `classify()` |
| `make_cleanup_branch.sh` | dry-run-by-default builder of the local `main-cleanup` branch (`CATS`, `BASE`, `BRANCH`, `APPLY` env vars); tested in a throwaway clone (144 files removed, `grn/` 296 -> 152 files) |
