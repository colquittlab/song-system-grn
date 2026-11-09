# Supplementary tables

`export_supplementary_tables.py` writes one `.xlsx` workbook per supplementary table listed in
the Google Sheet "adult-multiome_science_tables". Each workbook opens with a README sheet that
records its source files, method, thresholds and the commit that produced it.

```bash
conda run -n cpdb python tables/export_supplementary_tables.py            # all nine
conda run -n cpdb python tables/export_supplementary_tables.py --only 07  # one table
```

Output goes to `<SSG_DATA_ROOT>/tables/output/` (gitignored); override with `--out` or `SSG_TABLE_OUT`.
Inputs are read from `SSG_DATA_ROOT` (default: the main checkout, because most inputs are
gitignored or untracked analysis outputs that exist only there) and `SSG_SCENICPLUS_ROOT`.

| # | Table | Figure | Source |
|---|---|---|---|
| 01 | Finch-mouse integration labels | fig1-s3 | `snrna/integration/composite_scoring/results/yao_adult_hybrid/` |
| 02 | Finch-chicken integration labels | fig1-s2 | `snrna/integration/composite_scoring/results/gg_adult_hybrid/` |
| 03 | Glut differentially expressed genes | fig2-s1 | `snrna/deg/song-surround_deg_glut_hybrid/song_surround_paired_deg.csv` |
| 04 | Song system core genes | fig2-s1 | `snrna/deg/song-surround_deg_glut_hybrid/core_{no_lmansh_,}gene_sets.csv` |
| 05 | Glut finch-chicken label transfers | fig2-s3 | `snrna/integration/gg_glutonly_hybrid/composite/results/composite_plusSATB2_noPre/` |
| 06 | Song/non-song DARs | fig3-s1 | `multiome/archr/peaks_glut_hybrid/beds/` |
| 07 | eRegulons, full dataset | fig4-s1 | SCENIC+ `ra-arco-hvc-nc_hybrid` config37 |
| 08 | eRegulons, RA / C1H-1 | fig6 | SCENIC+ `glut-ra_glut-c1h-1_allruns_hybrid` config2 |
| 09 | eRegulons, GABA | fig6 | SCENIC+ `ra-arco-hvc-nc_gaba_hybrid` config2 |

Caveats on inputs:
- 01/02: the source CSVs keep raw GABA cluster IDs; the script renames them to the hybrid labels
  with `snrna/naming/hybrid_division_naming/hybrid_label_lookup.csv`, as the figures do.
- 05: the composite matching the figure (SCT SAMap, no HVCra-Int) is an uncommitted output in the
  main working tree; the copy tracked in git is an older version.
- 03/04: the producing notebook `snrna/deg/song-surround_deg_glut_hybrid.qmd` is on branch
  `worktree-hybrid-glut-labels`, not yet on main.
- 07-09: only +/+ and -/- eRegulons are exported; direct and extended are pooled (`is_extended`).
