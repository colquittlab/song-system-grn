# Figure 3 — standalone panel scripts

Scripts for the code-made parts of Figure 3 (`fig3_multiome_overview_v6.svg`). Shared input paths and helpers are in `common.R`, which every R panel
script sources; it in turn sources the repository's `config/paths.R` and `config/figure_theme.R`. Panels are written to `figures/figure3/output/`, which
is not tracked; set `FIG3_OUT` to write elsewhere.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| C | Chromatin accessibility at the SIX2 locus (genome-browser tracks) and SIX2 expression bars | `panel_C_six2_expression_bars.R` (bars only) | `panel_C_six2_expression_bars_no_points.svg` (+ `panel_C_six2_expression_bars.svg` with per-bird dots, `.csv`) |

## Running

```bash
Rscript figures/figure3/panel_C_six2_expression_bars.R    # ~2 min (reads the Seurat object)
```

## Inputs

| Input | Made by | Used by |
|---|---|---|
| Hybrid-labeled multiome Seurat object (`multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2`; override with `FIG3_SEURAT_OBJ`) | `multiome/seurat/reduction_viz/combined_all_umap_hybrid.R` | C |

## Notes on individual panels

- **C.** The accessibility tracks (pseudobulk ATAC coverage of RA, C1H-1, HVCra, D2-1, HVCx, D2-4 at chr3:16.505-16.532 Mb, with zoom insets of three
  regions) are a genome-browser export (`hgt_genome_*_six2*.pdf` beside the figure file), not scripted. In v6 the panel's layer was slimmed with
  `scripts/browser_pdf_to_svg.py --within layer2` (as for the Figure 4E tracks): the clipped copies of whole tracks behind the zoom insets were dropped and
  the bars of each track merged into one path, 45,691 -> 27 paths (file 17.8 -> 7.5 MB); renderings before and after differ only in anti-aliasing at
  bar edges. The expression bars are SCT log-normalized SIX2 expression (natural log), mean of per-bird means (bird = souporcell assignment within a
  library, >= 10 cells), one row per track in the track order and colors, each row labeled; same construction as the KCNC1 bars of Figure 4E. They sit
  under the gene model, level with the zoom-inset rows (translate(154.5, 40.5) mm in layer2). SIX2 is expressed in the three song neurons only:
  RA 0.14 (19 % of cells), HVCra 0.07 (7 %), HVCx 0.23 (29 %); C1H-1, D2-1, D2-4 <= 0.001 (<= 0.15 %).
