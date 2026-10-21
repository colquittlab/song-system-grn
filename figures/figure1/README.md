# Figure 1 — standalone panel scripts

One script per data panel of Figure 1 (`fig1_snrna_v10.svg`). Each regenerates its panel from the outputs of the upstream pipelines, without running
any other panel or notebook. Shared input paths and helpers are in `common.R`, which every panel script sources; it in turn sources the repository's
`config/paths.R` and `config/figure_theme.R`. The two helpers the original scripts took from the lab's shared utilities (`/opt/colquittlab/utils`) are
copied into `common.R`, so that directory is not needed. Panels are written to `figures/figure1/output/`, which is not tracked; set `FIG1_OUT` to write elsewhere.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| A | Sensorimotor and song circuits (photographs, schematic) | — (drawing) | — |
| B | Song system (schematic) | — (drawing) | — |
| C (top) | Microdissection, pooling, snRNA-seq and genotype demultiplexing (schematic) | — (drawing) | — |
| C (bottom right) | Xenium transcript density of one section (10 µm bins, 20 µm smoothing) | `panel_C_xenium_transcript_density.R` (+ `panel_C_bin_transcripts.py`) | `panel_C_transcript_density_OR52YW26_1_4_10um_smooth20.{pdf,png}` |
| D | snRNA-seq UMAP coloured by cell type | `panel_DE_snrna_umap.R` | `panel_D_umap_celltype.{pdf,png}` |
| E | snRNA-seq UMAP coloured by dissected region | `panel_DE_snrna_umap.R` | `panel_E_umap_position.{pdf,png}`, `panel_DE_umap_embedding.csv` |
| F | Xenium positions of GABA-MGE-PVALB-1/2/3 and Astro-1/2/3 cells (label transfer) | `panel_F_xenium_spatial_maps.R` | `panel_F_spatial_<type>_gte400.png` |
| G | Marker dot plot; relative fraction of each type per region; snRNA-seq and spatial regional specificity | `panel_G_markers_regional_specificity.R` | `panel_G_dotplot_manual_markers.pdf`, `panel_G_region_fraction_snrna_specificity.pdf`, `panel_G_spatial_specificity_aligned_bars_allcells.{pdf,png}`, `.csv` tables |

## Running

```bash
Rscript figures/figure1/panel_DE_snrna_umap.R       # run from anywhere inside the repository (here::here)
```

R ≥ 4.3 with the packages in `renv.lock` (Seurat 5, scCustomize, qs2, ComplexHeatmap, circlize, egg, cowplot, nanoparquet, Matrix, viridisLite).
Panel C needs Python 3 with numpy and pyarrow only if the binned transcript table is missing (`FIG1_PYTHON` names the interpreter).
Run times on the lab server: C 20 s; D/E 70 s; F 90 s; G 40 s. D/E and G read the 3 GB snRNA-seq object and need ~25 GB of memory.

## Inputs (upstream pipelines)

All paths are set in `common.R`; each can be overridden with the environment variable named there.

| Input | Made by | Used by |
|---|---|---|
| snRNA-seq Seurat object with the hybrid division-based cell-type labels (`snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2`; `FIG1_SNRNA_OBJ`) | `snrna/naming/hybrid_division_naming.qmd` | D, E, G |
| Proseg segmentation of each Xenium section (`<XENIUM_PROC_DIR>/proseg/<section>/`: `cell-metadata.parquet`, `expected-counts.csv.gz`, `transcript-metadata.parquet`; `FIG1_PROSEG_DIR`) | `xenium/preprocessing/hpc_proseg/` | C |
| Transcripts binned at 10 µm (`xenium/transcript_density/tx_density_OR52YW26_1_4_10um.csv.gz`; `FIG1_TX_DENSITY_DIR`) | `xenium/bin_transcripts.py` (copied here as `panel_C_bin_transcripts.py`, run by the panel C script if the table is missing) | C |
| RCTD label transfer onto the proseg cells (`<XENIUM_PROC_DIR>/hpc_rctd_proseg_hybrid/`: `rctd_all.csv.gz`, `proseg_cell_metadata.csv.gz`, `proseg_ncount.csv.gz`; `FIG1_RCTD_DIR`) | `xenium/label_transfer/hpc_rctd_proseg_hybrid/` | F, G |

## Notes on individual panels

- **A, B, top of C, left and middle of the bottom of C, schematics in F.** Drawings and photographs; no code.
- **C.** Same as `xenium/transcript_density_smoothed.R OR52YW26_1_4`, which also draws other bin sizes and a per-cell map. Transcripts with qv ≥ 20
  (proseg's observed, not repositioned, coordinates) are counted in 10 µm bins. Smoothing is a Gaussian with sd 20 µm, normalized by the
  smoothed occupancy so that tissue edges are not darkened. Bins with smoothed occupancy < 0.15 are dropped. The colour scale is clipped at the 2nd
  and 98th percentiles of in-tissue bins.
- **D, E.** Same calls as `snrna/reduction_viz/combined_all_umap_hybrid.R`, which fits a grid of UMAPs; the figure uses `dims40nn30mindist0.3`
  (PCA on the SCT assay; UMAP on PCs 1–40, 30 neighbours, min.dist 0.3). The script refits that embedding from the labelled object. The refit equals
  the one saved by the original script exactly (PCA and UMAP coordinates, maximum difference 0). Cell-type and region names were added by hand.
- **F.** Same as the "every cell call in the taxonomy (gte400)" section of `xenium/label_transfer/proseg_spatial_maps_hybrid.qmd`. Confident RCTD
  calls (singlet or confident doublet) with ≥ 400 transcripts are shown in black over all cells in grey. Each section is rotated/flipped to
  posterior-left, dorsal-up. The script draws all ten sections. The figure crops OR52YW26_1_4 (1.2 mm from the midline) for the PVALB types and
  OR52YW26_1_7 (2.2 mm) for the astrocyte types.
- **G.**
  - Dot plot: same as `snrna/deg/dotplot_manual_markers_hybrid.qmd`.
  - Region fraction and snRNA-seq specificity: same as `snrna/trees/celltypes_hclust_all_hybrid.qmd` (`cluster_proportions_vert_with-bar_ordered.pdf`).
    The region × cell-type count table is normalized to region totals, then to cell-type totals. Specificity is the KL divergence of each type's
    normalized profile from uniform.
  - Spatial specificity ("spatial trans."): same as `xenium/label_transfer/spatial_variation_by_type.R 0` plus `spatial_variation_aligned.R 0`
    (`spatial_variation_aligned_bars_allcells.pdf`). All confident RCTD calls are used, with no transcript floor. Cells are binned on a 400 µm grid.
    For each type and section, the score is the Jensen–Shannon divergence of 100 subsampled cells from all cells, minus the same for 100 random
    cells (30 draws each, seed 1). It is averaged over the sections with ≥ 100 cells of the type. Glut-DACH2-HVCra-Int reaches 100 cells in one
    section only, so its score rests on that section.
  - The figure takes the heatmap and snRNA-seq bars from the second file and the spatial bars from the third. The legends and the "rel. fraction in
    region" label were restyled by hand.

## Transcript floor of the Xenium analyses (panels F and G)

The two Xenium analyses use different cell sets, as in the figure:
- the panel F maps keep confident RCTD calls with ≥ 400 transcripts (392,675 cells);
- the panel G spatial specificity keeps every confident call (921,376 cells; none lacks a transcript count).

`FIG1_MIN_TRANSCRIPTS` overrides both. Outputs are tagged `_gte<N>` or `_allcells`. An earlier version of the figure used ≥ 400 for panel G too:

```bash
FIG1_MIN_TRANSCRIPTS=400 Rscript figures/figure1/panel_G_markers_regional_specificity.R   # earlier panel G (Glut-DACH2-HVCra-Int N.D.)
FIG1_MIN_TRANSCRIPTS=0 Rscript figures/figure1/panel_F_xenium_spatial_maps.R               # maps of every confident call
```

How the all-cells version compares with the ≥ 400 version:

- **F.** The all-cells maps are pixel-identical to the no-floor maps of `proseg_spatial_maps_hybrid.qmd` (`spatial_all_<type>.png`).
- **G, scores.** The spatial-specificity scores of the 46 types scored in both versions agree closely (Spearman ρ = 0.98, Pearson r = 0.98).
  Without the floor, scores are about 10 % lower (median ratio 0.90).
- **G, correlation with snRNA-seq.** It is unchanged: ρ = 0.78, p = 1.5 × 10⁻¹⁰ (n = 47), against ρ = 0.77, p = 2.8 × 10⁻¹⁰ (n = 46) with the floor.
- **G, by class.** Glutamatergic types still score highest (mean 0.28) and GABAergic lowest (0.06); glutamatergic vs GABAergic Wilcoxon p = 1.7 × 10⁻⁷.
  Astrocyte types drop the most (mean 0.17 → 0.11, mainly Astro-1 and Astro-2), and Endo drops from 0.075 to 0.021.
- **G, Glut-DACH2-HVCra-Int.** It is scored in this version (0.15), but from a single section that has exactly the minimum of 100 cells, so its
  score is much less certain than the others'.

## Check against the figure files

Each script was run and its output compared with the file placed in `fig1_snrna_v10.svg`. Raster outputs were compared pixel by pixel, PDFs
rasterized at 100 dpi and their text extracted:

| Panel | Result |
|---|---|
| C | Pixel-identical (PNG and PDF); same text |
| D | Pixel-identical (PNG and PDF) |
| E | Pixel-identical (PNG and PDF) |
| F | All six maps pixel-identical |
| G, dot plot | The file in the figure predates the renaming of Glut-DACH2-HVCra-Pre to Glut-DACH2-HVCra-Int; the figure shows the new name, as does this output. Every dot has the same colour, and the same size to within cairo's 1/64 pt coordinate grid. The plot sits 2.8 pt further right. |
| G, region fraction + snRNA-seq specificity | Same; only ComplexHeatmap's automatic legend name differs (`matrix_1` here, `matrix_2` in the original, which drew another heatmap first) |
| G, spatial specificity | Pixel-identical (PNG and PDF) to `spatial_variation_aligned_bars_allcells` from `xenium/label_transfer/spatial_variation_aligned.R 0`; per-type scores identical to `spatial_localisation_by_type_allcells.csv` (47 types). Compared by eye with the bars placed in the figure: same bars. The ≥ 400 version (`FIG1_MIN_TRANSCRIPTS=400`) is pixel-identical to the earlier figure's `spatial_variation_aligned_bars.pdf`. |
