# Figure 2 — standalone panel scripts

One script per data panel of Figure 2 (`fig2_snrna_glut_v10.svg`). Each regenerates its panel from the outputs of the upstream pipelines, without
running any other panel or notebook. Shared input paths and helpers are in `common.R` (R panels A–C) and `common.py` (Python panels E–G); each panel
script uses only its own common file from this folder. The R side sources the repository's `config/paths.R` and `config/figure_theme.R`.

Code from outside this repository is copied into the common files, so neither source is needed:
- the dissection-position colours from the lab's shared utilities (`/opt/colquittlab/utils`);
- the reciprocal transfer score from `finch-integration-toolkit` (a private lab repository).

Panels are written to `figures/figure2/output/`, which is not tracked; set `FIG2_OUT` to write elsewhere.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| A | Trees of the nidopallial (DACH2) and arcopallial (CACNA1H) glutamatergic types after removing song-type markers, with bootstrap support | `panel_A_glut_trees.R` | `panel_A_tree_nido.pdf`, `panel_A_tree_arco.pdf`, `panel_A_tree_support.csv`, `panel_A_features.csv` |
| B (UMAPs) | Song and non-song glutamatergic types embedded without and with CCA integration | `panel_B_song_nonsong_integration.R` | `panel_B_umap_by_type_bare.{pdf,png,svg}` (+ `_songlabels` with titles), `panel_B_umap_coordinates.csv.gz` |
| B (heatmap) | "Integration partner": reciprocal transfer score between each song type and each non-song type | `panel_B_song_nonsong_integration.R` | `panel_B_reciprocal_score_bare.{pdf,png,svg}` (+ version with colour key), `panel_B_reciprocal_scores.csv` |
| C | Xenium positions of the five song types and their non-song partners | `panel_C_xenium_spatial_maps.R` | `panel_C_spatial_<type>_gte400.png` |
| D | Song-system circuit with the partner types (schematic) | — (drawing) | — |
| E | Best reciprocal transfer score of each finch type to any chicken type, across integration strengths (k.anchor) | `panel_E_rpca_sweep_heatmap.py` | `panel_E_rpca_sweep_heatmap.{pdf,png,csv}` |
| F | Finch × chicken reciprocal transfer scores at k.anchor 20 | `panel_F_transfer_heatmap.py` | `panel_F_transfer_heatmap_ka20.{pdf,png,csv}` |
| G | Joint finch–chicken RPCA embedding at k.anchor 20, by finch region and chicken class | `panel_G_rpca_umap.py` | `panel_G_rpca_umap_ka20.{pdf,png,csv}` |
| H | Cell-type evolution model (schematic) | — (drawing) | — |

## Running

```bash
Rscript figures/figure2/panel_A_glut_trees.R      # R panels: run from anywhere inside the repository (here::here)
python figures/figure2/panel_E_rpca_sweep_heatmap.py
```

R needs ≥ 4.3 with the packages in `renv.lock`: Seurat 5 (with DESeq2, used by `FindAllMarkers`), qs2, pvclust, ape, ggtree, ggrepel, ggrastr,
egg, cowplot, ragg, svglite, RColorBrewer.

Python needs Python 3 with matplotlib, seaborn 0.13, pandas, numpy and the Arial font. Text measurement depends on the FreeType library matplotlib
is built with, which moves legend rows by up to 2 pt:
- E and G were drawn with a pip-wheel matplotlib (bundled FreeType 2.6.1);
- F was drawn with a conda-forge matplotlib (system FreeType 2.14).

Run times on the lab server: A 5.5 min (~18 GB); B 10 min (~19 GB); C 2 min; E–G a few seconds each.

## Inputs (upstream pipelines)

All paths are set in `common.R` / `common.py` or at the top of the panel script; each can be overridden with the environment variable named there.

| Input | Made by | Used by |
|---|---|---|
| snRNA-seq object with the hybrid cell-type labels (`snrna/naming/hybrid_division_naming/obj_hybrid_labels.qs2`; `FIG2_SNRNA_OBJ`) | `snrna/naming/hybrid_division_naming.qmd` | A |
| Glutamatergic subset of that object (`snrna/trees/celltypes_hclust_glut_robustness_hybrid/obj_clustered_glutamatergic.qs2`; `FIG2_GLUT_OBJ`) | `snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd` | B |
| Optional: the saved pvclust run of the figure (`snrna/trees/celltypes_hclust_glut_robustness_hybrid/correlation_ward.D/pvclust_sweep.qs2`; `FIG2_PVCLUST_SWEEP`) | `snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd` | A (see notes) |
| RCTD label transfer onto the proseg cells (`<XENIUM_PROC_DIR>/hpc_rctd_proseg_hybrid/`; `FIG2_RCTD_DIR`) | `xenium/label_transfer/hpc_rctd_proseg_hybrid/` | C |
| Finch × chicken label-transfer matrices (`*_matrix.csv`) and joint embedding (`*_umap.csv`) of the RPCA runs, tracked in `snrna/integration/gg_glutonly_hybrid/rpca/results/cca/` (`FIG2_RPCA_RESULTS`) | `R/cca_transfer.R` and `R/cca_transfer_umap.R` of `finch-integration-toolkit` (private; run on the cluster with its `run_cca.sbatch` / `run_cca_umap.sbatch`), on the subsets made by `snrna/integration/gg_glutonly_hybrid/prepare_satb2_test_subsets.py` | E, F, G |

## Notes on individual panels

- **A.** Same computation and drawing as `snrna/trees/celltypes_hclust_glut_robustness_hybrid.qmd` at its manuscript threshold (page 4 of
  `correlation_ward.D/trees_{nido,arco}_thresholds.pdf`).
  - Within each division (DACH2: nidopallial; CACNA1H: arcopallial), markers are called on pseudobulk samples (cell type × bird; DESeq2).
  - Song-type markers with adjusted p < 0.1 are removed. The top 100 up- and down-regulated markers per type (adjusted p < 0.1) are clustered: gene
    specificity across types, Ward's method (ward.D) on correlation distance.
  - Node labels are pvclust AU support (%, 1,000 replicates). Tips are coloured by the dissected region in the type's original cluster name.
  - The notebook's bootstrap ran in parallel without a seed, so its support values cannot be regenerated. The script seeds the bootstrap
    (`iseed = 1`); on identical trees, its values differ from the figure's by up to 5.6 points. Setting `FIG2_PVCLUST_SWEEP` to the saved run draws
    the figure's values, after checking that the saved tree equals the recomputed one.
- **B.** Song (HVC, LMAN, RA) and non-song glutamatergic neurons are treated as two datasets. Each is normalized with SCTransform, with 3,000
  integration features. This repeats the "ALL" set of `snrna/integration/song_nonsong_integration_hybrid.qmd`.
  - UMAPs (`song_nonsong_umap_hybrid.qmd`): left, the two datasets merged without correction; right, after CCA integration (k.anchor 20, dims 1–30,
    IntegrateData). PCA and UMAP on PCs 1–30 with Seurat defaults. Song types are coloured by nucleus, non-song types grey.
  - Heatmap (`song_nonsong_reciprocal_transfer_hybrid.qmd`): label transfer in both directions (CCA, k.anchor 20, k.weight 50). In each direction,
    the score is the share of one type's cells assigned to the other type × their mean prediction score. The reciprocal score is the square root of
    forward × reverse, so it is high only when each type's cells vote for the other.
  - The script combines the three notebooks. It starts from the glutamatergic object and reproduces the notebooks' saved prep, UMAP coordinates
    and scores exactly.
  - Labels on the UMAPs and the heatmap's colour key were restyled by hand.
- **C.** Same as the "every cell call in the taxonomy (gte400)" section of `xenium/label_transfer/proseg_spatial_maps_hybrid.qmd`. Confident RCTD
  calls (singlet or confident doublet) with ≥ 400 transcripts are shown in black over all cells in grey. Each section is rotated/flipped to
  posterior-left, dorsal-up.
  - The script draws all ten sections. The figure crops OR52YW26_1_7 (2.2 mm from the midline) for the HVC and RA rows and OR52YW26_1_4 (1.2 mm)
    for the LMAN rows (identified by eye), and the boxes were drawn by hand.
- **E.** RPCA label transfer (SCTransform, 40 dims) between the finch glutamatergic types and the chicken excitatory types of Zaremba et al. 2025,
  at k.anchor 2–100.
  - Finch side: Glut-SATB2-1 kept, Glut-DACH2-HVCra-Int left out. Chicken side: Ex_Pre_KCNH7 and Ex_Pre_SATB2 removed.
  - The reciprocal score is the geometric mean of the two transfer directions. Each cell is a finch type's best reciprocal score over the chicken types.
  - Same as the plusSATB2_noPre call of `snrna/integration/gg_glutonly_hybrid/rpca/plot_sweep_heatmap_satb2_tests.py`. The k.anchor 20 box and the
    labels were added by hand.
- **F.** The full finch × chicken reciprocal-score matrix at k.anchor 20. Rows are in panel E's order; chicken columns are grouped by name marker
  (DACH2, TSHZ2, CACNA1H, SATB2, KIAA1217, BCL6).
  - Same as `snrna/integration/gg_glutonly_hybrid/rpca/plot_transfer_heatmap_plusSATB2.py` in its fixed-column-order version.
  - The column order is `chicken_col_order()` of `snrna/integration/position_palette.py`; earlier versions of the script ordered the columns by
    clustering.
  - Row labels were dropped in the figure.
- **G.** UMAP of the k.anchor 20 RPCA joint embedding. Left: finch cells by the dissected region of their type, chicken cells grey. Right: chicken
  cells by name-marker family, finch cells grey.
  - Same as `make_finch_chicken_class_figure()` in `snrna/integration/gg_glutonly_hybrid/rpca/plot_rpca_umap_satb2_tests.py`.
  - Titles, legends and labels were redrawn by hand.

## Check against the figure files

Each script was run and its output compared with the file placed in `fig2_snrna_glut_v10.svg`. Raster outputs were compared pixel by pixel, PDFs
rasterized at 100 dpi and their text extracted:

| Panel | Result |
|---|---|
| A | The script recomputes the song-marker removal, the clustered genes (nidopallial: 688 removed, 458 clustered; arcopallial: 620 removed, 210 clustered) and both trees. All equal the notebook's saved run (merge order, tips, heights). With `FIG2_PVCLUST_SWEEP`, both trees are pixel-identical to page 4 of `trees_{nido,arco}_thresholds.pdf`, with identical text. With the seeded bootstrap (default), only the AU labels differ, by up to 5.6 points. |
| B, UMAPs | The per-dataset prep rebuilt from the glutamatergic object equals the notebook's saved prep (same features and cells; data and PCA differ by 0). UMAP coordinates equal the saved ones (34,490 points; maximum difference 1.8 × 10⁻¹⁵). Pixel-identical to `umap_all_by_type_bare`; same text. |
| B, heatmap | All 60 song × non-song scores (forward, reverse, reciprocal) equal `reciprocal_scores.csv` (ALL set, CCA, k.anchor 20). Same text and tile colours as `reciprocal_all_cca_k20_bare`; 1 pixel differs, from a 0.02 pt font-measurement rounding of the page size. |
| C | All ten maps pixel-identical |
| E | Pixel-identical, same text. Values equal `sweep_heatmap_SCT_plusSATB2_noPre_reciprocal.csv` (18 types × 8 k.anchor values). |
| F | Pixel-identical to `transfer_heatmap_plusSATB2_noPre_ka20_kf200_d40_SCT`, same text, including the fixed chicken-column order |
| G | Pixel-identical to `rpca_umap_finchpos_chickenclass_plusSATB2_noPre_ka20_kf200_d40_SCT`, same text |

Each output also matches its panel in `fig2_snrna_glut_v10.svg` by eye. The hand edits are listed in the notes above.

### Sources and scripts that differ from the first manifest

- **E** was made by `plot_sweep_heatmap_satb2_tests.py`, not `plot_sweep_heatmap_glutonly_hybrid_sct.py`. The latter draws the subset without
  Glut-SATB2-1, at k.filter 50.
- **F** uses the fixed chicken-column order of `plot_transfer_heatmap_plusSATB2.py` and `snrna/integration/position_palette.py` as committed with
  this directory; earlier committed versions ordered the columns by clustering.
- **G** was made by `plot_rpca_umap_satb2_tests.py`, not `plot_rpca_umap.py`.
- **B** also contains the "integration partner" heatmap, `song_nonsong_reciprocal_transfer_hybrid/reciprocal_all_cca_k20_bare.pdf`. Its three
  notebooks (`snrna/integration/song_nonsong_{integration,umap,reciprocal_transfer}_hybrid.qmd`) are committed with this directory.
