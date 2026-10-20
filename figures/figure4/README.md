# Figure 4 — standalone panel scripts

One script per panel of Figure 4 (`fig4_scenicplus_v9.svg`). Each regenerates its panel from the outputs of the upstream pipelines, without running
any other panel or notebook. Shared input paths and helpers are in `common.R`, which every R panel script sources; it in turn sources the repository's
`config/paths.R` and `config/figure_theme.R`. Panels are written to `figures/figure4/output/`, which is not tracked; set `FIG4_OUT` to write elsewhere.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| A | GRN activity (+/+ eRegulon target-gene AUC) per TF and cell type | `panel_A_grn_activity_heatmap.R` | `panel_A_positive_gene_auc_heatmap_names.pdf` (+ `_no-names.pdf`) |
| B | GRNs with differential activity between song and matched non-song neurons; MAFB module | `panel_B_grn_contrast_modules.R` | `panel_B_gene_auc_contrast_modules.{svg,pdf,png}`, `..._legend.*`, `.csv` |
| C | MAFB expression density on the UMAP | `panel_C_mafb_umap_density.R` | `panel_C_MAFB_density.pdf` |
| D | *In situ* hybridization (SLC17A6, MAFB, GAD1) in RA and arcopallium | — (microscopy; see below) | — |
| E | Chromatin accessibility at the KCNC1 locus: MAFB cistrome track and KCNC1 expression bars | `panel_E_kcnc1_track_annotations.R` | `panel_E_MAFB_cistrome_regions.bed`, `panel_E_kcnc1_expression_bars_no_points.svg` |
| F | Expression of predicted MAFB targets: RA vs C1H-1 against PVALB-1 vs LGE-1 | `panel_F_mafb_target_scatter.R` | `panel_F_mafb_onoff_scatter_PVALB1_vs_LGE1_absLFC1.{pdf,png,svg}` |
| G | Share of each eRegulon's targets higher in both contrasts, with permutation p | `panel_G_regulon_upper_right_share.R` | `panel_G_quadrant_share_size_n_color_p_flipped.{pdf,png,svg}`, `.csv` |
| H | Simulated MAFB knockout and overexpression (PCA flow field) | `panel_H_mafb_perturbation.py` (+ `panel_H/`) | `panel_H_pca_arrows_kd_MAFB.pdf`, `panel_H_pca_arrows_oe_MAFB.pdf` |

## Running

```bash
Rscript figures/figure4/panel_A_grn_activity_heatmap.R        # R panels: run from anywhere inside the repository (here::here)
python figures/figure4/panel_H_mafb_perturbation.py           # panel H: in the scenicplus conda environment
```

R ≥ 4.3 with the packages in `renv.lock` (Seurat 5, ComplexHeatmap, DESeq2, rhdf5, ArchR/SummarizedExperiment for the marker tables, Nebulosa,
scCustomize, ggrepel, egg, svglite). Panel B also needs Inkscape and Python 3 with lxml (`scripts/crop_svg.py` crops it to the drawing). Run times on
the lab server: A, B 3–5 min (reading the SCENIC+ MuData); C 2 min; E 2 min; F 5 min; G 10 min (DESeq2 + permutations); H ~25 min with its caches.

## Inputs (upstream pipelines)

All paths are set in `common.R` (R panels) or at the top of `panel_H_mafb_perturbation.py`; each can be overridden with the environment variable named there.

| Input | Made by | Used by |
|---|---|---|
| SCENIC+ MuData, hybrid-labeled multiome, all cell types, config37 (`ra-arco-hvc-nc_hybrid/scenicplus/config37/outs/scplusmdata[_slim].h5mu`) | `grn/ra-arco-hvc-nc/scenicplus/` (Snakemake; `config37`) and `assemble_slim_scplusmdata.py` | A, B, E, F, G |
| Hybrid-labeled multiome Seurat object with UMAPs (`multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2`) | `multiome/seurat/reduction_viz/combined_all_umap_hybrid.R` | C, E, F, G |
| ArchR marker tests on the SCENIC+ consensus regions (`ra-arco-hvc-nc_hybrid/archr_consensus/markers_Glut-CACNA1H-RA_VS_Glut-CACNA1H-1.rds`, `markers_GABA-MGE-PVALB-1_VS_GABA-LGE-1.rds`) | `grn/ra-arco-hvc-nc/pycisTopic/make_song_pair_dars_archr.R`, `make_interneuron_dars_archr.R` | G |
| SCENIC+ run on glutamatergic neurons only, config3, with `scplusmdata_perturb.h5mu` (`ra-arco-hvc-nc_glut_hybrid/results/config3/outs/`) | `grn/ra-arco-hvc-nc_glut_hybrid/scenicplus/` (branch `grn-glut-hybrid`) and `panel_H/assemble_perturbation_h5mu.py` | H |

## Notes on individual panels

- **A.** Same computation and drawing as `plot_auc_heatmaps()` in `grn/analysis/scenicplus_hybrid_all_template.qmd` (rendered by `grn/analysis/all.qmd`):
  mean gene-based AUC per cell type, direct and extended eRegulons of a TF pooled, rows scaled 0–1, Ward clustering on Euclidean distances. The only
  difference from the original file is the legend title (the original shows ComplexHeatmap's default name, `matrix_2`).
- **B.** Same as `grn/analysis/scenicplus_gene_auc_cluster_means.R` + `scenicplus_gene_auc_contrast_clusters.R`. A contrast (HVCra vs D2-1, HVCx vs D2-4,
  RA vs C1H-1) is called when the per-bird mean difference in scaled AUC is ≥ 0.15, positive in ≥ 80 % of the paired birds, and the song neuron's
  scaled AUC is ≥ 0.4. This is an effect-size rule, not a test. The legend is written separately (`..._legend`).
- **C.** Same call as `plot_eregulon_umaps()` in the template (Nebulosa density, SCT data, reduction `dims25nn30mindist0.3`).
- **D.** Fluorescent *in situ* hybridization images (`/nest/bcolquit/ish/231001_slc17a6_mafb_gad1/`); acquisition and display settings only, no analysis script.
- **E.** The accessibility tracks (pseudobulk ATAC coverage per cell type at the KCNC1 locus) are a genome-browser export, not scripted; the export
  (`genome_browser/KCNC1/kcnc1_overview.pdf`) was converted to a lean SVG with `scripts/browser_pdf_to_svg.py` (track bars merged, clipped copies
  trimmed). The script regenerates the two code-made parts: the MAFB cistrome track (BED, loaded into the browser) and the KCNC1 expression bars
  (SCT log-normalized expression, mean of per-bird means; RA, C1H-1, PVALB-1, LGE-1).
- **F.** Pseudobulk DESeq2 (library × bird × cell type, ≥ 10 cells, `~ library + group`) for RA vs C1H-1 and PVALB-1 vs LGE-1; MAFB +/+ targets with
  |log2FC| > 1 in at least one contrast. Same as `grn/analysis/scenicplus_pv1_vs_lge1_contrast.R` + `scenicplus_mafb_pvalb_vs_mafb_low_mge.R` (RA vs C1H-1)
  + `scenicplus_mafb_pvalb1_vs_lge1_scatters.R`.
- **G.** Same as `grn/analysis/scenicplus_regulon_quadrant_fractions.R` + `scenicplus_regulon_quadrant_permutation.R` +
  `scenicplus_regulon_quadrant_plot_size_n_color_p_flipped.R` (all with `PV_COMPARATOR=LGE1`). Gene log2FCs are rounded to 3 decimals and the share and
  p-value to 4 significant digits, as in the intermediate tables the figure was drawn from, so the output matches it exactly.
- **H.** `grn/ra-arco-hvc-nc_glut_hybrid/perturbation/perturbation_parallel_template.ipynb` (branch `grn-glut-hybrid`, commit 46dd519, the version that
  wrote the figure files), extracted to a script: the notebook's code cells are unchanged, preceded by the settings of the figure's run (MAFB knockout and
  overexpression, MAFB targets forced into the regressor set, quiver arrows, seed 0, output folder `perturbation_mafb-targets`). Regressors, simulations
  and arrow fields are cached in the output folder and reused.

## Check against the figure files

Each script was run and its output compared with the file placed in `fig4_scenicplus_v9.svg` (text labels extracted from the PDFs; renderings compared pixel by pixel at 100 dpi):

| Panel | Result |
|---|---|
| A | Same TF and cell-type order; only the legend title differs (see above) |
| B | Identical labels; same groups; same size (2.58 × 2.94 in) |
| C | Pixel-identical |
| E | Cistrome BED identical (2,636 regions, same scores); expression-bar SVGs byte-identical |
| F | Same genes, statistics (n = 219, r = 0.12, ρ = 0.08) and labels; only the drawing order of overlapping points differs |
| G | Identical per-eRegulon share, n and p-value (166 rows); pixel-identical |
| H | see the commit message / PR that added this folder |
