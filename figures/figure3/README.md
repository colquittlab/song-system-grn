# Figure 3 — standalone panel scripts

One script per data panel of Figure 3 (`fig3_multiome_overview_v6.svg`). Each regenerates its panel from the outputs of the upstream pipelines, without
running any other panel or notebook. Shared input paths and helpers are in `common.R`, which every panel script sources; it in turn sources the
repository's `config/paths.R` and `config/figure_theme.R`. Panels are written to `figures/figure3/output/`, which is not tracked; set `FIG3_OUT` to write
elsewhere.

| Panel | Content | Script | Output (`output/`) |
|---|---|---|---|
| A | Song-system circuit and the multiome experiment (schematic) | — (drawing) | — |
| B | UMAP of the multiome nuclei (all cell types), colored by cell type | `panel_B_multiome_umap.R` | `panel_B_umap_dims25nn30mindist0.3_cluster_hybrid_no-label.{pdf,png}` (+ `panel_B_cluster_colors.csv`) |
| C (tracks) | Chromatin accessibility at the SIX2 locus, with zoom insets | — (genome-browser export) | — |
| C (bars) | SIX2 expression in the six cell types of the tracks | `panel_C_six2_expression_bars.R` | `panel_C_six2_expression_bars_no_points.svg` (+ `_bars.svg` with per-bird dots, `_labeled.svg`, `.csv`) |
| D | Song vs. non-song DARs: row-scaled accessibility of 1,000 sampled DARs in 8 k-means groups | `panel_D_song_nonsong_dar_heatmap.R` | `panel_D_DAR_fdr0.01_lfc0.6_nk8_scale01.pdf` (+ `panel_D_DAR_rows.csv`, `panel_D_DAR_counts.csv`) |
| E | Motifs enriched in the song-up and song-down DARs (top 25 per pair) | `panel_E_song_nonsong_motif_heatmaps.R` | `panel_E_song-glut_motifs_{up,down}_top25_heatmap.pdf` (+ `.csv`) |
| F | Mean chromVAR motif enrichment per TF family, song minus non-song, compared between song neurons | `panel_F_motif_family_scatter.R` | `panel_F_points_fam_<A>_vs_<B>_norm.pdf` (×3), `panel_F_motif_family_song_minus_nonsong.csv` |
| G | Expression log2 fold-change (song vs. non-song) of the TFs of those families, compared between song neurons | `panel_G_tf_deg_scatter.R` | `panel_G_points_deg_<A>_vs_<B>.pdf` (×3), `panel_G_tf_log2fc.csv` |
| H | log2 fold-change, song vs. non-song, of the MAF, TCF, BHLH and AP-1 family members | `panel_H_tf_family_deg_heatmaps.R` | `panel_H_heatmap_deg_{MAF,TCF,BHLH,AP1}.pdf`, `panel_H_tf_family_log2fc.csv` |

## Running

```bash
Rscript figures/figure3/panel_B_multiome_umap.R                # run from anywhere inside the repository (here::here)
```

R ≥ 4.3 with the packages in `renv.lock` (Seurat 5, scCustomize, ArchR, ComplexHeatmap, circlize, cowplot, qs2). Run times and peak memory on the
lab server:

| Panel | Time | Memory |
|---|---|---|
| B | 40 s | 14 GB (reads the Seurat object) |
| C | 2 min | 14 GB |
| D | 4 min the first time (marker tests and cell-type means from the ArchR project), 1 min cached | 3.4 GB |
| E | 2 min the first time (motif enrichment), 30 s cached | 2.8 GB |
| F | 35 s | 1.4 GB (reads the cached chromVAR deviations) |
| G, H | 35 s each | 13 GB (read the Seurat object for its gene names) |

Values D and E compute from the shared ArchR project (marker tests when their saved file is absent, cell-type mean accessibility, motif enrichment)
are cached in `FIG3_CACHE` (default `output/cache/`). The project itself is only read: ArchR's logs and temporary files go into the cache.

## Inputs (upstream pipelines)

All paths are set in `common.R` (or at the top of panel F); each can be overridden with the environment variable named here. The untracked ones are
outputs of the pipelines below; in a fresh worktree, point the variables at the main checkout's copies.

| Input | Made by | Used by |
|---|---|---|
| Hybrid-labeled multiome Seurat object with the UMAP reductions (`multiome/seurat/reduction_viz/combined_all_umap_hybrid/obj_clustered.qs2`; `FIG3_SEURAT_OBJ`) | `multiome/seurat/reduction_viz/combined_all_umap_hybrid.R` | B, C |
| Shared ArchR project (`ARCHR_PROJ_DIR` in `config/paths.R`): PeakMatrix and MotifHS peak annotation; hybrid labels added in memory by `multiome/archr/hybrid_labels.R` | `multiome/archr/archr_processing.qmd` | D, E |
| ArchR marker tests of the three song vs. non-song pairs (`multiome/archr/peaks_glut_hybrid/peaks_nonsong_glut_se.rds`; `FIG3_DAR_MARKERS`; recomputed into `FIG3_CACHE` when absent) | `multiome/archr/peaks_glut.qmd` ("Find DARs") | D, E |
| Lambert et al. TF table (`data/lambert_tfs/tfs.csv`, tracked) | — | E |
| Glutamatergic-subset ArchR project with the chromVAR deviations (`MotifHSMatrix`, CIS-BP human motifs) of the six song / non-song clusters (`multiome/archr/chromvar_hybrid/proj_glut_sub`, 20 GB; `FIG3_CHROMVAR_PROJ`). Read only; computing it takes ~94 min | `multiome/archr/chromvar.qmd` (chromVAR chunk) | F |
| Pseudobulk song vs. non-song DE table (`<MULTIOME_SEURAT_DIR>/deg/song-surround_deg_glut/marker_genes_cluster_SCT.rds`; `FIG3_DEG_TABLE`); pre-hybrid labels mapped to the hybrid names with `hybrid_cluster_map()` | `multiome/seurat/deg/song-surround_deg_glut.qmd` | G, H |
| Hybrid-labeled multiome object, for its gene names only (`multiome/seurat/label_transfer_hybrid/obj_clustered_hybrid.qs2`; `FIG3_HYBRID_LABEL_OBJ`) | `multiome/seurat/label_transfer_hybrid.qmd` | G, H |

Two of these were lost and restored on 2026-10-23, when the `archr-hybrid` worktree that held them was deleted:
- `proj_glut_sub` was moved from that worktree into the main checkout before the deletion.
- `peaks_nonsong_glut_se.rds` was recomputed by `song_pair_markers()`. Before the deletion it was checked against the saved file: same DAR counts, a
  byte-identical motif enrichment, identical panel D and E output. `peaks_glut.qmd`'s other saved objects (`peaks_nonsong_glut.rds`, `kmeans_8.rds`)
  were not restored; no panel uses them.

## Notes on individual panels

- **B.** Same drawing as the `_no-label` plot of `multiome/seurat/reduction_viz/combined_all_umap_hybrid.R`. It shows the saved reduction
  `dims25nn30mindist0.3` (Harmony-integrated SCT PCA, 25 dimensions, 30 neighbors, min.dist 0.3), which is read from the object, not recomputed.
  16,588 nuclei, 27 cell types (`cluster_hybrid`).
  - Drawn with scCustomize `DimPlot_scCustom` (default point size, shuffled drawing order, seed 1) and the "varibow" palette; axes and legend
    removed; 7 × 7 in. The session theme is cowplot's, as in the source script; it sets the plot margin and so the panel size.
  - Hand edits in Inkscape: the PDF was scaled to 33 % keeping the stroke width, so the dots are about 3× the stroke of a proportional scale; the
    `cluster_hybrid` title was dropped. Added: the cell-type and group labels (Arial 6 pt, in the cluster colors), the three gray dashed regions
    (nidopallial glutamatergic, nidopallial neurogenesis, arcopallial glutamatergic), the UMAP1/UMAP2 arrows and "n = 16,588". Endo (44 nuclei)
    and GABA-MGE-LAMP5 are not labeled.
- **C.** The accessibility tracks (pseudobulk ATAC coverage at chr3:16.505–16.532 Mb, with zoom insets of three regions) are a genome-browser export
  (`hgt_genome_*_six2*.pdf` beside the figure file), not scripted. In v6 the panel's layer was slimmed with
  `scripts/browser_pdf_to_svg.py --within layer2`, as for the Figure 4E tracks: the clipped copies of whole tracks behind the zoom insets were dropped
  and the bars of each track merged into one path, 45,691 → 27 paths (file 17.8 → 7.5 MB). Renderings before and after differ only in anti-aliasing
  at bar edges.
  - Bars: SCT log-normalized SIX2 expression (natural log, log1p of corrected counts), mean of per-bird means (bird = souporcell assignment within a
    library, ≥ 10 cells), one row per track in the track order and colors; same construction as the KCNC1 bars of Figure 4E. They stand to the right
    of the tracks, so they carry no row labels. SIX2 is expressed in the three song neurons only: RA 0.14 (19 % of cells), HVCra 0.07 (7 %), HVCx 0.23
    (29 %); C1H-1, D2-1, D2-4 ≤ 0.001 (≤ 0.15 %).
  - Hand edits: the zoom wedges were drawn by hand.
- **D.** Same computation and drawing as the nk = 8 heatmap chunk of `multiome/archr/peaks_glut.qmd` (commit 112f60f), which wrote
  `peaks_glut_hybrid/DAR_fdr0.01_lfc0.6_nk8_scale01.pdf`.
  - DARs: FDR < 0.01 and |log2FC| ≥ 0.6 in any of the three ArchR marker tests (RA vs. C1H-1, HVCra vs. D2-1, HVCx vs. D2-4; PeakMatrix, Wilcoxon,
    matched on TSS enrichment and log10 fragments), n = 97,420. Value: mean accessibility per cell type (getGroupSE), min–max scaled per row.
  - A random 1,000 regions (seed 42) are split into 8 k-means groups (25 repeats, seed 42) and ordered within groups by Ward clustering on correlation
    distance. The notebook ran under ArchR's L'Ecuyer-CMRG generator (set by `addArchRThreads()`); the script sets it explicitly, otherwise the sample
    and groups differ.
  - Hand edits: only the raster heatmap body (the 8 slices, at the PDF's own scale) was placed, as group `g887492` at the page root. Dendrogram and
    slice numbers removed; column labels rewritten without the `Glut-` prefix and colored by pair; legend redrawn (min–max, "row-scaled
    accessibility"); "DARs" and "n = 97,420" added.
- **E.** Same as the "Plot up and down separately" chunk of `peaks_glut.qmd`.
  - ArchR `peakAnnoEnrichment` (MotifHS, hypergeometric against all peaks) in the DARs with FDR ≤ 0.01 and log2FC ≥ 0.6 (up) or ≤ −0.6 (down) of
    each pair; the 25 motifs named by an Ensembl gene ID are dropped. Rows: the union of each pair's top 25 by −log10 adjusted p (33 up, 31 down).
    Value: fold enrichment, sign flipped for the down set. Rows clustered by ComplexHeatmap defaults; viridis "C" (up) and reversed viridis "E" (down).
  - Hand edits: columns renamed (ra, hvc_1, hvc_2 → RA, HVCra, HVCx for the up set; CACNA1H-1, DACH2-1, DACH2-4 for the down set), dendrograms
    removed, legends redrawn horizontally as "motif enrichment". The down legend reads 0–4 (magnitude) where the PDF reads −4–0. Motif names keep
    their CIS-BP suffixes (e.g. FOSL1_142).
- **F.** Same computation and drawing as the "Remove redundant TFs" section of `multiome/archr/chromvar.qmd` (`points_fam_<A>_vs_<B>_norm.pdf`).
  - chromVAR deviation z-scores on the cells of RA, C1H-1, HVCra, D2-1, HVCx and D2-4; the 84 motifs whose names match the notebook's TF list are
    averaged per cluster and pooled into 29 families (mean of the motif means). Each song neuron has its partner's value subtracted (RA − C1H-1,
    HVCra − D2-1, HVCx − D2-4). Pearson R over the 29 families. Per-cell means; there is no test.
  - Hand edits: "R = " removed, leaving the number; family labels (FOS, NFE2, NR3 [= NR3C], AR, MAF, NFI, MEF2, RFX) added, orange above zero and
    teal below.
- **G.** Same as the `points_deg` chunk of `chromvar.qmd`.
  - `avg_log2FC` of each song neuron vs. its partner for the 80 genes of the multiome object that match the TF-family patterns (`TF_TERMS`); a gene
    missing from a contrast's table counts as 0. Pearson R over the 80 genes.
  - The DE table is `FindMarkers(test.use = "DESeq2")` on `AggregateExpression` (SCT) pseudobulks of cluster × souporcell `assignment`, made under the
    pre-hybrid labels. `assignment` indexes birds *within* a library, so each pseudobulk "sample" pools different birds across libraries.
  - Hand edits: "R = " removed; gene labels (ATOH8, FOSL2, AR, SIX2, NFIB, JUN, MEF2C) added.
- **H.** Same as the `heatmap_deg_<family>` loop of `chromvar.qmd`; the notebook writes all 22 families, the figure stacks MAF, TCF, BHLH and AP1.
  - Values are the panel G log2FCs. Rows clustered with ComplexHeatmap defaults (complete linkage, Euclidean); columns RA, HVCra, HVCx.
  - Shared color scale: PRGn, green at the minimum and purple at the maximum over all 80 TF genes (−2.8 … +3.0), white at 0.
  - The legend titles are ComplexHeatmap's automatic `matrix_<n>` names of the notebook session, kept so the output matches the original files.
  - Hand edits: column titles relabeled ("RA vs. CACNA1H-1", "HVCra vs. DACH2-1", "HVCx vs. DACH2-4"); the four legends replaced by one color bar,
    "log2 fold-change expr. song versus non-song", labeled −3 … 3 (until 2026-10-23 it read −4 … 4, beyond the data range).

## Check against the figure files

Each script was run and its output compared with the upstream file placed in `fig3_multiome_overview_v6.svg` (PDFs rasterized at 100 dpi and
compared pixel by pixel with `compare -metric AE`, text extracted with pdftotext) and with the panel in the figure:

| Panel | Result |
|---|---|
| B | Pixel-identical to `combined_all_umap_hybrid/umap_dims25nn30mindist0.3_cluster_hybrid_no-label.pdf` and `.png`; same text. That PDF is the file placed: its 16,588 dots match group `g66474` of layer "panel b" in color, drawing order and coordinates. |
| C | The bars placed in v6 have the script's widths (4.191, 0.008, 2.057, 0, 6.762, 0.033 mm), title and ticks. |
| D | Same DARs (n = 97,420; per comparison, down/up: RA 10,038/17,913, HVCra 28,141/45,389, HVCx 4,152/20,564) and same k-means groups (97, 153, 86, 101, 104, 138, 142, 179 rows); the 8 heatmap strips are pixel-identical to the ones placed in the figure. Elsewhere the page is shifted 0.34 pt (ComplexHeatmap centers the layout on text widths measured on the notebook session's device). |
| E | Pixel-identical to both PDFs; identical labels |
| F | Pixel-identical to `points_fam_*_norm.pdf` (all three); same text and R (0.93, 0.95, 0.99) as the figure |
| G | Pixel-identical to `points_deg_*.pdf` as re-rendered on 2026-10-23 (all three); same R as the figure (0.60, 0.51, 0.73) |
| H | Pixel-identical to `heatmap_deg_{MAF,TCF,BHLH,AP1}.pdf`; same text (rows: MAF 6, TCF 9, BHLH 5, AP1 8) |

### Differences from the manifest

Found in the manifest of 2026-10-23 (beside the figure file) and corrected there the same day, which now also names each panel's standalone script:

- **B** names `umap_dims25nn30mindist0.3_cluster_hybrid_no-label.png`; the figure uses the vector `.pdf` of the same name. The reduction comes from
  the dims 15–30 sweep of `combined_all_umap_hybrid.R`, which was uncommitted until 2026-10-23 (the committed sweep could not produce it).
- **C** is listed as "schematic"; it is the genome-browser export plus `panel_C_six2_expression_bars.R`.
- **D–H**: sources and scripts are right. The F and G source lists are not in figure order (the figure runs RA vs. HVCra, RA vs. HVCx, HVCra vs. HVCx).
- In v6 the OPC and Oligo labels of panel B were swapped; corrected on 2026-10-23.
