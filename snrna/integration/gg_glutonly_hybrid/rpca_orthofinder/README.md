# RPCA on OrthoFinder orthologs (finch Glut × chicken Ex, plusSATB2_noPre)

Re-run of the finch × chicken RPCA label-transfer sweep (Figure 2E–G, the GlutIntegration supplementary figure) with the gene set
changed from reciprocal best BLASTP hits to the OrthoFinder 2.5.4 one-to-one orthologs that GSI uses, so that RPCA and GSI share
one ortholog source. The original analysis in `../rpca/` and `../composite/` is untouched; everything here is written to this folder.

Everything except the genes is held fixed:
- the same cells (the median-subsampled lists of the original run), the same labels;
- the same, unchanged `finch-integration-toolkit/R/cca_transfer.R`, with SCTransform per species, RPCA, dims 1:40, k.filter 200 (nominal under SCT),
  max.features 200, k.score 30, and k.anchor 2, 5, 10, 20, 30, 50, 75, 100.

| Step | Script | Output |
|---|---|---|
| Inputs: 8,726 OrthoFinder genes × 12,848 finch / 16,214 chicken cells | `prepare_inputs.R` | `data/seurat/*_{query,ref}_*` (`.mtx` not tracked) |
| RPCA sweep (8 runs, both directions) | `run_rpca_sweep.sh` | `results/cca/` |
| RBH vs OrthoFinder comparison; composite rebuilt with each | `compare_rbh_orthofinder.py` | `results/*.csv`, `results/composite_plusSATB2_noPre_{rbh_rebuilt,orthofinder}/` |
| Song vs non-song tests (same R test as the figure) | `song_vs_nonsong_stats.R` | `results/song_vs_nonsong_stats_rbh_vs_orthofinder.csv` |

```bash
Rscript --no-init-file snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/prepare_inputs.R     # ~5 min, ~40 GB
bash snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/run_rpca_sweep.sh                      # ~15 min, 8 runs in parallel
FINCH_INTEGRATION_TOOLKIT=<toolkit> COMPOSITE_REFERENCE=<current composite_confidence_matrix.csv> \
  python snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/compare_rbh_orthofinder.py
Rscript --no-init-file snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/song_vs_nonsong_stats.R
```

## Inputs and checks

- **Finch counts:** `RNA` assay counts (CellBender-corrected) of `obj_hybrid_labels.qs2`, all 17,029 genes. The local `obj_hybrid_labels.h5ad` is
  not usable: it is the later MAGIC/cellxgene export (15,523 genes, X not raw counts).
- **Chicken counts:** `RNA` counts of Zaremba et al. 2025's `Gg_adult_snRNA_seq_srt.rds`, rounded, with Ensembl 87 IDs translated to gene names
  as `build_gg_adult_h5ad.py` does. The ID → name map (`data/ens87_gg_geneid_to_name.tsv`) was rebuilt from Ensembl's release-87
  Gallus_gallus-5.0 GTF; the cluster's copy is not available here.
- **Check:** on GSI's subsampled cells (identical to the RPCA cells), the genes detected in ≥ 5 cells equal GSI's SCTransform gene list exactly,
  for finch (8,219/8,219, `RAW` counts give 8,526) and chicken (8,555/8,555). This confirms the counts, the ID translation and the cell sets.
- **Reference composite:** the rebuilt RBH composite equals the current `composite_plusSATB2_noPre/composite_confidence_matrix.csv`
  exactly (max difference 0) and reproduces the manuscript's numbers. The version committed in git at the time is older (19 rows,
  Glut-DACH2-HVCra-Int still scored, older SAMap input), hence `COMPOSITE_REFERENCE`.
- Seurat 5.5.0 (local renv) here; the original sweep ran on the cluster (comments cite 5.5.1 source).

## Results

Song (5) vs non-song (13) types, best reciprocal score per finch type; two-sided Mann–Whitney U.

| k.anchor | RBH song / non-song | RBH δ, P | OrthoFinder song / non-song | OrthoFinder δ, P | Same top chicken match | Best-score Spearman ρ |
|---|---|---|---|---|---|---|
| 2 | 0.001 / 0.223 | −0.91, 0.0014 | 0.006 / 0.222 | −0.85, 0.0042 | 13/18 | 0.97 |
| 5 | 0.009 / 0.311 | −0.94, 0.0009 | 0.026 / 0.309 | −0.86, 0.0030 | 14/18 | 0.97 |
| 10 | 0.015 / 0.333 | −0.91, 0.0016 | 0.072 / 0.328 | −0.78, 0.0098 | 12/18 | 0.94 |
| 20 | 0.039 / 0.349 | −0.85, 0.0044 | 0.113 / 0.341 | −0.63, 0.046 | 15/18 | 0.91 |
| 30 | 0.055 / 0.347 | −0.85, 0.0044 | 0.141 / 0.339 | −0.63, 0.046 | 16/18 | 0.91 |
| 50 | 0.075 / 0.342 | −0.85, 0.0044 | 0.152 / 0.336 | −0.63, 0.046 | 15/18 | 0.88 |
| 75 | 0.083 / 0.336 | −0.85, 0.0044 | 0.156 / 0.330 | −0.60, 0.059 | 15/18 | 0.90 |
| 100 | 0.094 / 0.330 | −0.82, 0.0068 | 0.147 / 0.326 | −0.60, 0.059 | 16/18 | 0.90 |

Composite (GSI + SAMap + RPCA at k.anchor 20): RBH 0.40 vs 0.65 (δ −0.72, P = 0.019); OrthoFinder 0.43 vs 0.65 (δ −0.69, P = 0.026).
All three methods agree on a reciprocal top-2 match for 1/5 song vs 9/13 non-song types with either ortholog set; at least two agree
for 2/5 vs 11/13 (RBH) and 1/5 vs 11/13 (OrthoFinder).

Over all finch × chicken pairs, the two reciprocal-score matrices agree closely (Pearson r 0.95–0.98). The change is concentrated in one
type: **Glut-DACH2-HVCx**, whose best score at k.anchor 20 rises from 0.069 (Ex_KIAA1217_BCL6) to 0.326 (Ex_KIAA1217), about the non-song
average. The other song types move little (HVCra 0.048 → 0.068, LMANco 0.019 → 0.099, LMANsh 0.027 → 0.053, RA 0.034 → 0.020), and no
non-song type changes by more than 0.1.
