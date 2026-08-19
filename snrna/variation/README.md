# Cross-individual variation

`cross_individual_variation.R` asks how much of the snRNA expression variation *within* a cell type
is between birds, using the genotype-derived individual labels from
`snrna/clustering/snrna_souporcell_clustering.qmd`.

Run `snrna/reduction_viz/combined_all_umap_hybrid.R` first: it joins those labels onto
`obj_clustered.qs2` as `individual` / `soup_batch` / `soup_library` and saves the object back, and
this script asserts they are present rather than repeating the join.

## The design, and what it permits

Six birds in two batches that share no individual:

| batch | libraries                | birds       |
|-------|--------------------------|-------------|
| 1     | ra, arco, hvc, nc        | bird1-3     |
| 2     | lman, nr                 | bird4-6     |

A bird term is only estimable within a batch, so nothing is pooled across the two. Bird and library
are crossed within a batch, which is what lets a bird effect be separated from a region effect.

`position` is *not* the library: the preprocessing notebook re-derives position for the arco library
from cluster identity (`Glut-RA -> "ra"`), so ~500 cells sequenced in `ra_run1` carry position
`arco`. Everything here keys on `soup_library`.

## Two analyses, because one design does not cover the question

1. **Crossed decomposition** (`variance_summary_by_celltype.csv`). Per gene per cell type, sequential
   sums of squares for depth, library, individual, residual, with individual entered last. Only
   works for cell types present in >= 2 libraries: **16 of 64 (batch, cell type) combinations**.
   That excludes every song-nucleus type, which is exactly the wrong set to lose.

2. **Split-fold replicates** (`split_fold_summary.csv`). Each bird's cells within one library and
   cell type are split into 3 folds, pseudobulked separately. Folds of one bird differ only by which
   cells were drawn, which measures the sampling floor instead of assuming it, and the test is
   whether birds differ by more than their own folds do. Applies to single-library cell types:
   **28 (library, cell type) combinations**, song nuclei included.

   Three folds and not two because the null is a permutation of bird labels over folds: two folds x
   three birds has only 15 distinct labellings, flooring p at 1/15 regardless of effect size. Three
   folds gives 280 and a floor of 0.0036.

## Results

**Composition.** Cell-type proportions are fairly stable between birds -- median CV 0.22 across cell
types above 1% abundance. The exceptions are large: `Glut-CACNA1H-4` in ra (CV 1.26, 14-fold range),
`Glut-SATB2-1` in nr (40-fold range between birds). Treat those cell types' expression results as
composition results first.

**Expression.** Every one of the 28 split-fold combinations exceeds its sampling floor (whole-profile
permutation test, BH < 0.05), so between-bird differences are real everywhere, not only in a few cell
types. Effect size varies ~3x across cell types, largest in the song-nucleus projection neurons:

| library | cell type | F ratio (all cells) | F ratio (60 cells/bird) |
|---------|-----------|--------------------|--------------------------|
| ra      | Glut-CACNA1H-RA  | 2.93 | 1.58 |
| hvc     | Glut-DACH2-HVCra | 2.35 | 1.21 |
| arco    | Glut-CACNA1H-1   | 2.10 | 1.41 |
| nc      | Glut-DACH2-1     | 1.73 | 1.16 |

The two columns differ because effect size tracks group size (rho = 0.68 between F ratio and cell
count): more cells means a tighter sampling floor, so the same difference reads larger. Subsampling
every bird to 60 cells removes that (rho = 0.24, n.s.) and preserves the ordering (rank rho = 0.72),
with 27 of 28 still significant. **Glut-CACNA1H-RA stays first either way; Glut-DACH2-HVCra's rank
does not survive matching** -- most of its headline value was its 3017 cells.

`split_fold_matched_ranking.pdf` ranks the matched ratio by cell type, coloured by class. The
ordering is broadly by class: glutamatergic projection neurons at the top, GABAergic and
non-neuronal types in the middle, and the two neurogenic populations (`Glut-NSC`, `Glut-NB`) at the
bottom, `Glut-NB` sitting essentially on the sampling floor at 1.02.

**Genes.** 2061 genes are individually significant in at least one combination; only 54 in three or
more. The recurrent ones (DLGAP2, FRMD8, PDE10A, PHACTR1, CACNA2D3, HOMER1, KCNJ6, NTRK2) are the
ones to be most careful with: a gene that separates birds in five different cell types is a property
of the animal or its library -- genotype, or state at sacrifice, `HOMER1` being the obvious tell --
not of the cell type.

## Caveats

The split-fold floor is a *sampling* floor. Folds share a bird, a dissection and a library, so
"between-bird" here still bundles genotype, behavioural state, and dissection variability; nothing in
this design separates those. With three birds per batch, the per-gene tests are weak by construction
-- the whole-profile test is the one with power, and disagreement between them is expected, not a
contradiction.
