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

## Residency: which library a cell type belongs to

RA sits inside the arcopallium, so the ra dissection also catches surrounding arcopallium:
`Glut-CACNA1H-1` has 891 cells in arco and 277 in ra. Those 277 are real cells but a spillover
fringe, and between-bird differences measured on them are partly differences in how much
arcopallium each dissection happened to include.

A (library, cell type) pair is kept only if the cell type's **density** in that library -- its
proportion within the library, divided by its highest such proportion across the batch's libraries
-- is at least `min_rel_density` (0.5). Density rather than share of the cell type, because share
cannot distinguish a spillover fringe from a genuinely widespread type: an interneuron class spread
over four libraries has ~25% of itself in each, which looks identical to contamination.

The threshold sits in a gap in the data rather than being picked: among the pairs entering the
split-fold analysis, one is at 0.32 (ra / `Glut-CACNA1H-1`), then nothing until 0.61, then a run to
1.0. The three in 0.61-0.75 are `Astro-2`, `GABA-1-1`, `Glut-DACH2-2` at similar density in hvc and
nc -- what a widespread type looks like -- and are kept.

This drops 51 pairs holding >= 20 cells, leaving 29,245 of 33,636 cells. The full table, including
what was dropped and why, is `celltype_library_residency.csv`; raise or lower `min_rel_density` there
if a cell type is being judged wrongly. Note it also drops ubiquitous glia from libraries where they
are much less dense than at their peak (`OPC` outside ra, `Astro-1` in hvc), which is stricter than
the anatomical spillover argument requires.

## Two analyses, because one design does not cover the question

1. **Crossed decomposition** (`variance_summary_by_celltype.csv`). Per gene per cell type, sequential
   sums of squares for depth, library, individual, residual, with individual entered last. Needs the
   cell type in >= 2 libraries, which after the residency filter means >= 2 libraries where it is
   genuinely resident: **9 of 64 (batch, cell type) combinations**, down from 16 before the filter.
   Nothing here reaches significance (smallest BH-adjusted global p = 0.13).

2. **Split-fold replicates** (`split_fold_summary.csv`). Each bird's cells within one library and
   cell type are split into 3 folds, pseudobulked separately. Folds of one bird differ only by which
   cells were drawn, which measures the sampling floor instead of assuming it, and the test is
   whether birds differ by more than their own folds do. Applies to single-library cell types:
   **27 (library, cell type) combinations**, song nuclei included. This is the analysis with power.

   Three folds and not two because the null is a permutation of bird labels over folds: two folds x
   three birds has only 15 distinct labellings, flooring p at 1/15 regardless of effect size. Three
   folds gives 280 and a floor of 0.0036.

## Results

**Composition.** Cell-type proportions are fairly stable between birds -- median CV 0.22 across cell
types above 1% abundance. The exceptions are large: `Glut-CACNA1H-4` in ra (CV 1.26, 14-fold range),
`Glut-SATB2-1` in nr (40-fold range between birds). Treat those cell types' expression results as
composition results first. Composition is computed before the residency filter, so it describes
every cell.

**Expression.** All 27 split-fold combinations exceed their sampling floor (whole-profile permutation
test, BH < 0.05), so between-bird differences are real everywhere, not only in a few cell types.
Effect size varies ~3x, largest in the song-nucleus projection neurons:

| library | cell type | F ratio (all cells) | F ratio (60 cells/bird) |
|---------|-----------|--------------------|--------------------------|
| ra      | Glut-CACNA1H-RA  | 3.00 | 1.54 |
| hvc     | Glut-DACH2-HVCra | 2.45 | 1.21 |
| arco    | Glut-CACNA1H-1   | 2.11 | 1.38 |
| nr      | Glut-DACH2-6     | 1.71 | 1.31 |

The two columns differ because effect size tracks group size (rho = 0.73 between F ratio and cell
count): more cells means a tighter sampling floor, so the same difference reads larger. Subsampling
every bird to 60 cells removes that (rho = 0.25, n.s.) and largely preserves the ordering (rank rho =
0.74), with 25 of 27 still significant -- the two that are not are `Glut-NSC` and `Glut-NB`.
**Glut-CACNA1H-RA stays first either way; Glut-DACH2-HVCra's rank does not survive matching** -- most
of its headline value was its 3017 cells.

The matched ratio is the median over `n_matched_repeats` (5) independent subsamples, not one draw:
a single draw at 20 cells per fold moves the ratio by ~0.08 (max 0.16), which is enough to reorder
the middle of the ranking. Rank correlation between repeats is 0.86 (range 0.69-0.93).

`split_fold_matched_ranking.pdf` ranks the matched ratio by cell type, coloured by class, with the
range over repeats as a line. **Those ranges overlap heavily: neighbouring cell types are not
separable, only the ends of the ranking are.** With that caveat, the top is glutamatergic --
`Glut-CACNA1H-RA` (ra) at 1.54, then `GABA-1-1` (nc) and `Glut-CACNA1H-1` (arco) -- and the bottom
is the two neurogenic populations, `Glut-NB` lowest at 1.08 and not significant.

**Genes.** 2232 genes are individually significant in at least one combination; only 73 in three or
more. The recurrent ones (DSCAM, PDE10A, ROBO2, ABCF1, ASTN2, CACNA2D3, FRMD8, KCNJ6) are the ones to
be most careful with: a gene that separates birds in five different cell types is a property of the
animal or its library -- genotype, or state at sacrifice -- not of the cell type.

## Caveats

The split-fold floor is a *sampling* floor. Folds share a bird, a dissection and a library, so
"between-bird" here still bundles genotype, behavioural state, and dissection variability; nothing in
this design separates those. With three birds per batch, the per-gene tests are weak by construction
-- the whole-profile test is the one with power, and disagreement between them is expected, not a
contradiction.

The neurogenic populations sitting at the bottom of the ranking is not established as biology: less
between-bird difference in newborn neurons is plausible, but a less complex transcriptome dominated
by a shared proliferative programme would compress the same statistic, and both are among the
smallest source populations here.
