library(Seurat)
library(tidyverse)
library(qs2)
library(cowplot)
library(scCustomize)
library(ComplexHeatmap)
library(circlize)
library(here)
theme_set(theme_cowplot())

source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "scRNA.R"))
source(file.path(Sys.getenv("COLQUITTLAB_UTILS", "/opt/colquittlab/utils"), "R", "common_aesthetics.R"))

## How much of the expression variation within a cell type is between *birds*?
##
## Input is the object written by snrna/reduction_viz/combined_all_umap_hybrid.R, which carries the
## genotype-derived individual labels from snrna/clustering/snrna_souporcell_clustering.qmd. Run
## that script first; this one asserts the labels are present rather than re-joining them, so there
## is one place where cell -> bird is defined.
##
## What the design allows, and what it does not:
##
##   * Six birds in two batches that share no individual -- bird1-3 in ra/arco/hvc/nc, bird4-6 in
##     lman/nr. A bird term is only estimable *within* a batch, so every model below is fit per
##     batch and the two are never pooled. batch1 (4 libraries x 3 birds) is the informative one;
##     batch2 (2 x 3) is reported alongside it as a replication check, not as an equal partner.
##   * Bird and library are crossed, not nested: each library holds all three of its batch's birds,
##     each bird appears in every library of its batch. That is what makes the question answerable
##     at all -- a bird effect can be separated from a region effect -- but it also means the
##     library term has to be in the model. Region differences are far larger than bird differences
##     here, so a model without it would hand most of the region variance to whichever bird happened
##     to contribute more cells.
##   * `position` is not the library. The preprocessing notebook re-derived position for the arco
##     library from cluster identity (Glut-RA -> "ra"), so ~500 ra_run1 cells carry position "arco".
##     Library identity (`soup_library`) is what was sequenced, and is what the model uses.
##
## The estimand is a variance decomposition per gene per cell type: sequencing depth first (a
## nuisance the pseudobulk normalisation does not fully remove), then library, then individual, then
## residual, as sequential sums of squares in that order. Individual is entered *last*, so it is
## credited only with variation that region cannot already explain -- the conservative direction.
## Its null is not zero: with 12 pseudobulks and 2 df, a bird term absorbs variance by construction,
## so the same decomposition is run on bird labels permuted within library, and the permuted
## distribution is what the observed values are read against.

# Directories -------------------------------------------------------------

viz_dir = here::here("snrna/reduction_viz/combined_all_umap_hybrid")
data_fname = file.path(viz_dir, "obj_clustered.qs2")

script_name = "cross_individual_variation"
out_dir = here::here("snrna/variation", script_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res_to_use = "celltype_hybrid"

# Parameters --------------------------------------------------------------

min_cells_per_pseudobulk = 20L   # below this a pseudobulk is a handful of cells' noise
min_pseudobulks_per_ct = 6L      # a cell type needs enough samples for the model to have df
min_libraries_per_ct = 2L
min_individuals_per_ct = 3L
min_df_residual = 2L             # residual df the bird term is tested against
min_cpm = 1                      # gene detection floor, in CPM
min_frac_detected = 0.5          # ... in this fraction of the cell type's pseudobulks
n_perm = 50L                     # permutations of bird label within library
seed = 10

set.seed(seed)

# Load --------------------------------------------------------------------

obj = qs_read(data_fname, nthreads = 8)

stopifnot(
  "no `individual` metadata -- run snrna/reduction_viz/combined_all_umap_hybrid.R first" =
    all(c("individual", "soup_batch", "soup_library") %in% colnames(obj@meta.data))
)

md = obj@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(celltype = .data[[res_to_use]],
         individual = as.character(individual),
         soup_batch = as.character(soup_batch)) %>%
  select(cell, celltype, individual, soup_batch, soup_library, position)

## The design, from the data.
design_tbl = md %>% count(soup_batch, soup_library, individual)
write_csv(design_tbl, file.path(out_dir, "design_cells_per_library_individual.csv"))
print(design_tbl %>% pivot_wider(names_from = soup_library, values_from = n, values_fill = 0))

# Cell-type composition across individuals --------------------------------

## Before expression: does each bird even contribute the same *mix* of cells? This is the cheaper
## question and the one that most often explains an apparent expression difference -- a cell type
## whose abundance swings between birds is also the one whose pseudobulk is most easily contaminated
## by its neighbours. Proportions are taken within (library, individual), because that is the unit
## that was dissected and captured.
comp = md %>%
  count(soup_batch, soup_library, individual, celltype, name = "n_cells") %>%
  complete(nesting(soup_batch, soup_library, individual), celltype, fill = list(n_cells = 0)) %>%
  group_by(soup_batch, soup_library, individual) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

comp_summary = comp %>%
  group_by(soup_batch, soup_library, celltype) %>%
  summarise(mean_prop = mean(prop),
            sd_prop = sd(prop),
            cv_prop = sd(prop) / mean(prop),
            min_prop = min(prop),
            max_prop = max(prop),
            fold_range = (max(prop) + 1e-6) / (min(prop) + 1e-6),
            n_individuals = n(),
            .groups = "drop")
write_csv(comp_summary, file.path(out_dir, "composition_variation_by_celltype.csv"))
write_csv(comp, file.path(out_dir, "composition_proportions.csv"))

## Ordered by how variable the cell type is between birds, so the top of the panel is the answer.
comp_order = comp_summary %>%
  group_by(celltype) %>%
  summarise(cv = median(cv_prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(cv))

gg = comp %>%
  mutate(celltype = factor(celltype, levels = comp_order$celltype)) %>%
  ggplot(aes(celltype, prop, colour = individual)) +
  geom_point(size = 1.2) +
  facet_wrap(~ soup_library, ncol = 2, scales = "free_y") +
  scale_y_log10() +
  labs(x = NULL, y = "proportion of library's cells (log10)",
       title = "Cell-type composition, by bird, within library") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "composition_by_individual.pdf"), gg, base_height = 10, base_asp = 1.3)
save_plot(file.path(out_dir, "composition_by_individual.png"), gg, base_height = 10, base_asp = 1.3)

# Pseudobulk --------------------------------------------------------------

## Counts, not the SCT residuals: the decomposition below is over library-size-normalised counts,
## and SCT's Pearson residuals are already a variance-stabilised quantity whose scale depends on the
## model fit per cell, which is not what should be summed across cells.
## Group by an opaque id rather than by a pasted "celltype|library|bird" string: AggregateExpression
## rewrites characters in group names (underscores become dashes, among others), and cell type names
## here contain both dashes and digits, so parsing the identity back out of the column name is a
## silent-mislabelling risk. The id maps back through pb_meta.
pb_meta = md %>%
  count(celltype, soup_batch, soup_library, individual, name = "n_cells") %>%
  mutate(pb_id = sprintf("pb%04d", row_number()))

pb_group = md %>%
  left_join(pb_meta %>% select(celltype, soup_library, individual, pb_id),
            by = c("celltype", "soup_library", "individual"))
stopifnot(identical(pb_group$cell, colnames(obj)))
obj$pb_group = pb_group$pb_id

pb_counts_fname = file.path(out_dir, "pseudobulk_counts.qs2")
if (!file.exists(pb_counts_fname)) {
  pb = AggregateExpression(obj, assays = "RNA", group.by = "pb_group")$RNA
  qs_save(pb, pb_counts_fname, nthreads = 8)
} else {
  pb = qs_read(pb_counts_fname, nthreads = 8)
}

stopifnot("pseudobulk columns do not match the group table" =
            setequal(colnames(pb), pb_meta$pb_id))
pb_meta = pb_meta %>% mutate(col = pb_id)
pb = pb[, pb_meta$col]

pb_meta = pb_meta %>% mutate(total_counts = colSums(pb))
write_csv(pb_meta %>% select(-col), file.path(out_dir, "pseudobulk_metadata.csv"))

## Drop thin pseudobulks. A cell type present in one bird at 5 cells and another at 500 is a
## composition result (above), not an expression one, and keeping the 5-cell sample here would make
## it look like the latter.
keep_pb = pb_meta$n_cells >= min_cells_per_pseudobulk
message(sprintf("pseudobulks: %s total, %s with >= %s cells",
                nrow(pb_meta), sum(keep_pb), min_cells_per_pseudobulk))
pb_meta = pb_meta[keep_pb, ]
pb = pb[, pb_meta$col]

## CPM on the pseudobulk's own total, then log. The depth term in the model below is what catches
## the part of depth that this does not remove.
cpm = t(t(pb) / colSums(pb)) * 1e6
## Dense: ~15.5k genes x a few hundred pseudobulks is tens of MB, and every consumer below (qr,
## cor, apply) wants a base matrix anyway.
logcpm = as.matrix(log2(cpm + 1))

# Variance decomposition --------------------------------------------------

## Sequential (type I) sums of squares from nested QR fits, all genes at once: with one design
## matrix per cell type, qr.resid() over the whole gene matrix costs one decomposition rather than
## one lm() per gene, which is the difference between minutes and hours over ~40 cell types x 2
## batches x 50 permutations.
rss = function(qr_obj, Y) colSums(qr.resid(qr_obj, Y)^2)

decompose = function(Y, meta, use_depth = TRUE) {
  ## Y: genes x samples (rows are genes). Transposed here because qr.resid wants samples in rows,
  ## and densified: qr.resid has no sparse method, and a cell type's slice is a few hundred columns.
  Yt = t(as.matrix(Y))
  meta$depth = scale(log10(meta$total_counts))[, 1]

  q0 = qr(model.matrix(~ 1, meta))
  ## The depth term costs a df, and batch2 has only six pseudobulks per cell type against a model
  ## that already spends four. Where it does not fit, it is dropped rather than the cell type being
  ## discarded -- CPM has already removed the first-order depth effect, so this is the term most
  ## affordable to lose, and `used_depth` records which cell types lost it.
  q1 = if (use_depth) qr(model.matrix(~ depth, meta)) else q0
  q2 = if (use_depth) qr(model.matrix(~ depth + soup_library, meta)) else qr(model.matrix(~ soup_library, meta))
  q3 = if (use_depth) qr(model.matrix(~ depth + soup_library + individual, meta)) else
    qr(model.matrix(~ soup_library + individual, meta))

  ss_total = rss(q0, Yt)
  r1 = rss(q1, Yt); r2 = rss(q2, Yt); r3 = rss(q3, Yt)

  df_res = nrow(meta) - q3$rank
  df_ind = q3$rank - q2$rank

  tibble(
    gene = colnames(Yt),
    var_depth = (ss_total - r1) / ss_total,
    var_library = (r1 - r2) / ss_total,
    var_individual = (r2 - r3) / ss_total,
    var_residual = r3 / ss_total,
    ss_total = ss_total,
    ## `if`, not `ifelse`: the condition is scalar (a property of the design, not of the gene), and
    ## ifelse() would return a length-1 result and silently recycle one gene's F over all of them.
    F_individual = if (df_res > 0 && df_ind > 0) ((r2 - r3) / df_ind) / (r3 / df_res) else NA_real_,
    df_individual = df_ind,
    df_residual = df_res
  )
}

## Which (batch, cell type) combinations the design can actually support. Stated in df rather than
## in a flat sample-count cut, because the two batches have different geometry: batch1 offers up to
## 12 pseudobulks per cell type against a rank-6 model, batch2 only 6 against rank-4, so one
## threshold on n_pb would either exclude batch2 entirely or admit rank-deficient batch1 fits.
ct_design = pb_meta %>%
  group_by(soup_batch, celltype) %>%
  summarise(n_pb = n(),
            n_lib = n_distinct(soup_library),
            n_ind = n_distinct(individual),
            n_cells = sum(n_cells),
            .groups = "drop") %>%
  mutate(rank_reduced = 1 + (n_lib - 1) + (n_ind - 1),
         df_residual = n_pb - rank_reduced,
         use_depth = n_pb - rank_reduced - 1 >= min_df_residual,
         testable = n_lib >= min_libraries_per_ct &
           n_ind >= min_individuals_per_ct &
           n_pb >= min_pseudobulks_per_ct &
           df_residual >= min_df_residual)
write_csv(ct_design, file.path(out_dir, "celltype_design.csv"))
message(sprintf("testable cell type x batch combinations: %s of %s",
                sum(ct_design$testable), nrow(ct_design)))

run_one = function(batch, celltype_cur, use_depth) {
  idx = which(pb_meta$soup_batch == batch & pb_meta$celltype == celltype_cur)
  meta_cur = pb_meta[idx, ] %>%
    mutate(soup_library = factor(soup_library), individual = factor(individual))
  Y = logcpm[, idx, drop = FALSE]

  ## Detection filter inside the cell type: a gene off in this cell type contributes only zeros,
  ## and its variance fractions are then noise divided by noise.
  detected = rowMeans(cpm[, idx, drop = FALSE] >= min_cpm) >= min_frac_detected
  Y = Y[detected, , drop = FALSE]
  if (nrow(Y) < 100) return(NULL)

  obs = decompose(Y, meta_cur, use_depth = use_depth) %>%
    mutate(soup_batch = batch, celltype = celltype_cur, n_genes = nrow(Y), used_depth = use_depth)

  ## Null: bird labels shuffled within library, which preserves the library term and the depth
  ## profile and destroys only the bird assignment. Genes are pooled across permutations; the
  ## comparison used below is between the observed and null *distributions* of var_individual.
  perm = map_dfr(seq_len(n_perm), function(p) {
    meta_p = meta_cur %>%
      group_by(soup_library) %>%
      mutate(individual = sample(individual)) %>%
      ungroup()
    decompose(Y, meta_p, use_depth = use_depth) %>%
      transmute(gene, perm = p, var_individual = var_individual, F_individual = F_individual)
  })

  list(obs = obs, perm = perm)
}

todo = ct_design %>% filter(testable)
results = pmap(list(todo$soup_batch, todo$celltype, todo$use_depth), function(b, ct, ud) {
  message(sprintf("  %s / %s (depth term: %s)", b, ct, ud))
  run_one(b, ct, use_depth = ud)
})
names(results) = paste(todo$soup_batch, todo$celltype, sep = "|")
results = compact(results)

obs_all = map_dfr(results, "obs")
perm_all = map_dfr(results, "perm", .id = "key") %>%
  separate(key, into = c("soup_batch", "celltype"), sep = "\\|")

## Empirical p per gene: how often a permuted bird assignment reaches this gene's F, using the
## cell type's whole permuted F distribution (all genes x all permutations) as the null. Per-gene
## permutation alone would give at best 1/n_perm resolution; pooling over genes within the cell type
## buys resolution at the cost of assuming genes share a null shape, which is the usual trade and is
## the reason the summary below leads with the distributional comparison rather than gene counts.
null_F = perm_all %>% select(soup_batch, celltype, F_individual) %>% filter(is.finite(F_individual))
obs_all = obs_all %>%
  group_by(soup_batch, celltype) %>%
  group_modify(function(d, key) {
    nulls = sort(null_F$F_individual[null_F$soup_batch == key$soup_batch & null_F$celltype == key$celltype])
    d$p_perm = if (length(nulls) == 0) NA_real_ else
      (length(nulls) - findInterval(d$F_individual, nulls) + 1) / (length(nulls) + 1)
    d
  }) %>%
  ungroup() %>%
  group_by(soup_batch, celltype) %>%
  mutate(padj_perm = p.adjust(p_perm, method = "BH")) %>%
  ungroup()

write_csv(obs_all, file.path(out_dir, "variance_partition_by_gene.csv.gz"))

# Summaries ---------------------------------------------------------------

perm_summary = perm_all %>%
  group_by(soup_batch, celltype) %>%
  summarise(median_var_individual_null = median(var_individual, na.rm = TRUE),
            q95_var_individual_null = quantile(var_individual, 0.95, na.rm = TRUE),
            .groups = "drop")

## The same whole-cell-type test as in the split-fold section: each permutation gives one median
## variance fraction, and the observed median is read against those 50 numbers. This is the test
## with power at n = 3 birds; the per-gene p-values are the ones that will mostly come back empty.
perm_medians = perm_all %>%
  group_by(soup_batch, celltype, perm) %>%
  summarise(median_var_individual = median(var_individual, na.rm = TRUE), .groups = "drop")

global_test = obs_all %>%
  group_by(soup_batch, celltype) %>%
  summarise(obs_median = median(var_individual, na.rm = TRUE), .groups = "drop") %>%
  left_join(perm_medians, by = c("soup_batch", "celltype"), relationship = "one-to-many") %>%
  group_by(soup_batch, celltype, obs_median) %>%
  summarise(p_global = (sum(median_var_individual >= dplyr::first(obs_median)) + 1) / (n() + 1),
            median_of_perm_medians = median(median_var_individual),
            .groups = "drop") %>%
  select(-obs_median)

ct_summary = obs_all %>%
  group_by(soup_batch, celltype) %>%
  summarise(n_genes = dplyr::first(n_genes),
            median_var_individual = median(var_individual, na.rm = TRUE),
            median_var_library = median(var_library, na.rm = TRUE),
            median_var_depth = median(var_depth, na.rm = TRUE),
            median_var_residual = median(var_residual, na.rm = TRUE),
            n_sig = sum(padj_perm < 0.05, na.rm = TRUE),
            frac_sig = mean(padj_perm < 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(perm_summary, by = c("soup_batch", "celltype")) %>%
  left_join(global_test, by = c("soup_batch", "celltype")) %>%
  mutate(excess_var_individual = median_var_individual - median_var_individual_null) %>%
  left_join(ct_design, by = c("soup_batch", "celltype")) %>%
  mutate(padj_global = p.adjust(p_global, method = "BH")) %>%
  arrange(desc(excess_var_individual))
write_csv(ct_summary, file.path(out_dir, "variance_summary_by_celltype.csv"))
print(ct_summary %>% select(soup_batch, celltype, n_cells, median_var_individual,
                            median_var_individual_null, excess_var_individual,
                            p_global, padj_global, frac_sig), n = 30)

## Genes that come up in many cell types. A bird difference confined to one cell type is a
## cell-type-specific result; one that appears in most cell types of a bird is a property of the
## animal or its library -- genotype (cis-eQTL-like), sex, or state at sacrifice -- and reading it
## as cell biology is the main way this analysis goes wrong.
recurrent = obs_all %>%
  filter(padj_perm < 0.05) %>%
  count(soup_batch, gene, name = "n_celltypes") %>%
  left_join(obs_all %>% group_by(soup_batch, gene) %>%
              summarise(median_var_individual = median(var_individual), .groups = "drop"),
            by = c("soup_batch", "gene")) %>%
  arrange(desc(n_celltypes))
write_csv(recurrent, file.path(out_dir, "recurrent_individual_variable_genes.csv"))
print(head(recurrent, 25))

top_genes = obs_all %>%
  group_by(soup_batch, celltype) %>%
  slice_max(var_individual, n = 25) %>%
  ungroup() %>%
  arrange(soup_batch, celltype, desc(var_individual))
write_csv(top_genes, file.path(out_dir, "top_individual_variable_genes.csv"))

# Split-fold replicates: between-bird vs sampling noise ---------------------

## The decomposition above can only see cell types that appear in two or more libraries, which
## excludes exactly the cell types this project is about: Glut-DACH2-HVCra, -HVCx, -LMANsh,
## Glut-CACNA1H-RA and the rest are one-library types by definition, so a model with a library term
## has nothing to estimate and they drop out (see celltype_design.csv -- Glut-DACH2-HVCra is 3017
## cells and n_lib = 1). They are also the types where a bird difference would matter most.
##
## So for those, ask the question a different way. Within one library and one cell type, split each
## bird's cells at random into `n_folds` folds and pseudobulk each fold separately. Folds of the
## same bird differ only by which cells were drawn -- that is the sampling floor, measured rather
## than assumed -- and the question becomes whether birds differ by more than their own folds do. It
## is a one-way design with real replicates, needs no library term, and applies to every cell type
## with enough cells, single-library or not.
##
## Three folds rather than two, because the null is a permutation of bird labels over the folds and
## the *number of distinct labellings* is what floors the p-value. Two folds x three birds gives 15
## distinct ways to group six samples into three pairs, so no result can be reported below p = 1/15
## no matter how large the effect -- the test would be structurally unable to reach 0.05. Three
## folds gives 9!/(3!^3 3!) = 280 groupings and a floor of 0.0036, and takes the residual df from 3
## to 6, which is also what makes the per-gene F usable. The cost is a higher cell requirement
## (n_folds x min cells per fold per bird), which drops the thinnest groups.
##
## What it cannot do: folds share a bird, a dissection and a library, so this floor is the
## *sampling* floor only. A difference that clears it is between-bird, but "between-bird" here still
## bundles genotype, state at sacrifice, and anything that varied in the dissection of that animal.

n_folds = 3L
min_cells_per_fold = 20L
min_cells_split = n_folds * min_cells_per_fold  # per bird per library
min_birds_split = 3L
n_perm_split = 500L

set.seed(seed)
split_md = md %>%
  group_by(celltype, soup_library, individual) %>%
  filter(n() >= min_cells_split) %>%
  ## rep_len then shuffle, so folds are within one cell of equal size rather than multinomial --
  ## an unequal fold is a noisier pseudobulk and would widen the floor it is supposed to measure.
  mutate(fold = sample(rep_len(as.character(seq_len(n_folds)), n()))) %>%
  ungroup()

split_meta = split_md %>%
  count(celltype, soup_batch, soup_library, individual, fold, name = "n_cells") %>%
  mutate(sp_id = sprintf("sp%04d", row_number()))

split_groups = split_md %>%
  left_join(split_meta %>% select(celltype, soup_library, individual, fold, sp_id),
            by = c("celltype", "soup_library", "individual", "fold"))

## Cells not in a qualifying group get NA and are dropped from this aggregation.
sp_vec = setNames(rep(NA_character_, ncol(obj)), colnames(obj))
sp_vec[split_groups$cell] = split_groups$sp_id
obj$split_group = unname(sp_vec)

sp_counts_fname = file.path(out_dir, "split_fold_counts.qs2")
if (!file.exists(sp_counts_fname)) {
  obj_sp = subset(obj, cells = colnames(obj)[!is.na(obj$split_group)])
  sp = AggregateExpression(obj_sp, assays = "RNA", group.by = "split_group")$RNA
  qs_save(sp, sp_counts_fname, nthreads = 8)
  rm(obj_sp)
} else {
  sp = qs_read(sp_counts_fname, nthreads = 8)
}
stopifnot("split-fold columns do not match the group table" = setequal(colnames(sp), split_meta$sp_id))
sp = sp[, split_meta$sp_id]
sp_cpm = t(t(sp) / colSums(sp)) * 1e6
sp_logcpm = as.matrix(log2(sp_cpm + 1))

split_design = split_meta %>%
  group_by(soup_batch, soup_library, celltype) %>%
  summarise(n_birds = n_distinct(individual),
            n_folds_total = n(),
            n_cells = sum(n_cells),
            .groups = "drop") %>%
  ## Every fold of every bird, or the within-bird floor is estimated from a different set of birds
  ## than the between-bird term, and the permutation's group sizes stop being exchangeable.
  filter(n_folds_total == n_folds * n_birds, n_birds >= min_birds_split)
write_csv(split_design, file.path(out_dir, "split_fold_design.csv"))
message(sprintf("split-fold combinations: %s", nrow(split_design)))

split_one = function(lib, ct, meta_tbl = split_meta, cpm_mat = sp_cpm, logcpm_mat = sp_logcpm) {
  idx = which(meta_tbl$soup_library == lib & meta_tbl$celltype == ct)
  meta_cur = meta_tbl[idx, ] %>% mutate(individual = factor(individual))
  Y = logcpm_mat[, idx, drop = FALSE]
  detected = rowMeans(cpm_mat[, idx, drop = FALSE] >= min_cpm) >= min_frac_detected
  Y = Y[detected, , drop = FALSE]
  if (nrow(Y) < 100) return(NULL)

  Yt = t(Y)
  q0 = qr(model.matrix(~ 1, meta_cur))
  qb = qr(model.matrix(~ individual, meta_cur))
  df_b = qb$rank - q0$rank
  df_w = nrow(meta_cur) - qb$rank

  ss0 = rss(q0, Yt); ssw = rss(qb, Yt)
  Fobs = ((ss0 - ssw) / df_b) / (ssw / df_w)

  ## Null: bird labels permuted over the folds. Group sizes are preserved (every bird keeps
  ## n_folds samples), so the permuted statistic asks exactly "how far apart would arbitrary groups
  ## of fold-samples be".
  null_F = map(seq_len(n_perm_split), function(p) {
    meta_p = meta_cur %>% mutate(individual = sample(individual))
    qp = qr(model.matrix(~ individual, meta_p))
    sswp = rss(qp, Yt)
    ((ss0 - sswp) / df_b) / (sswp / df_w)
  })
  null_vec = sort(unlist(null_F))

  p_perm = (length(null_vec) - findInterval(Fobs, null_vec) + 1) / (length(null_vec) + 1)

  ## Per-gene significance is the wrong resolution for this design: with 2 and 3 df the permuted F
  ## has a tail that swallows almost any single gene, so a cell type whose whole profile shifts with
  ## bird still shows ~0 significant genes. The statistic that has power here is the *median* F over
  ## genes, tested against the median F of each permutation -- one number per permutation, asking
  ## whether the cell type as a whole distinguishes birds. Both are reported; they answer different
  ## questions, and disagreement between them ("no gene significant, profile clearly bird-dependent")
  ## is the expected outcome at n = 3, not a contradiction.
  perm_medians = map_dbl(null_F, median)
  obs_median_F = median(Fobs)
  p_global = (sum(perm_medians >= obs_median_F) + 1) / (length(perm_medians) + 1)

  tibble(
    soup_library = lib,
    celltype = ct,
    gene = colnames(Yt),
    var_between_birds = (ss0 - ssw) / ss0,
    F_between = Fobs,
    p_perm = p_perm,
    padj_perm = p.adjust(p_perm, method = "BH"),
    median_null_F = median(null_vec),
    obs_median_F = obs_median_F,
    median_of_perm_medians = median(perm_medians),
    p_global = p_global,
    n_birds = nlevels(meta_cur$individual)
  )
}

split_res = map2_dfr(split_design$soup_library, split_design$celltype, function(lib, ct) {
  message(sprintf("  split-fold: %s / %s", lib, ct))
  split_one(lib, ct)
})

write_csv(split_res, file.path(out_dir, "split_fold_variance_by_gene.csv.gz"))

split_summary = split_res %>%
  group_by(soup_library, celltype) %>%
  summarise(n_genes = n(),
            n_birds = dplyr::first(n_birds),
            median_var_between_birds = median(var_between_birds),
            median_F = median(F_between),
            median_null_F = dplyr::first(median_null_F),
            median_of_perm_medians = dplyr::first(median_of_perm_medians),
            F_ratio = dplyr::first(obs_median_F) / dplyr::first(median_of_perm_medians),
            p_global = dplyr::first(p_global),
            n_sig = sum(padj_perm < 0.05),
            frac_sig = mean(padj_perm < 0.05),
            .groups = "drop") %>%
  left_join(split_design %>% select(soup_batch, soup_library, celltype, n_cells),
            by = c("soup_library", "celltype")) %>%
  mutate(padj_global = p.adjust(p_global, method = "BH")) %>%
  arrange(desc(F_ratio))
write_csv(split_summary, file.path(out_dir, "split_fold_summary.csv"))
print(split_summary, n = 40)

## The effect size has a built-in confound worth measuring rather than warning about: more cells per
## fold means a tighter sampling floor, so the same true between-bird difference reads as a larger
## F_ratio in a large group than in a small one. If F_ratio tracks group size, the ordering of cell
## types is partly an ordering of how many cells each has. Reported here so the reader can discount
## it; the plot below shows the same thing per point.
size_conf = cor.test(log10(split_summary$n_cells), split_summary$F_ratio, method = "spearman")
message(sprintf("F_ratio vs log10(n_cells): rho = %.2f, p = %.3g", size_conf$estimate, size_conf$p.value))
write_csv(tibble(rho = size_conf$estimate, p = size_conf$p.value, n = nrow(split_summary)),
          file.path(out_dir, "split_fold_size_confound.csv"))

split_top = split_res %>%
  group_by(soup_library, celltype) %>%
  slice_max(F_between, n = 25) %>%
  ungroup() %>%
  arrange(soup_library, celltype, desc(F_between))
write_csv(split_top, file.path(out_dir, "split_fold_top_genes.csv"))

## Genes that separate birds in many (library, cell type) combinations at once. Same caution as
## above, and more force here: with the sampling floor measured, a gene that is bird-variable
## everywhere is a property of the animal or the library, not of a cell type.
split_recurrent = split_res %>%
  filter(padj_perm < 0.05) %>%
  count(gene, name = "n_combinations") %>%
  arrange(desc(n_combinations))
write_csv(split_recurrent, file.path(out_dir, "split_fold_recurrent_genes.csv"))
print(head(split_recurrent, 25))

gg = split_summary %>%
  mutate(label = paste(soup_library, celltype),
         significant = padj_global < 0.05) %>%
  ggplot(aes(n_cells, F_ratio, colour = soup_batch, shape = significant)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = label), size = 2, max.overlaps = 15) +
  scale_x_log10() +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = "cells in the (library, cell type) group (log10)",
       y = "median F between birds / median F under permuted bird labels",
       title = "Between-bird differences against a split-fold sampling floor",
       subtitle = "filled = whole-profile permutation test, BH < 0.05; dashed line = sampling floor") +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "split_fold_frac_significant.pdf"), gg, base_height = 6, base_asp = 1.5)
save_plot(file.path(out_dir, "split_fold_frac_significant.png"), gg, base_height = 6, base_asp = 1.5)

gg = split_res %>%
  mutate(label = paste(soup_library, celltype)) %>%
  ggplot(aes(fct_reorder(label, var_between_birds), var_between_birds)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2, fill = "grey90") +
  coord_flip() +
  labs(x = NULL, y = "fraction of gene variance between birds (rest is split-fold sampling)") +
  theme(axis.text.y = element_text(size = 5),
        plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "split_fold_variance_by_group.pdf"), gg, base_height = 10, base_asp = 0.8)
save_plot(file.path(out_dir, "split_fold_variance_by_group.png"), gg, base_height = 10, base_asp = 0.8)

# Split-fold, cell-matched -------------------------------------------------

## The confound measured just above (rho ~ 0.7 between F_ratio and group size) has a fix: give every
## group the same number of cells. Each bird is subsampled to exactly n_folds x min_cells_per_fold
## cells before folding, so the sampling floor is identical everywhere and F_ratio can be compared
## between cell types rather than only against 1. The price is that every group is now as thin as
## the thinnest one -- 20 cells per fold -- so this version is *less* powerful in absolute terms and
## is here for the ranking, not for the significance calls. Read the unmatched version for "is there
## a bird effect in this cell type", this one for "which cell types have the largest one".

set.seed(seed + 1)
matched_md = md %>%
  group_by(celltype, soup_library, individual) %>%
  filter(n() >= min_cells_split) %>%
  slice_sample(n = min_cells_split) %>%
  mutate(fold = sample(rep_len(as.character(seq_len(n_folds)), n()))) %>%
  ungroup()

matched_meta = matched_md %>%
  count(celltype, soup_batch, soup_library, individual, fold, name = "n_cells") %>%
  mutate(sp_id = sprintf("ms%04d", row_number()))

matched_groups = matched_md %>%
  left_join(matched_meta %>% select(celltype, soup_library, individual, fold, sp_id),
            by = c("celltype", "soup_library", "individual", "fold"))

ms_vec = setNames(rep(NA_character_, ncol(obj)), colnames(obj))
ms_vec[matched_groups$cell] = matched_groups$sp_id
obj$matched_group = unname(ms_vec)

ms_counts_fname = file.path(out_dir, "split_fold_matched_counts.qs2")
if (!file.exists(ms_counts_fname)) {
  obj_ms = subset(obj, cells = colnames(obj)[!is.na(obj$matched_group)])
  ms = AggregateExpression(obj_ms, assays = "RNA", group.by = "matched_group")$RNA
  qs_save(ms, ms_counts_fname, nthreads = 8)
  rm(obj_ms)
} else {
  ms = qs_read(ms_counts_fname, nthreads = 8)
}
stopifnot("matched columns do not match the group table" = setequal(colnames(ms), matched_meta$sp_id))
ms = ms[, matched_meta$sp_id]
ms_cpm = t(t(ms) / colSums(ms)) * 1e6
ms_logcpm = as.matrix(log2(ms_cpm + 1))

matched_design = matched_meta %>%
  group_by(soup_batch, soup_library, celltype) %>%
  summarise(n_birds = n_distinct(individual), n_folds_total = n(), n_cells = sum(n_cells),
            .groups = "drop") %>%
  filter(n_folds_total == n_folds * n_birds, n_birds >= min_birds_split)

matched_res = map2_dfr(matched_design$soup_library, matched_design$celltype, function(lib, ct) {
  message(sprintf("  matched: %s / %s", lib, ct))
  split_one(lib, ct, meta_tbl = matched_meta, cpm_mat = ms_cpm, logcpm_mat = ms_logcpm)
})

matched_summary = matched_res %>%
  group_by(soup_library, celltype) %>%
  summarise(n_genes = n(),
            median_var_between_birds = median(var_between_birds),
            median_F = median(F_between),
            F_ratio = dplyr::first(obs_median_F) / dplyr::first(median_of_perm_medians),
            p_global = dplyr::first(p_global),
            n_sig = sum(padj_perm < 0.05),
            .groups = "drop") %>%
  left_join(matched_design %>% select(soup_batch, soup_library, celltype), by = c("soup_library", "celltype")) %>%
  mutate(padj_global = p.adjust(p_global, method = "BH")) %>%
  arrange(desc(F_ratio))
write_csv(matched_summary, file.path(out_dir, "split_fold_matched_summary.csv"))
print(matched_summary, n = 40)

matched_conf = cor.test(log10(split_summary$n_cells[match(paste(matched_summary$soup_library, matched_summary$celltype),
                                                          paste(split_summary$soup_library, split_summary$celltype))]),
                        matched_summary$F_ratio, method = "spearman")
message(sprintf("matched F_ratio vs log10(original n_cells): rho = %.2f, p = %.3g",
                matched_conf$estimate, matched_conf$p.value))

## Does equalising cell number change the ranking? If the two orderings agree, the size confound was
## inflating the numbers without reordering the cell types, and the unmatched ranking can be read as
## biology; if they disagree, it could not.
rank_cmp = split_summary %>%
  select(soup_library, celltype, n_cells, F_ratio_unmatched = F_ratio) %>%
  inner_join(matched_summary %>% select(soup_library, celltype, F_ratio_matched = F_ratio),
             by = c("soup_library", "celltype"))
rank_rho = cor(rank_cmp$F_ratio_unmatched, rank_cmp$F_ratio_matched, method = "spearman")
message(sprintf("unmatched vs matched F_ratio rank correlation: rho = %.2f (n = %s)",
                rank_rho, nrow(rank_cmp)))
write_csv(rank_cmp, file.path(out_dir, "split_fold_matched_vs_unmatched.csv"))

gg = rank_cmp %>%
  mutate(label = paste(soup_library, celltype)) %>%
  ggplot(aes(F_ratio_unmatched, F_ratio_matched)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(aes(size = n_cells), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2, max.overlaps = 15) +
  labs(x = "F ratio, all cells", y = sprintf("F ratio, %s cells per bird", min_cells_split),
       title = sprintf("Cell-number matching, rank rho = %.2f", rank_rho)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "split_fold_matched_vs_unmatched.pdf"), gg, base_height = 6, base_asp = 1.4)
save_plot(file.path(out_dir, "split_fold_matched_vs_unmatched.png"), gg, base_height = 6, base_asp = 1.4)

## Ranking of the matched F ratio, which is the version comparable across cell types.
##
## Class comes from the explicit precursor set in snrna/naming/hybrid_division_naming.qmd
## (`precursor_overrides`), not from a name pattern: "-Pre" would misfile Glut-DACH2-HVCra-Pre, which
## is a subtype of the HVCra projection neurons rather than a precursor.
neurogenic_types = c("Glut-NSC", "Glut-NB", "Glut-Im", "GABA-Im")
cell_class = function(x) {
  case_when(x %in% neurogenic_types ~ "Neurogenic",
            str_starts(x, "Glut-") ~ "Glut",
            str_starts(x, "GABA-") ~ "GABA",
            TRUE ~ "Non-neuronal")
}
class_colors = c(Glut = "#c0392b", GABA = "#2c7fb8", Neurogenic = "#7b3294", `Non-neuronal` = "#7f8c8d")

gg = matched_summary %>%
  mutate(label = sprintf("%s (%s)", celltype, sub("_run1$", "", soup_library)),
         class = factor(cell_class(celltype), levels = names(class_colors))) %>%
  ggplot(aes(F_ratio, fct_reorder(label, F_ratio), colour = class)) +
  ## The sampling floor: birds no further apart than an arbitrary regrouping of their own folds.
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = class_colors, name = NULL, drop = FALSE) +
  ## Two rows: four classes in one row is wider than the panel here and the last label is cut off.
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = sprintf("between-bird F ratio (%s cells per bird)", min_cells_split), y = NULL) +
  ## The axis title and legend are sized down with the tick labels rather than left at the cowplot
  ## default: at this width (base_asp 0.6) a 14 pt title and legend run past the panel and are
  ## clipped, which is not visible in the object, only in the saved file.
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 7),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size = 6),
        legend.key.size = unit(8, "pt"),
        plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "split_fold_matched_ranking.pdf"), gg, base_height = 6, base_asp = 0.6)
save_plot(file.path(out_dir, "split_fold_matched_ranking.png"), gg, base_height = 6, base_asp = 0.6)

# Figures -----------------------------------------------------------------

## Observed against its own permuted null, per cell type. A cell type sitting on the diagonal has no
## bird effect beyond what shuffled labels produce, whatever its absolute value.
gg = ct_summary %>%
  ggplot(aes(median_var_individual_null, median_var_individual, colour = soup_batch)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(aes(size = n_cells)) +
  ggrepel::geom_text_repel(aes(label = celltype), size = 2, max.overlaps = 20) +
  scale_size_continuous(range = c(1, 5)) +
  labs(x = "median variance fraction, bird labels permuted within library",
       y = "median variance fraction attributed to bird",
       title = "Between-bird expression variation vs its permutation null") +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "variance_individual_vs_null.pdf"), gg, base_height = 6, base_asp = 1.4)
save_plot(file.path(out_dir, "variance_individual_vs_null.png"), gg, base_height = 6, base_asp = 1.4)

## The full decomposition, cell types ordered by the excess over null.
ct_levels = ct_summary %>% arrange(excess_var_individual) %>% pull(celltype) %>% unique()
gg = obs_all %>%
  select(soup_batch, celltype, var_depth, var_library, var_individual, var_residual) %>%
  pivot_longer(starts_with("var_"), names_to = "term", values_to = "frac") %>%
  mutate(term = factor(sub("^var_", "", term), levels = c("depth", "library", "individual", "residual")),
         celltype = factor(celltype, levels = ct_levels)) %>%
  ggplot(aes(celltype, frac, fill = term)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
  facet_wrap(~ soup_batch, ncol = 1) +
  coord_flip() +
  labs(x = NULL, y = "fraction of gene variance across pseudobulks") +
  theme(axis.text.y = element_text(size = 5),
        plot.background = element_rect(fill = "white", colour = NA))
save_plot(file.path(out_dir, "variance_partition_by_celltype.pdf"), gg, base_height = 12, base_asp = 0.8)
save_plot(file.path(out_dir, "variance_partition_by_celltype.png"), gg, base_height = 12, base_asp = 0.8)

## Sample-sample correlation over the pseudobulks of the most abundant cell types, annotated by bird
## and library. This is the sanity check behind the numbers: if the decomposition says bird matters
## little, these should cluster by library and not by bird.
cor_celltypes = ct_summary %>%
  group_by(celltype) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  slice_max(n_cells, n = 6) %>%
  pull(celltype)

for (ct in cor_celltypes) {
  idx = which(pb_meta$celltype == ct)
  if (length(idx) < 4) next
  Y = logcpm[, idx, drop = FALSE]
  v = apply(Y, 1, var)
  Y = Y[order(v, decreasing = TRUE)[1:min(2000, sum(v > 0))], , drop = FALSE]
  cm = cor(Y, method = "spearman")
  ann = HeatmapAnnotation(
    bird = pb_meta$individual[idx],
    library = pb_meta$soup_library[idx],
    col = list(bird = setNames(DiscretePalette_scCustomize(6, palette = "varibow"),
                               sprintf("bird%s", 1:6)))
  )
  hm = Heatmap(cm, name = "rho", top_annotation = ann,
               col = colorRamp2(quantile(cm, c(0.01, 0.5, 1)), c("#2166ac", "#f7f7f7", "#b2182b")),
               show_row_names = FALSE, show_column_names = FALSE,
               column_title = ct)
  pdf(file.path(out_dir, sprintf("pseudobulk_correlation_%s.pdf", gsub("[^A-Za-z0-9]", "-", ct))),
      width = 6, height = 5)
  draw(hm)
  dev.off()
}

message("done: ", out_dir)
