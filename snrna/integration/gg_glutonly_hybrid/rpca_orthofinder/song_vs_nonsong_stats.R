## Song vs non-song tests for the RBH and OrthoFinder RPCA sweeps and composites, with the
## same R test as the published figure (song_nonsong_helpers.R: two-sided wilcox.test,
## Cliff's delta), so the P values are directly comparable with
## song_vs_nonsong_plusSATB2_noPre_stats.csv. Reads the per-type best scores and composite
## confidence matrices written by compare_rbh_orthofinder.py.
##
##   Rscript --no-init-file snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/song_vs_nonsong_stats.R
renv_lib <- Sys.getenv("RENV_LIB", "/ssd/brad/rstudio/multiome/song-system-grn/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))
suppressMessages({library(dplyr); library(tidyr); library(purrr)})
here_dir <- "snrna/integration/gg_glutonly_hybrid/rpca_orthofinder"
res <- file.path(here_dir, "results")

SONG <- c("Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA")
cliffs_delta <- function(x, y) mean(sign(outer(x, y, "-")))
test <- function(d) {
  s <- d$score[d$group == "song"]; n <- d$score[d$group == "non-song"]
  mw <- suppressWarnings(wilcox.test(s, n, alternative = "two.sided"))
  tibble(n_song = length(s), n_nonsong = length(n), song_mean = mean(s), nonsong_mean = mean(n),
         cliffs_delta = cliffs_delta(s, n), U = unname(mw$statistic), p_mannwhitney = mw$p.value)
}

best <- read.csv(file.path(res, "best_match_by_type_rbh_vs_orthofinder.csv"), check.names = FALSE) |>
  pivot_longer(c(best_rbh, best_orthofinder), names_to = "orthologs", values_to = "score") |>
  mutate(orthologs = recode(orthologs, best_rbh = "RBH", best_orthofinder = "OrthoFinder"))
rpca <- best |> group_by(method = "RPCA", k_anchor, orthologs) |> group_modify(\(d, k) test(d)) |> ungroup()

comp_best <- function(dir, name) {
  m <- as.matrix(read.csv(file.path(res, dir, "composite_confidence_matrix.csv"), row.names = 1, check.names = FALSE))
  tibble(finch_type = rownames(m), score = apply(m, 1, max), orthologs = name,
         group = if_else(rownames(m) %in% SONG, "song", "non-song"))
}
comp <- bind_rows(comp_best("composite_plusSATB2_noPre_rbh_rebuilt", "RBH"),
                  comp_best("composite_plusSATB2_noPre_orthofinder", "OrthoFinder")) |>
  group_by(method = "Composite", orthologs) |> group_modify(\(d, k) test(d)) |> ungroup()

stats <- bind_rows(rpca, comp) |> arrange(method, k_anchor, desc(orthologs))
write.csv(stats, file.path(res, "song_vs_nonsong_stats_rbh_vs_orthofinder.csv"), row.names = FALSE)
print(as.data.frame(stats |> mutate(across(where(is.double), \(x) signif(x, 3)))))
