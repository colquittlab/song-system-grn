## Summary figure: one sagittal section, shown two ways.
##   A. per-cell transcript number, cells drawn small to limit overlap
##   B. cell-free transcript density, 20 um bins, no segmentation involved
##
## Section OR52YW26_1_4: most cells (181,627) among the sections that need no
## rotation -- it is already posterior-left / dorsal-up in raw coordinates, so
## nothing here depends on the orientation table used by the spatial maps.
##
## The binned density comes from xenium/bin_transcripts.py: nanoparquet
## segfaults on the 2 GB transcript-metadata files, so the 99M transcripts are
## streamed and binned in pyarrow and only the 75k occupied bins reach R.
## Regenerate the input (needs an env with pyarrow; `bgr` has it):
##   ~/miniforge3/envs/bgr/bin/python xenium/bin_transcripts.py \
##     <proseg>/<SECTION>/transcript-metadata.parquet out.csv 20 20
##   gzip -c out.csv > xenium/transcript_density/tx_density_<SECTION>_20um.csv.gz

suppressMessages({ library(tidyverse); library(Matrix); library(here) })
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))
select <- dplyr::select

SECTION <- "OR52YW26_1_4"
BIN_UM  <- 20
out_dir <- here::here("xenium", "transcript_density")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sec_dir <- file.path(path.expand(XENIUM_PROSEG_DIR), SECTION)

## Sequential ramp: ONE hue, light to dark (never a rainbow). Built around
## FIG_PAL's blue so it sits in the project's palette.
FIG_SEQ <- c("#f4f7fd", "#bcd3f0", "#6fa3e0", "#2a78d6", "#14406f")

## --- per-cell transcript number ---------------------------------------------
## proseg's cell-metadata carries no count column, so total assigned
## transcripts per cell is the column sum of the expected-counts matrix. This
## is every cell proseg called, with no RCTD/UMI filter applied.
cm <- nanoparquet::read_parquet(file.path(sec_dir, "cell-metadata.parquet"))
m  <- Matrix::readMM(gzfile(file.path(sec_dir, "expected-counts.csv.gz")))   # cells x genes
stopifnot(nrow(m) == nrow(cm))
cells <- tibble(x = cm$centroid_x, y = cm$centroid_y,
                n_tx = Matrix::rowSums(m))
cat("cells:", nrow(cells), "  median transcripts/cell:", round(median(cells$n_tx), 1),
    "  range:", paste(round(range(cells$n_tx), 1), collapse = "-"), "\n")
## Cells with zero assigned transcripts go to -Inf on a log scale and would be
## silently dropped by the scale. Drop them explicitly and say how many.
n_zero <- sum(cells$n_tx <= 0)
cells <- cells |> filter(n_tx > 0)
cat("zero-transcript cells excluded from the log-scaled panel:", n_zero,
    sprintf("(%.2f%%)\n", 100 * n_zero / (n_zero + nrow(cells))))

## --- cell-free binned density ----------------------------------------------
dens <- read_csv(file.path(out_dir, paste0("tx_density_", SECTION, "_", BIN_UM, "um.csv.gz")),
                 show_col_types = FALSE)
cat("occupied", BIN_UM, "um bins:", nrow(dens),
    "  transcripts:", format(sum(dens$n), big.mark = ","), "\n")

## --- shared spatial panel styling -------------------------------------------
## Axes carry no information on a tissue section, so they come off entirely and
## a 1 mm scale bar goes in instead. theme_fig() still supplies type, colours
## and margins.
## NOTE: blanking the PARENT element is not enough. theme_fig() sets
## axis.title.x/.y and panel.grid.major/.minor explicitly, and those specific
## settings survive a blanked `axis.title` / `panel.grid` -- the stray "x"/"y"
## titles and grid lines came back until each child was blanked by name.
theme_section <- function() {
  theme_fig() +
    theme(axis.title = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.55, "cm"),
          legend.title = element_text(size = FIG_PT_AXIS_TITLE),
          legend.text = element_text(size = FIG_PT_AXIS_TEXT),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(0.1, "cm"),
          legend.position = "right")
}

xr <- range(c(cells$x, dens$x)); yr <- range(c(cells$y, dens$y))
bar_um <- 1000
bar <- tibble(x = xr[1] + 0.03 * diff(xr),
              xend = xr[1] + 0.03 * diff(xr) + bar_um,
              y = yr[1] + 0.05 * diff(yr))
add_scalebar <- function(p) {
  p +
    annotate("segment", x = bar$x, xend = bar$xend, y = bar$y, yend = bar$y,
             colour = FIG_INK_SECONDARY, linewidth = 0.5) +
    ## size in a geom is MILLIMETRES, not points -- fig_pt() or this is ~17 pt
    annotate("text", x = (bar$x + bar$xend) / 2, y = bar$y + 0.035 * diff(yr),
             label = "1 mm", size = fig_pt(5), colour = FIG_INK_SECONDARY)
}

## --- A: per-cell -------------------------------------------------------------
## log10 colour: transcripts/cell spans ~20 to >10,000, so a linear ramp would
## put almost every cell in the bottom fifth of the scale.
## Limits at the 1st-99th percentile with out-of-range squished: on the full
## range a handful of very bright cells compress everything else into one flat
## mid-blue, which is what the first version looked like.
lim_a <- quantile(cells$n_tx, c(0.01, 0.99))
pa <- ggplot(cells |> arrange(n_tx), aes(x, y, colour = n_tx)) +
  geom_point(size = 0.035, shape = 16) +
  scale_colour_gradientn(colours = FIG_SEQ, trans = "log10",
                         limits = lim_a, oob = scales::squish,
                         breaks = c(100, 300, 900),
                         labels = c("100", "300", "900"),
                         name = "transcripts\nper cell") +
  coord_equal() +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  theme_section()
pa <- add_scalebar(pa)
fig_save(pa, file.path(out_dir, paste0("per_cell_transcripts_", SECTION)), width = 3.4, height = 2.1)

## --- B: cell-free density ----------------------------------------------------
## geom_tile, not geom_raster: only occupied bins are stored, so the spacing in
## the data is uneven and geom_raster would shift pixels to fake a regular grid.
lim_b <- unname(quantile(dens$n[dens$n >= 100], c(0.05, 0.95)))
cat("density scale limits (tissue bins, 5th-95th pct):", round(lim_b), "\n")
pb <- ggplot(dens, aes(x, y, fill = n)) +
  geom_tile(width = BIN_UM, height = BIN_UM) +
  ## Linear, not log, and scaled to the TISSUE bins only. Within tissue the
  ## range is just ~4x (5th-95th pct 600-2250, median 1393), so a log scale
  ## spends almost the whole ramp on near-empty edge bins and renders the
  ## section as one flat block of dark blue -- which is what the first version
  ## did. Limits from bins with >=100 transcripts, out-of-range squished.
  scale_fill_gradientn(colours = FIG_SEQ,
                       limits = lim_b, oob = scales::squish,
                       breaks = c(800, 1400, 2000),
                       labels = c("800", "1,400", "2,000"),
                       name = paste0("transcripts\nper ", BIN_UM, " \u00b5m bin")) +
  coord_equal() +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  theme_section()
pb <- add_scalebar(pb)
fig_save(pb, file.path(out_dir, paste0("transcript_density_", SECTION, "_", BIN_UM, "um")),
         width = 3.4, height = 2.1)

## --- B2: coarser bins ---------------------------------------------------
## 20 um bins are below what the printed panel can resolve, so their Poisson
## noise shows up as speckle rather than signal. 40 um quadruples the count per
## bin and halves the relative noise. Aggregated from the 20 um grid by summing
## 2x2 blocks -- exact, and avoids re-streaming 99M transcripts.
BIN2 <- BIN_UM * 2
## Recover the INTEGER bin index with round() before pairing, rather than
## doing floor() arithmetic on the centre coordinates. The centres are written
## to one decimal, so (x - min(x)) / BIN2 carries float error and floor() lands
## a bin off at particular columns -- that put 5-9 source bins into some
## targets instead of 4 and drew two dark seam lines across the section.
dens2 <- dens |>
  mutate(ix = round((x - min(x)) / BIN_UM), iy = round((y - min(y)) / BIN_UM),
         x2 = (ix %/% 2) * BIN2 + min(x) + BIN_UM / 2,
         y2 = (iy %/% 2) * BIN2 + min(y) + BIN_UM / 2) |>
  group_by(x2, y2) |> summarize(n = sum(n), n_src = n(), .groups = "drop") |>
  rename(x = x2, y = y2)
stopifnot("more than 4 source bins per 40 um bin" = max(dens2$n_src) <= 4)
lim_b2 <- unname(quantile(dens2$n[dens2$n >= 400], c(0.05, 0.95)))
cat("occupied", BIN2, "um bins:", nrow(dens2), "  scale limits:", round(lim_b2), "\n")

pb2 <- ggplot(dens2, aes(x, y, fill = n)) +
  geom_tile(width = BIN2, height = BIN2) +
  scale_fill_gradientn(colours = FIG_SEQ, limits = lim_b2, oob = scales::squish,
                       breaks = c(3000, 5500, 8000),
                       labels = c("3,000", "5,500", "8,000"),
                       name = paste0("transcripts\nper ", BIN2, " \u00b5m bin")) +
  coord_equal() +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  theme_section()
pb2 <- add_scalebar(pb2)
fig_save(pb2, file.path(out_dir, paste0("transcript_density_", SECTION, "_", BIN2, "um")),
         width = 3.4, height = 2.1)

## small tracked summary so the numbers are checkable without re-running
write_csv(tibble(section_id = SECTION,
                 n_cells = nrow(cells) + n_zero,
                 median_tx_per_cell = median(cells$n_tx),
                 mean_tx_per_cell = round(mean(cells$n_tx), 1),
                 n_zero_tx_cells = n_zero,
                 q05 = quantile(cells$n_tx, .05), q95 = quantile(cells$n_tx, .95),
                 n_transcripts = sum(dens$n),
                 bin_um = BIN_UM, n_occupied_bins = nrow(dens)),
          file.path(out_dir, "summary_stats.csv"))
cat("DONE\n")
