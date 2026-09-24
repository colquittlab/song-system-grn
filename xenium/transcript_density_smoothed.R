## Summary figures for one section: per-cell transcript number, and cell-free
## transcript density at several bin sizes, Gaussian-smoothed.
##
## Companion to transcript_density.R (section OR52YW26_1_4, single-hue ramp,
## unsmoothed). This one is parameterised and uses viridis.
##
## NOTE ON COLOUR: viridis is multi-hue, which departs from the sequential
## single-hue rule in CLAUDE.md / figure_theme.R. Used here because it was
## asked for explicitly; FIG_SEQ is still available via PALETTE = "seq".
##
## Regenerate the binned input (needs an env with pyarrow; `bgr` has it):
##   ~/micromamba/envs/bgr/bin/python xenium/bin_transcripts.py \
##     <proseg>/<SECTION>/transcript-metadata.parquet out.csv 10 20
##   gzip -c out.csv > xenium/transcript_density/tx_density_<SECTION>_10um.csv.gz

suppressMessages({ library(tidyverse); library(Matrix); library(here) })
source(here::here("config/paths.R"))
source(here::here("config/figure_theme.R"))
select <- dplyr::select

SECTION  <- "OR52YW26_1_7"
BASE_BIN <- 10                 # bin size of the CSV on disk
BINS     <- c(10, 20, 40)      # 20 and 40 are aggregated from BASE_BIN
SIGMA_UM <- 30                 # Gaussian sd for the smoothed panels, in microns
PALETTE  <- "viridis"
out_dir  <- here::here("xenium", "transcript_density")
sec_dir  <- file.path(path.expand(XENIUM_PROSEG_DIR), SECTION)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal_cols <- if (PALETTE == "viridis") viridisLite::viridis(256) else FIG_SEQ

## Orientation: unlike OR52YW26_1_4 (which needs none), this section is ccw90
## in the table the spatial maps use, so both the cells and the transcript bins
## are rotated into the shared posterior-left / dorsal-up frame.
## ccw90 is (x, y) -> (-y, x).
rot_ccw90 <- function(df, xcol = "x", ycol = "y") {
  x <- df[[xcol]]; y <- df[[ycol]]
  df[[xcol]] <- -y; df[[ycol]] <- x
  df
}

## --- per-cell transcript number ---------------------------------------------
cm <- nanoparquet::read_parquet(file.path(sec_dir, "cell-metadata.parquet"))
m  <- Matrix::readMM(gzfile(file.path(sec_dir, "expected-counts.csv.gz")))
stopifnot(nrow(m) == nrow(cm))
cells <- tibble(x = cm$centroid_x, y = cm$centroid_y, n_tx = Matrix::rowSums(m)) |>
  rot_ccw90()
n_zero <- sum(cells$n_tx <= 0)
cells <- cells |> filter(n_tx > 0)
cat("cells:", nrow(cells) + n_zero, " median tx/cell:", round(median(cells$n_tx), 1),
    " zero-tx excluded:", n_zero, "\n")

## --- binned density, aggregation and smoothing ------------------------------
dens0 <- read_csv(file.path(out_dir, paste0("tx_density_", SECTION, "_", BASE_BIN, "um.csv.gz")),
                  show_col_types = FALSE)
cat("occupied", BASE_BIN, "um bins:", nrow(dens0),
    " transcripts:", format(sum(dens0$n), big.mark = ","), "\n")

## Aggregate via the INTEGER bin index (round() first). Doing floor() arithmetic
## on the centre coordinates -- which are written to one decimal -- carries
## float error and silently merges 5-9 source bins into some targets, which
## draws dark seam lines across the section.
aggregate_bins <- function(df, from_bin, to_bin) {
  if (to_bin == from_bin) return(df |> mutate(n_src = 1L))
  f <- to_bin / from_bin
  stopifnot(f == round(f))
  x0 <- min(df$x); y0 <- min(df$y)
  out <- df |>
    mutate(ix = round((x - x0) / from_bin) %/% f,
           iy = round((y - y0) / from_bin) %/% f) |>
    group_by(ix, iy) |>
    summarize(n = sum(n), n_src = n(), .groups = "drop") |>
    mutate(x = ix * to_bin + x0 + from_bin / 2,
           y = iy * to_bin + y0 + from_bin / 2) |>
    select(x, y, n, n_src)
  stopifnot("aggregation merged more bins than it should" = max(out$n_src) <= f^2)
  out
}

## Normalised 2D Gaussian smoothing on the bin grid.
## Normalised, i.e. the value grid and an occupancy mask are smoothed with the
## same kernel and divided: a plain blur would drag tissue values toward zero
## wherever the kernel overlaps empty space, darkening every edge and the
## ventricle rim. Bins whose smoothed occupancy stays below MASK_MIN are left
## NA so the blur cannot invent tissue outside the section.
MASK_MIN <- 0.15
smooth_grid <- function(df, bin, sigma_um) {
  x0 <- min(df$x); y0 <- min(df$y)
  ix <- round((df$x - x0) / bin); iy <- round((df$y - y0) / bin)
  nx <- max(ix) + 1L; ny <- max(iy) + 1L
  V <- matrix(0, nx, ny); W <- matrix(0, nx, ny)
  V[cbind(ix + 1L, iy + 1L)] <- df$n
  W[cbind(ix + 1L, iy + 1L)] <- 1
  s <- sigma_um / bin
  half <- max(1L, ceiling(3 * s))
  k <- dnorm(seq(-half, half), sd = s); k <- k / sum(k)
  conv_rows <- function(M) {                       # convolve along columns of M
    pad <- (length(k) - 1L) / 2L
    apply(M, 2, function(v) {
      z <- stats::filter(c(rep(0, pad), v, rep(0, pad)), k, sides = 2)
      as.numeric(z)[(pad + 1L):(pad + length(v))]
    })
  }
  Vs <- conv_rows(t(conv_rows(t(V))))              # separable: x then y
  Ws <- conv_rows(t(conv_rows(t(W))))
  stopifnot(dim(Vs) == c(nx, ny))
  val <- ifelse(Ws > MASK_MIN, Vs / pmax(Ws, 1e-9), NA_real_)
  ## Ordering matters and is easy to get backwards: expand_grid() varies its
  ## LAST variable fastest, so with (iy, ix) the ix index moves fastest -- which
  ## matches as.vector() of an [nx, ny] matrix (column-major, first index
  ## fastest). Using t(val) here instead scrambles the grid into stripes.
  stopifnot(length(val) == nx * ny)
  expand_grid(iy = seq_len(ny) - 1L, ix = seq_len(nx) - 1L) |>
    mutate(x = ix * bin + x0, y = iy * bin + y0, n = as.vector(val)) |>
    filter(!is.na(n)) |> select(x, y, n)
}

theme_section <- function() {
  theme_fig() +
    theme(axis.title = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.55, "cm"),
          legend.title = element_text(size = FIG_PT_AXIS_TITLE),
          legend.text = element_text(size = FIG_PT_AXIS_TEXT),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(0.1, "cm"))
}

## scale bar, from the rotated cell extent so every panel shares it
xr <- range(cells$x); yr <- range(cells$y); bar_um <- 1000
bx <- xr[1] + 0.03 * diff(xr); by <- yr[1] + 0.04 * diff(yr)
add_scalebar <- function(p) p +
  annotate("segment", x = bx, xend = bx + bar_um, y = by, yend = by,
           colour = FIG_INK_SECONDARY, linewidth = 0.5) +
  annotate("text", x = bx + bar_um / 2, y = by + 0.028 * diff(yr),
           label = "1 mm", size = fig_pt(5), colour = FIG_INK_SECONDARY)

## The ccw90 rotation puts this section in the same posterior-left/dorsal-up
## frame as the others, which is LANDSCAPE (7898 x 5369 um after rotation) --
## the raw section is portrait, so it is the rotated extent that sets the panel.
PW <- 3.4; PH <- 2.1

## --- A: per-cell -------------------------------------------------------------
lim_a <- quantile(cells$n_tx, c(0.01, 0.99))
pa <- ggplot(cells |> arrange(n_tx), aes(x, y, colour = n_tx)) +
  geom_point(size = 0.035, shape = 16) +
  scale_colour_gradientn(colours = pal_cols, trans = "log10",
                         limits = lim_a, oob = scales::squish,
                         breaks = c(100, 300, 900), labels = c("100", "300", "900"),
                         name = "transcripts\nper cell") +
  coord_equal() +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  theme_section()
fig_save(add_scalebar(pa), file.path(out_dir, paste0("per_cell_transcripts_", SECTION)),
         width = PW, height = PH)

## --- B: density, raw and smoothed, at each bin size --------------------------
stats_rows <- list()
for (b in BINS) {
  d <- aggregate_bins(dens0, BASE_BIN, b) |> rot_ccw90()
  floor_n <- 100 * (b / 20)^2                       # "in tissue" scales with bin area
  lim <- unname(quantile(d$n[d$n >= floor_n], c(0.05, 0.95)))
  brk <- signif(seq(lim[1], lim[2], length.out = 3), 2)

  p_raw <- ggplot(d, aes(x, y, fill = n)) +
    geom_tile(width = b, height = b) +
    scale_fill_gradientn(colours = pal_cols, limits = lim, oob = scales::squish,
                         breaks = brk, labels = scales::comma(brk),
                         name = paste0("transcripts\nper ", b, " µm bin")) +
    coord_equal() +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
    theme_section()
  fig_save(add_scalebar(p_raw),
           file.path(out_dir, paste0("transcript_density_", SECTION, "_", b, "um")),
           width = PW, height = PH)

  ds <- smooth_grid(d, b, SIGMA_UM)
  lim_s <- unname(quantile(ds$n[ds$n >= floor_n], c(0.02, 0.98)))
  brk_s <- signif(seq(lim_s[1], lim_s[2], length.out = 3), 2)
  p_sm <- ggplot(ds, aes(x, y, fill = n)) +
    geom_raster() +
    scale_fill_gradientn(colours = pal_cols, limits = lim_s, oob = scales::squish,
                         breaks = brk_s, labels = scales::comma(brk_s),
                         name = paste0("transcripts\nper ", b, " µm bin\n(",
                                       SIGMA_UM, " µm smooth)")) +
    coord_equal() +
    scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
    theme_section()
  fig_save(add_scalebar(p_sm),
           file.path(out_dir, paste0("transcript_density_", SECTION, "_", b, "um_smooth", SIGMA_UM)),
           width = PW, height = PH)

  cat("bin", b, "um: raw bins", nrow(d), " smoothed bins", nrow(ds),
      " limits", round(lim), "\n")
  stats_rows[[as.character(b)]] <- tibble(section_id = SECTION, bin_um = b,
    n_bins_raw = nrow(d), n_bins_smoothed = nrow(ds),
    median_n = median(d$n[d$n >= floor_n]), lo = lim[1], hi = lim[2],
    sigma_um = SIGMA_UM)
}

write_csv(bind_rows(stats_rows) |>
            mutate(n_cells = nrow(cells) + n_zero,
                   median_tx_per_cell = round(median(cells$n_tx), 1),
                   n_transcripts = sum(dens0$n)),
          file.path(out_dir, paste0("summary_stats_", SECTION, ".csv")))
cat("DONE\n")
