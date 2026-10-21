#!/usr/bin/env Rscript
# Figure 1, panel C (right) -- cell-free transcript density of one Xenium section.
#
# Xenium transcripts (proseg's observed coordinates, qv >= 20) of section OR52YW26_1_4 counted in 10 um bins, smoothed with a normalized Gaussian
# (sd 20 um; the value grid and an occupancy mask are smoothed with the same kernel and divided, so tissue edges are not darkened, and bins whose
# smoothed occupancy is < 0.15 are dropped), viridis colour scale clipped at the 2nd and 98th percentiles of in-tissue bins. Same computation and
# drawing as xenium/transcript_density_smoothed.R (Rscript xenium/transcript_density_smoothed.R OR52YW26_1_4), which wrote
#   xenium/transcript_density/transcript_density_OR52YW26_1_4_10um_smooth20.png
# The panel extent and scale bar come from the proseg cell centroids (cells with >= 1 transcript), as in that script.
#
# Inputs: proseg output of the section (cell-metadata.parquet, expected-counts.csv.gz; transcript-metadata.parquet if the binned table must be made)
#         and the binned table tx_density_<section>_10um.csv.gz (common.R). If the table is missing it is made with panel_C_bin_transcripts.py,
#         which needs Python 3 with numpy and pyarrow (set FIG1_PYTHON to that interpreter).
# Output: FIG1_OUT/panel_C_transcript_density_OR52YW26_1_4_10um_smooth20.{pdf,png}
#
#   Rscript figures/figure1/panel_C_xenium_transcript_density.R
source(here::here("figures/figure1/common.R"))
suppressMessages(library(Matrix))
SECTION <- "OR52YW26_1_4"
BIN <- 10          # um
QV_MIN <- 20       # transcript quality floor used when binning
SIGMA_UM <- 20     # Gaussian sd, um
MASK_MIN <- 0.15   # minimum smoothed occupancy kept
sec_dir <- file.path(PROSEG_DIR, SECTION)

# ---- binned transcripts -----------------------------------------------------------------------------------------------------------------------------
dens_csv <- file.path(TX_DENSITY_DIR, paste0("tx_density_", SECTION, "_", BIN, "um.csv.gz"))
if (!file.exists(dens_csv)) {
  dens_csv <- file.path(FIG1_OUT, basename(dens_csv))
  tmp <- tempfile(fileext = ".csv")
  st <- system2(Sys.getenv("FIG1_PYTHON", "python3"), c(here::here("figures/figure1/panel_C_bin_transcripts.py"),
                                                         file.path(sec_dir, "transcript-metadata.parquet"), tmp, BIN, QV_MIN))
  stopifnot("transcript binning failed" = st == 0)
  stopifnot(system2("gzip", c("-c", tmp), stdout = dens_csv) == 0)
}
dens <- read_csv(dens_csv, show_col_types = FALSE) %>% orient_xy(SECTION)
message(sprintf("%s: %d occupied %d um bins, %s transcripts", SECTION, nrow(dens), BIN, format(sum(dens$n), big.mark = ",")))

# ---- cells (panel extent and scale bar) -------------------------------------------------------------------------------------------------------------
cm <- nanoparquet::read_parquet(file.path(sec_dir, "cell-metadata.parquet"))
m <- Matrix::readMM(gzfile(file.path(sec_dir, "expected-counts.csv.gz")))
stopifnot(nrow(m) == nrow(cm))
cells <- tibble(x = cm$centroid_x, y = cm$centroid_y, n_tx = Matrix::rowSums(m)) %>% orient_xy(SECTION) %>% filter(n_tx > 0)

# ---- normalized Gaussian smoothing on the bin grid --------------------------------------------------------------------------------------------------
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
  conv_cols <- function(M) {                       # convolve along the columns of M
    pad <- (length(k) - 1L) / 2L
    apply(M, 2, function(v) {
      z <- stats::filter(c(rep(0, pad), v, rep(0, pad)), k, sides = 2)
      as.numeric(z)[(pad + 1L):(pad + length(v))]
    })
  }
  Vs <- conv_cols(t(conv_cols(t(V))))              # separable: x then y
  Ws <- conv_cols(t(conv_cols(t(W))))
  val <- ifelse(Ws > MASK_MIN, Vs / pmax(Ws, 1e-9), NA_real_)
  # expand_grid() varies its last variable fastest, matching the column-major order of the [nx, ny] matrix
  expand_grid(iy = seq_len(ny) - 1L, ix = seq_len(nx) - 1L) %>%
    mutate(x = ix * bin + x0, y = iy * bin + y0, n = as.vector(val)) %>%
    filter(!is.na(n)) %>% select(x, y, n)
}
ds <- smooth_grid(dens, BIN, SIGMA_UM)
floor_n <- 100 * (BIN / 20)^2                      # "in tissue": >= 100 transcripts per 20 um bin, scaled to the bin area
lim <- unname(quantile(ds$n[ds$n >= floor_n], c(0.02, 0.98)))
brk <- signif(seq(lim[1], lim[2], length.out = 3), 2)

# ---- plot -------------------------------------------------------------------------------------------------------------------------------------------
xr <- range(cells$x); yr <- range(cells$y); bar_um <- 1000
bx <- xr[1] + 0.03 * diff(xr); by <- yr[1] + 0.04 * diff(yr)
PW <- 3.4; PH <- max(1.5, PW * diff(yr) / diff(xr) * 1.08)
p <- ggplot(ds, aes(x, y, fill = n)) +
  geom_tile(width = BIN, height = BIN) +            # geom_tile, not geom_raster: the masked grid has gaps
  scale_fill_gradientn(colours = viridisLite::viridis(256), limits = lim, oob = scales::squish, breaks = brk, labels = scales::comma(brk),
                       name = paste0("transcripts\nper ", BIN, " µm bin\n(", SIGMA_UM, " µm smooth)")) +
  coord_equal() +
  scale_x_continuous(expand = fig_expand()) + scale_y_continuous(expand = fig_expand()) +
  theme_fig() +
  theme(axis.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.55, "cm"),
        legend.title = element_text(size = FIG_PT_AXIS_TITLE), legend.text = element_text(size = FIG_PT_AXIS_TEXT),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.1, "cm")) +
  annotate("segment", x = bx, xend = bx + bar_um, y = by, yend = by, colour = FIG_INK_SECONDARY, linewidth = 0.5) +
  annotate("text", x = bx + bar_um / 2, y = by + 0.028 * diff(yr), label = "1 mm", size = fig_pt(5), colour = FIG_INK_SECONDARY)
out <- file.path(FIG1_OUT, paste0("panel_C_transcript_density_", SECTION, "_", BIN, "um_smooth", SIGMA_UM))
fig_save(p, out, width = PW, height = PH)
message("saved ", out)
