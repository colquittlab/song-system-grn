## ---------------------------------------------------------------------------
## Shared figure conventions for this project.
##
##   source(here::here("config/figure_theme.R"))
##
## Same convention as config/paths.R: sourced read-only, defines values and
## functions, does nothing on its own.
##
## The point of this file is that the figure standards are CODE, not prose in a
## style guide that drifts. Anything enforceable lives here; anything that still
## needs judgement is in CLAUDE.md.
##
## The standard:
##   - Arial throughout, genuinely embedded (not a metric substitute).
##   - Axis tick labels 5 pt; axis titles 6 pt.
##   - Whitespace tight: small plot margins, ~2% panel expansion, no floating
##     annotations that force the panel open past the data.
##   - Colours from a CVD-validated categorical palette, assigned in fixed slot
##     order and never cycled.
## ---------------------------------------------------------------------------

suppressMessages({
  library(ggplot2)
  library(grid)
})

## --- Type -------------------------------------------------------------------

FIG_FONT <- "Arial"
FIG_PT_AXIS_TEXT <- 5   # tick labels
FIG_PT_AXIS_TITLE <- 6  # axis titles

## Fail loudly rather than let fontconfig silently substitute Liberation Sans.
## It is metrically identical to Arial, so nothing in the output would look
## wrong -- but it is not Arial, and a journal that asks for Arial means Arial.
## Verify with:  pdffonts <file>.pdf   (expect ArialMT / Arial-BoldMT, emb=yes)
fig_check_font <- function(font = FIG_FONT) {
  ok <- tryCatch(any(systemfonts::system_fonts()$family == font),
                 error = function(e) NA)
  if (isFALSE(ok)) {
    warning("Font '", font, "' not installed; fontconfig will substitute and ",
            "the output will not actually be ", font, ".", call. = FALSE)
  }
  invisible(ok)
}

## THE UNIT TRAP, written down because it has bitten this project already:
## element_text(size = ) is in POINTS, geom_text()/annotate(size = ) is in
## MILLIMETRES. Writing size = 6 in a geom gives ~17 pt, not 6 pt. Use fig_pt().
fig_pt <- function(pt) pt / .pt

## --- Colour -----------------------------------------------------------------
## Categorical slots in fixed order, never cycled. The first three clear the
## all-pairs colour-vision-deficiency and normal-vision separation floors, so
## any chart where every series can touch every other (scatter, bubble, small
## multiples) is capped at three; past that, fold to "Other" or facet.
FIG_PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")
FIG_PAL_ALLPAIRS_MAX <- 3

FIG_INK_PRIMARY <- "#0b0b0b"
FIG_INK_SECONDARY <- "#52514e"
FIG_INK_MUTED <- "#898781"
FIG_GRID <- "#e1e0d9"
FIG_AXIS <- "#c3c2b7"
FIG_SURFACE <- "#fcfcfb"
FIG_NEUTRAL_FILL <- "#c9c8c2"

## Diverging: two hues either side of a NEUTRAL GRAY midpoint. Never a hue at
## the midpoint, never a rainbow.
FIG_DIV_LOW <- "#2a78d6"
FIG_DIV_MID <- "#f0efec"
FIG_DIV_HIGH <- "#e34948"

## --- Theme ------------------------------------------------------------------

theme_fig <- function(base_size = 7, base_family = FIG_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = FIG_SURFACE, colour = NA),
      panel.background = element_rect(fill = FIG_SURFACE, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = FIG_GRID, linewidth = 0.25),
      axis.line = element_line(colour = FIG_AXIS, linewidth = 0.3),
      axis.ticks = element_line(colour = FIG_AXIS, linewidth = 0.3),
      axis.text = element_text(colour = FIG_INK_MUTED, size = FIG_PT_AXIS_TEXT),
      axis.title = element_text(colour = FIG_INK_SECONDARY, size = FIG_PT_AXIS_TITLE),
      axis.title.x = element_text(margin = margin(t = 1.5)),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      plot.title = element_text(colour = FIG_INK_PRIMARY, face = "bold",
                                size = base_size + 1, margin = margin(b = 1.5)),
      plot.subtitle = element_text(colour = FIG_INK_SECONDARY,
                                   size = base_size - 0.5, margin = margin(b = 2.5)),
      plot.caption = element_text(colour = FIG_INK_MUTED, size = base_size - 2.4,
                                  hjust = 0, margin = margin(t = 1.5)),
      plot.tag = element_text(colour = FIG_INK_PRIMARY, face = "bold",
                              size = base_size + 3),
      legend.text = element_text(colour = FIG_INK_SECONDARY),
      legend.title = element_text(colour = FIG_INK_SECONDARY),
      plot.margin = margin(2, 4, 2, 2)
    )
}

## Tight default expansion, so panels do not carry a ring of dead space.
fig_expand <- function(mult = 0.02) ggplot2::expansion(mult = mult)

## --- Saving -----------------------------------------------------------------
## Writes the PDF (vector, for the manuscript) and a high-dpi PNG (for looking
## at) from one call, so the two cannot drift apart. `family` is passed to
## cairo_pdf explicitly -- the device default is otherwise "sans", which would
## put a non-Arial face on anything the theme does not cover.
fig_save <- function(plot, path_noext, width, height, dpi = 600) {
  fig_check_font()
  ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf, family = FIG_FONT)
  ggsave(paste0(path_noext, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = FIG_SURFACE)
  invisible(paste0(path_noext, c(".pdf", ".png")))
}

## Geoms do not inherit the theme's family; set their default on source.
update_geom_defaults("text", list(family = FIG_FONT))
update_geom_defaults("label", list(family = FIG_FONT))
