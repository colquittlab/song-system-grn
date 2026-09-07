# song-system-grn

## Figure conventions

**The enforceable parts live in `config/figure_theme.R`, not here.** Source it and
use it rather than re-declaring type, colour or theme in a plotting script:

```r
source(here::here("config/figure_theme.R"))

p <- ggplot(...) + theme_fig()
fig_save(p, file.path(out_dir, "my_figure"), width = 3, height = 2.04)
```

`fig_save()` writes the PDF and a 600-dpi PNG together so the two cannot drift.

The standard it encodes:

- **Arial throughout**, genuinely embedded. Real Arial is installed here
  (msttcorefonts). If it were not, fontconfig would silently substitute
  Liberation Sans -- metrically identical, so nothing would *look* wrong -- which
  is why `fig_check_font()` warns instead. Verify any figure with
  `pdffonts <file>.pdf`; expect `ArialMT` / `Arial-BoldMT`, `emb=yes`.
- **Axis tick labels 5 pt, axis titles 6 pt.**
- **Whitespace is kept tight.** Small plot margins, ~2% panel expansion
  (`fig_expand()`), and no floating annotations that force the panel open past
  the data. Figures are small: a single compact panel is ~3 x 2 in.
- **Colour by job.** Categorical hues from `FIG_PAL` in fixed slot order, never
  cycled. Charts where every series can touch every other (scatter, bubble,
  small multiples) are capped at `FIG_PAL_ALLPAIRS_MAX` = 3 slots, which is the
  number that clears the colour-vision-deficiency separation floors on the
  all-pairs list. Sequential = one hue light-to-dark; diverging =
  `FIG_DIV_LOW`/`FIG_DIV_MID`/`FIG_DIV_HIGH` with a **neutral gray** midpoint.
  Never a rainbow, never a hue at a diverging midpoint.
- **Identity is never carried by colour alone** -- direct labels or a legend as
  well.

### The unit trap

`element_text(size = )` is in **points**; `geom_text()`/`annotate(size = )` is in
**millimetres**. Writing `size = 6` in a geom gives ~17 pt. Use `fig_pt(6)`.

### When shrinking a figure

Text does not scale with the canvas. Cutting width without cutting the title,
subtitle and caption to the new measure just clips them off the right edge, and
`ggsave()` reports nothing. Re-budget the text, then re-render and *look at the
output*.

## Working conventions

- **Notebooks read their inputs read-only and write their own copy.** Nothing
  writes back to an upstream object.
- **Outputs are regenerated, not tracked.** `.gitignore` drops figures and large
  objects; small `.csv` summaries stay tracked so numbers are checkable without
  re-running anything.
- **Shared helpers are sourced from `config/`** (`paths.R`, `figure_theme.R`) or
  `scripts/`, via `source(here::here(...))`.
- Cluster labels are defined in one place,
  `snrna/naming/hybrid_division_naming.qmd`. Renaming one means propagating to
  the downstream `ct_order` lists and result tables -- and to the saved objects,
  which are not in git.
- **Cells are not independent samples.** Each dissection is a single library
  containing a few birds (souporcell `assignment` indexes birds *within* a
  library and is not comparable across libraries). Per-cell n inflates any test;
  report an effect size alongside a p-value.
