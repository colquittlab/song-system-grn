## Turn the Area X index focus into an outline usable for a schematic.
##
## For each section whose top-index striatal cells are spatially clustered
## (nn_ratio < NN_CUT from areax_index.R), fit a 2D density to those cells,
## take the contour at CONTOUR_FRAC of peak density, and keep the single
## polygon containing the peak. That polygon is the Area X outline.
##
## Outputs both oriented plotting coordinates (x_plot/y_plot, the common
## posterior-left/dorsal-up frame used by the spatial maps) and RAW Xenium
## coordinates (x_centroid/y_centroid), since overlaying on the original images
## needs the untransformed frame.

suppressMessages({ library(tidyverse); library(MASS); library(cowplot); library(here) })
source(here::here("config/paths.R")); select = dplyr::select

hpc_dir = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid")
dir_in  = here::here("xenium/areax", "areax_index")
out_dir = dir_in
NN_CUT = 0.80; CONTOUR_FRAC = 0.5; QTOP = 0.80

lge  = read_csv(file.path(dir_in,"areax_index_lge_cells.csv.gz"), show_col_types=FALSE)
conc = read_csv(file.path(dir_in,"areax_top_concentration.csv"), show_col_types=FALSE)
raw  = read_csv(file.path(hpc_dir,"proseg_cell_metadata.csv.gz"), show_col_types=FALSE) %>%
  select(cell, x_centroid, y_centroid)
lge = lge %>% left_join(raw, by="cell")

keep = conc %>% filter(!is.na(nn_ratio), nn_ratio < NN_CUT) %>% pull(section_id)
cat("sections with a clustered Area X focus (nn_ratio <", NN_CUT, "):\n  ",
    paste(keep, collapse=", "), "\n")
cat("sections without one (no Area X at that medio-lateral level):\n  ",
    paste(setdiff(conc$section_id, keep), collapse=", "), "\n\n")

outline_one = function(sid) {
  h = lge %>% filter(section_id==sid)
  tp = h %>% filter(areax_index >= quantile(areax_index, QTOP))
  k = kde2d(tp$x_plot, tp$y_plot, n=200)
  pk = which(k$z == max(k$z), arr.ind=TRUE)[1,]
  px = k$x[pk[1]]; py = k$y[pk[2]]
  cl = contourLines(k$x, k$y, k$z, levels = CONTOUR_FRAC*max(k$z))
  ## keep the polygon that encloses the density peak
  inside = map_lgl(cl, ~ sp::point.in.polygon(px, py, .x$x, .x$y) > 0)
  if (!any(inside)) return(NULL)
  poly = cl[inside][[which.max(map_dbl(cl[inside], ~ length(.x$x)))]]
  pol = tibble(section_id=sid, x_plot=poly$x, y_plot=poly$y)
  ## cells inside the outline, and their raw coordinates
  inx = sp::point.in.polygon(h$x_plot, h$y_plot, poly$x, poly$y) > 0
  list(poly = pol,
       cells = h[inx, ],
       summ = tibble(section_id=sid, n_cells_in_X=sum(inx),
                     peak_x_plot=px, peak_y_plot=py,
                     centroid_x_raw=mean(h$x_centroid[inx]),
                     centroid_y_raw=mean(h$y_centroid[inx]),
                     mean_index_in=mean(h$areax_index[inx]),
                     mean_index_out=mean(h$areax_index[!inx]),
                     width_um=diff(range(poly$x)), height_um=diff(range(poly$y))))
}

res = map(keep, outline_one) %>% compact()
polys = map_dfr(res, "poly"); summ = map_dfr(res, "summ")
cells_in = map_dfr(res, "cells")

## outline in raw coordinates too: invert the orientation transform
st = tribble(~section_id,~rotate,~flip_h,
  "OR52YW26_1_4","none",FALSE, "OR52YW26_1_7","ccw90",FALSE,
  "OR52YW26_2_2","ccw90",TRUE, "OR52YW26_2_4","ccw90",FALSE,
  "OR52YW26_2_7","cw90",TRUE,  "OR69PU4_1_4","none",FALSE,
  "OR69PU4_1_7","ccw90",TRUE,  "OR69PU4_2_2","cw90",FALSE,
  "OR69PU4_2_4","ccw90",FALSE, "OR69PU4_2_7","ccw90",FALSE)
polys_raw = polys %>% left_join(st, by="section_id") %>%
  mutate(x_un = if_else(flip_h, -x_plot, x_plot), y_un = y_plot,
         x_centroid = case_when(rotate=="none" ~ x_un, rotate=="ccw90" ~  y_un, rotate=="cw90" ~ -y_un),
         y_centroid = case_when(rotate=="none" ~ y_un, rotate=="ccw90" ~ -x_un, rotate=="cw90" ~  x_un)) %>%
  select(section_id, x_centroid, y_centroid)

write_csv(polys,     file.path(out_dir,"areax_outline_plotcoords.csv"))
write_csv(polys_raw, file.path(out_dir,"areax_outline_rawcoords.csv"))
write_csv(summ,      file.path(out_dir,"areax_outline_summary.csv"))
write_csv(cells_in %>% select(cell, section_id, first_type, areax_index,
                              x_centroid, y_centroid, x_plot, y_plot),
          file.path(out_dir,"areax_cells_in_outline.csv.gz"))

cat("=== Area X outline per section ===\n")
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~round(.,1)))), row.names=FALSE)

## verify the outline against the transform by round-tripping one section
chk = polys_raw %>% filter(section_id==keep[1]) %>% slice(1)
cat("\nraw-coordinate round-trip check for", keep[1], ": x=", round(chk$x_centroid,1),
    " y=", round(chk$y_centroid,1), "\n")

theme_spatial = theme_cowplot() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.line=element_blank(), axis.ticks=element_blank(),
        plot.title=element_text(size=11,hjust=0.5))
panels = map(keep, function(sid) {
  h = lge %>% filter(section_id==sid); pol = polys %>% filter(section_id==sid)
  s = summ %>% filter(section_id==sid)
  ggplot() +
    geom_point(data=h, aes(x_plot,y_plot,colour=areax_index), size=0.5, alpha=0.9) +
    scale_colour_gradient2(low="#2166AC", mid="grey92", high="#B2182B", midpoint=0, guide="none") +
    geom_path(data=pol, aes(x_plot,y_plot), colour="black", linewidth=0.7) +
    coord_equal() +
    labs(title=paste0(sid, "  (Area X: ", s$n_cells_in_X, " cells, ",
                      round(s$width_um), "x", round(s$height_um), " um)")) +
    theme_spatial
})
ggsave(file.path(out_dir,"areax_outline.png"),
       cowplot::plot_grid(plotlist=panels, ncol=2), width=14, height=11, dpi=300, bg="white")
cat("\nDONE\n")
