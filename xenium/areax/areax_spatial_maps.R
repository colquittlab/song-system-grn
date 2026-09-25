## Area X plotted in the same style as xenium/label_transfer/proseg_spatial_maps*.qmd:
## all ten sections, common orientation, every profiled cell in grey behind the
## highlight, one fixed highlight colour, 5-column grid.
## Conventions copied from that notebook: grey85/size 0.03/alpha 0.4 background,
## size 0.3/alpha 0.85 highlight, single_type_colour = black, 20 x 9.5 at 250 dpi.

suppressMessages({ library(tidyverse); library(cowplot); library(here) })
source(here::here("config/paths.R")); select = dplyr::select

hpc_dir = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid")
dir_in  = here::here("xenium/areax", "areax_index")
out_dir = dir_in
section_order = c("OR52YW26_1_4","OR52YW26_1_7","OR52YW26_2_2","OR52YW26_2_4","OR52YW26_2_7",
                  "OR69PU4_1_4","OR69PU4_1_7","OR69PU4_2_2","OR69PU4_2_4","OR69PU4_2_7")
single_type_colour = "#000000"

md = read_csv(file.path(hpc_dir,"proseg_cell_metadata.csv.gz"), show_col_types=FALSE) %>%
  select(cell, section_id, x_centroid, y_centroid)
st = tribble(~section_id,~rotate,~flip_h,
  "OR52YW26_1_4","none",FALSE, "OR52YW26_1_7","ccw90",FALSE,
  "OR52YW26_2_2","ccw90",TRUE, "OR52YW26_2_4","ccw90",FALSE,
  "OR52YW26_2_7","cw90",TRUE,  "OR69PU4_1_4","none",FALSE,
  "OR69PU4_1_7","ccw90",TRUE,  "OR69PU4_2_2","cw90",FALSE,
  "OR69PU4_2_4","ccw90",FALSE, "OR69PU4_2_7","ccw90",FALSE)
d = md %>% left_join(st, by="section_id") %>%
  mutate(section_id=factor(section_id, levels=section_order),
         x_rot=case_when(rotate=="none"~x_centroid,rotate=="ccw90"~-y_centroid,rotate=="cw90"~y_centroid),
         y_rot=case_when(rotate=="none"~y_centroid,rotate=="ccw90"~x_centroid,rotate=="cw90"~-x_centroid),
         x_plot=if_else(flip_h,-x_rot,x_rot), y_plot=y_rot) %>%
  select(-rotate,-flip_h,-x_rot,-y_rot)

lge    = read_csv(file.path(dir_in,"areax_index_lge_cells.csv.gz"), show_col_types=FALSE)
inX    = read_csv(file.path(dir_in,"areax_cells_in_outline.csv.gz"), show_col_types=FALSE)
poly   = read_csv(file.path(dir_in,"areax_outline_plotcoords.csv"), show_col_types=FALSE)

theme_spatial = theme_cowplot() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.line=element_blank(), axis.ticks=element_blank(),
        plot.title=element_text(size=10,hjust=0.5))

section_plot = function(sid, hl, point_colour=single_type_colour, show_outline=FALSE) {
  bg = d %>% filter(section_id==sid)
  hl_sid = hl %>% filter(section_id==sid)
  p = ggplot() +
    geom_point(data=bg, aes(x_plot,y_plot), colour="grey85", size=0.03, alpha=0.4) +
    coord_equal() + labs(title=sid) + theme_spatial
  if (nrow(hl_sid))
    p = p + geom_point(data=hl_sid, aes(x_plot,y_plot), colour=point_colour, size=0.3, alpha=0.85)
  if (show_outline) {
    pl = poly %>% filter(section_id==sid)
    if (nrow(pl)) p = p + geom_path(data=pl, aes(x_plot,y_plot), colour="#B2182B", linewidth=0.5)
  }
  p
}

save_set = function(hl, fname, label, show_outline=FALSE) {
  panels = map(section_order, section_plot, hl=hl, show_outline=show_outline)
  p = cowplot::plot_grid(plotlist=panels, ncol=5) +
    cowplot::draw_figure_label(label, position="top", size=14)
  ggsave(file.path(out_dir, fname), p, width=20, height=9.5, dpi=250, bg="white")
  cat("wrote:", fname, "\n")
}

## 1. the Area X cells themselves
save_set(inX, "spatial_AreaX.png",
         paste0("Area X  (n = ", nrow(inX), " cells inside the delineated outline, 4 of 10 sections)"))
## 2. same, with the outline drawn
save_set(inX, "spatial_AreaX_outline.png",
         paste0("Area X with delineated outline  (n = ", nrow(inX), " cells)"),
         show_outline=TRUE)
## 3. the striatal field the index was computed within, for context
save_set(lge, "spatial_striatum_LGE.png",
         paste0("Striatal calls used for the index: GABA-LGE-2 + Glut-GABA  (n = ", nrow(lge), ")"))
cat("DONE\n")
