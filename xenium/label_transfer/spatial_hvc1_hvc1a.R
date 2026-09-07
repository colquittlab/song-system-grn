## Combined spatial map of Glut-HVC-1 and Glut-HVC-1a, no-Arco-4 RCTD reference.
##
## Same data, orientation, background and theme as
## proseg_spatial_maps_no_arco4.qmd; written as a standalone script so the
## figures can be added without re-rendering that whole notebook. Output goes
## into the notebook's own directory, under names it does not use.
##
## Two views of "combined": the two types overlaid in distinct colours (which
## one is where), and the union pooled into a single colour (the HVC field as
## one population). Each in the permissive `confident` set and in the
## nCount >= 400 depth-filtered set.

library(tidyverse)
library(cowplot)
library(here)

script_name = "proseg_spatial_maps_no_arco4"
hpc_dir = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_no_arco4")
out_dir = here::here("xenium/label_transfer", script_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

focus_types = c("Glut-HVC-1", "Glut-HVC-1a")

section_order = c("OR52YW26_1_4", "OR52YW26_1_7", "OR52YW26_2_2", "OR52YW26_2_4", "OR52YW26_2_7",
                  "OR69PU4_1_4", "OR69PU4_1_7", "OR69PU4_2_2", "OR69PU4_2_4", "OR69PU4_2_7")
confident_classes = c("singlet", "doublet_certain")
deep_min_counts = 400

## ---- load -----------------------------------------------------------------

rctd = read_csv(file.path(hpc_dir, "rctd_all.csv.gz"), show_col_types = FALSE) %>%
  dplyr::select(cell, spot_class, first_type)
md = read_csv(file.path(hpc_dir, "proseg_cell_metadata.csv.gz"), show_col_types = FALSE) %>%
  dplyr::select(cell, section_id, x_centroid, y_centroid)
ncount = read_csv(file.path(hpc_dir, "proseg_ncount.csv.gz"), show_col_types = FALSE)

d = md %>%
  left_join(rctd, by = "cell") %>%
  left_join(ncount, by = "cell") %>%
  mutate(section_id = factor(section_id, levels = section_order),
         confident = spot_class %in% confident_classes)

## ---- common orientation (verbatim from the notebook) ----------------------

section_transform = tribble(
  ~section_id,     ~rotate,  ~flip_h,
  "OR52YW26_1_4",  "none",   FALSE,
  "OR52YW26_1_7",  "ccw90",  FALSE,
  "OR52YW26_2_2",  "ccw90",  TRUE,
  "OR52YW26_2_4",  "ccw90",  FALSE,
  "OR52YW26_2_7",  "cw90",   TRUE,
  "OR69PU4_1_4",   "none",   FALSE,
  "OR69PU4_1_7",   "ccw90",  TRUE,
  "OR69PU4_2_2",   "cw90",   FALSE,
  "OR69PU4_2_4",   "ccw90",  FALSE,
  "OR69PU4_2_7",   "ccw90",  FALSE,
)
stopifnot(setequal(section_transform$section_id, section_order))

d = d %>%
  left_join(section_transform, by = "section_id") %>%
  mutate(
    x_rot = case_when(rotate == "none"  ~ x_centroid,
                      rotate == "ccw90" ~ -y_centroid,
                      rotate == "cw90"  ~ y_centroid),
    y_rot = case_when(rotate == "none"  ~ y_centroid,
                      rotate == "ccw90" ~ x_centroid,
                      rotate == "cw90"  ~ -x_centroid),
    x_plot = if_else(flip_h, -x_rot, x_rot),
    y_plot = y_rot) %>%
  dplyr::select(-rotate, -flip_h, -x_rot, -y_rot)

stopifnot(all(focus_types %in% unique(d$first_type)))

## ---- palette / theme ------------------------------------------------------
## Glut-HVC-1 takes the lab's standard HVC colour, `position_colors["hvc"]`
## (orange, "#ff7f00") from /opt/colquittlab/utils/R/common_aesthetics.R, so
## this figure reads consistently alongside anything using that convention.
## Glut-HVC-1a is black -- unused elsewhere in this figure set's palette and
## clearly separable from the orange against the grey background.

type_colors = setNames(c("#FF7F00", "#000000"), focus_types)
pooled_colour = "#000000"

theme_spatial = theme_cowplot() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5))

section_plot = function(sid, hl, point_colour = NULL) {
  bg = d %>% filter(section_id == sid)
  hl_sid = hl %>% filter(section_id == sid)

  p = ggplot() +
    geom_point(data = bg, aes(x_plot, y_plot),
               colour = "grey85", size = 0.03, alpha = 0.4) +
    coord_equal() + labs(title = sid) + theme_spatial

  if (nrow(hl_sid) == 0) return(p)

  if (is.null(point_colour)) {
    ## Explicit draw order rather than one geom_point: with two intermingled
    ## types a single layer hides whichever rows happen to come last in the
    ## data. Glut-HVC-1a goes on top; flip the two calls to see the reverse.
    p +
      geom_point(data = filter(hl_sid, first_type == "Glut-HVC-1"),
                 aes(x_plot, y_plot, colour = first_type), size = 0.25, alpha = 0.85) +
      geom_point(data = filter(hl_sid, first_type == "Glut-HVC-1a"),
                 aes(x_plot, y_plot, colour = first_type), size = 0.25, alpha = 0.85) +
      scale_colour_manual(values = type_colors, name = NULL, limits = focus_types,
                          guide = "none")
  } else {
    p + geom_point(data = hl_sid, aes(x_plot, y_plot),
                   colour = point_colour, size = 0.3, alpha = 0.85)
  }
}

legend_plot = ggplot(tibble(first_type = factor(focus_types, levels = focus_types)),
                     aes(1, 1, colour = first_type)) +
  geom_point(size = 3) +
  scale_colour_manual(values = type_colors, name = NULL)
legend = cowplot::get_legend(legend_plot)

## ---- figures --------------------------------------------------------------

highlight = d %>% filter(confident, first_type %in% focus_types) %>%
  mutate(first_type = factor(first_type, levels = focus_types))

message("confident cell counts:")
print(highlight %>% count(first_type))
print(highlight %>% filter(nCount >= deep_min_counts) %>% count(first_type))

make_two_colour = function(hl, fname, label) {
  panels = map(section_order, section_plot, hl = hl)
  grid = cowplot::plot_grid(plotlist = panels, ncol = 5)
  p = cowplot::plot_grid(grid, legend, ncol = 2, rel_widths = c(1, 0.12)) +
    cowplot::draw_figure_label(label, position = "top", size = 14)
  ggsave(file.path(out_dir, fname), p, width = 22, height = 9.5, dpi = 250, bg = "white")
}

make_pooled = function(hl, fname, label) {
  panels = map(section_order, section_plot, hl = hl, point_colour = pooled_colour)
  p = cowplot::plot_grid(plotlist = panels, ncol = 5) +
    cowplot::draw_figure_label(label, position = "top", size = 14)
  ggsave(file.path(out_dir, fname), p, width = 20, height = 9.5, dpi = 250, bg = "white")
}

n_lab = function(hl) {
  cnt = hl %>% count(first_type)
  paste0(cnt$first_type, " n = ", cnt$n, collapse = ";  ")
}

make_two_colour(
  highlight, "spatial_Glut_HVC_1_and_1a.png",
  paste0("Glut-HVC-1 (orange) + Glut-HVC-1a (black) — ", n_lab(highlight)))

make_pooled(
  highlight, "spatial_Glut_HVC_1_and_1a_pooled.png",
  paste0("Glut-HVC-1 + Glut-HVC-1a pooled  (n = ", nrow(highlight), " confident cells)"))

highlight_deep = highlight %>% filter(nCount >= deep_min_counts)

make_two_colour(
  highlight_deep, "spatial_Glut_HVC_1_and_1a_gte400.png",
  paste0("Glut-HVC-1 (orange) + Glut-HVC-1a (black), nCount >= ", deep_min_counts,
         " — ", n_lab(highlight_deep)))

make_pooled(
  highlight_deep, "spatial_Glut_HVC_1_and_1a_pooled_gte400.png",
  paste0("Glut-HVC-1 + Glut-HVC-1a pooled  (n = ", nrow(highlight_deep),
         " confident cells, nCount >= ", deep_min_counts, ")"))

message("wrote:")
print(list.files(out_dir, pattern = "Glut_HVC_1_and_1a"))
