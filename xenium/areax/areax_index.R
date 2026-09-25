## Locate Area X in the ten sagittal Xenium sections, for schematic use.
##
## Index  : mean z(up markers) - mean z(down markers), ZEBrA Area X markers
##          intersected with the 425-gene panel (see README.md for provenance;
##          public portal gene list only, no other project's data).
## Scope  : computed ONLY within striatal (LGE) RCTD calls. ZEBrA's contrast is
##          Area X vs the tissue adjacent to it, so the markers separate X from
##          surrounding striatum and not striatum from everything else.
## z-score: within the LGE cells of each section, so the index ranks cells
##          against the striatum they sit in rather than against the whole brain.

suppressMessages({ library(tidyverse); library(Matrix); library(cowplot); library(here) })
source(here::here("config/paths.R")); select = dplyr::select

hpc_dir  = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid")
gene_dir = here::here("xenium/label_transfer/hpc_rctd_proseg_hybrid/gene_lists")
out_dir  = here::here("xenium/areax", "areax_index")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GTE400 = 400
## Striatal population, set from the marker check below rather than from the
## label names. GABA-LGE-1 is NOT striatal despite the LGE prefix (PPP1R1B
## 0.42 vs GABA-LGE-2's 2.31, DRD1 0.47 vs 1.47, PDYN 0.04 vs 0.42) and was
## diluting the z-scores. Glut-GABA carries the strongest striatal signature of
## any type here (PPP1R1B 2.52, PENK 1.87, DRD2 2.10) and reads as the
## D2/indirect-pathway MSNs, with GABA-LGE-2 as D1/direct (DRD1 1.47,
## PDYN 0.42, DRD2 0.20).
LGE_TYPES = c("GABA-LGE-2", "Glut-GABA")
LGE_VARIANTS = list(d1_d2 = c("GABA-LGE-2","Glut-GABA"), d1_only = "GABA-LGE-2")
MATURITY = c("NEFL","SV2B","UCHL1")   # neuronal maturity, not regional
section_order = c("OR52YW26_1_4","OR52YW26_1_7","OR52YW26_2_2","OR52YW26_2_4","OR52YW26_2_7",
                  "OR69PU4_1_4","OR69PU4_1_7","OR69PU4_2_2","OR69PU4_2_4","OR69PU4_2_7")

panel = read_csv(XENIUM_PANEL_CSV, show_col_types=FALSE)$Gene
ax = read_csv(here::here("xenium/areax","zebra_areax_markers.csv"), show_col_types=FALSE)
up_all = ax %>% filter(regulation=="Up")   %>% pull(gene) %>% unique() %>% intersect(panel)
dn_all = ax %>% filter(regulation=="Down") %>% pull(gene) %>% unique() %>% intersect(panel)
up = setdiff(up_all, MATURITY); dn = setdiff(dn_all, MATURITY)
cat("up markers on panel:", length(up_all), "-> after dropping maturity genes:", length(up), "\n")
cat("  ", paste(up, collapse=", "), "\n")
cat("down markers on panel:", length(dn), ":", paste(dn, collapse=", "), "\n")

## striatal identity check genes (independent of the ZEBrA list)
STRIATAL = intersect(c("PPP1R1B","PENK","PDYN","DRD1","DRD2","RGS9","SCN4B","FOXP1","FOXP2"), panel)
need = unique(c(up_all, dn_all, STRIATAL))

## ---- metadata + calls -----------------------------------------------------
rctd = read_csv(file.path(hpc_dir,"rctd_all.csv.gz"), show_col_types=FALSE) %>%
  select(cell, spot_class, first_type)
md = read_csv(file.path(hpc_dir,"proseg_cell_metadata.csv.gz"), show_col_types=FALSE) %>%
  select(cell, section_id, x_centroid, y_centroid)
nc = read_csv(file.path(hpc_dir,"proseg_ncount.csv.gz"), show_col_types=FALSE)
meta = md %>% left_join(rctd,by="cell") %>% left_join(nc,by="cell") %>%
  mutate(confident = spot_class %in% c("singlet","doublet_certain"), gte400 = nCount >= GTE400)

## ---- pull only the genes we need, per section ------------------------------
load_expr = function(sid) {
  d = file.path(path.expand(XENIUM_PROSEG_DIR), sid)
  g = readLines(file.path(gene_dir, paste0(sid,".txt"))); g = g[nzchar(g)]
  m = Matrix::readMM(gzfile(file.path(d,"expected-counts.csv.gz")))
  m = as(t(m), "CsparseMatrix"); rownames(m) = g
  m = m[intersect(need, g), , drop=FALSE]
  cm = nanoparquet::read_parquet(file.path(d,"cell-metadata.parquet"))
  colnames(m) = paste0(sid,"_proseg_",cm$cell)
  m
}
## Cached: pulling 10 MatrixMarket files to keep ~25 genes is minutes per run.
expr_cache = file.path(out_dir, "marker_expr_cache.rds")
if (file.exists(expr_cache)) {
  expr = readRDS(expr_cache); message("using cached marker expression")
} else {
  message("loading marker expression for ", length(section_order), " sections...")
  expr = map(section_order, load_expr); names(expr) = section_order
  saveRDS(expr, expr_cache)
}

## ---- sanity check: are the LGE calls actually striatal? --------------------
strcheck = map_dfr(section_order, function(sid) {
  m = expr[[sid]]; mm = meta %>% filter(section_id==sid, confident, gte400)
  cells = intersect(colnames(m), mm$cell)
  cpm = t(t(m[, cells, drop=FALSE]) / mm$nCount[match(cells, mm$cell)]) * 1000
  tibble(section_id=sid, cell=cells, first_type=mm$first_type[match(cells, mm$cell)]) %>%
    bind_cols(as_tibble(as.matrix(t(log1p(cpm[STRIATAL, , drop=FALSE])))))
})
str_summary = strcheck %>% group_by(first_type) %>%
  summarize(n=n(), across(all_of(STRIATAL), mean), .groups="drop") %>%
  arrange(desc(PPP1R1B))
cat("\n=== mean log1p(CPM) of striatal markers by cell type (top 12) ===\n")
print(as.data.frame(str_summary %>% slice_head(n=12) %>%
  mutate(across(where(is.numeric), ~round(.,2)))), row.names=FALSE)
write_csv(str_summary, file.path(out_dir,"striatal_marker_check_by_type.csv"))

## ---- Area X index within LGE cells ----------------------------------------
score_section = function(sid, up_set, dn_set) {
  m = expr[[sid]]
  mm = meta %>% filter(section_id==sid, confident, gte400, first_type %in% LGE_TYPES)
  cells = intersect(colnames(m), mm$cell)
  if (length(cells) < 50) return(tibble())
  mm = mm[match(cells, mm$cell), ]
  cpm = t(t(m[, cells, drop=FALSE]) / mm$nCount) * 1000
  lg = log1p(cpm)
  z = function(v) { s = sd(v); if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - mean(v))/s }
  zu = t(apply(lg[intersect(up_set, rownames(lg)), , drop=FALSE], 1, z))
  zd = t(apply(lg[intersect(dn_set, rownames(lg)), , drop=FALSE], 1, z))
  mm %>% mutate(up_score = colMeans(zu),
                down_score = if (nrow(zd)) colMeans(zd) else 0,
                areax_index = up_score - down_score)
}

lge = map_dfr(section_order, score_section, up_set=up, dn_set=dn)
cat("\nLGE cells scored per section:\n")
print(as.data.frame(lge %>% count(section_id)), row.names=FALSE)

## sensitivity: with the maturity genes included
lge_all = map_dfr(section_order, score_section, up_set=up_all, dn_set=dn_all) %>%
  select(cell, areax_index_all = areax_index)
sens = lge %>% select(cell, section_id, areax_index) %>% left_join(lge_all, by="cell")
cat(sprintf("\nprimary vs maturity-inclusive index: Pearson r = %.3f (n=%d)\n",
            cor(sens$areax_index, sens$areax_index_all), nrow(sens)))
write_csv(sens, file.path(out_dir,"areax_index_gene_sensitivity.csv"))

## ---- orientation (same transform as the spatial maps) ----------------------
st = tribble(~section_id,~rotate,~flip_h,
  "OR52YW26_1_4","none",FALSE, "OR52YW26_1_7","ccw90",FALSE,
  "OR52YW26_2_2","ccw90",TRUE, "OR52YW26_2_4","ccw90",FALSE,
  "OR52YW26_2_7","cw90",TRUE,  "OR69PU4_1_4","none",FALSE,
  "OR69PU4_1_7","ccw90",TRUE,  "OR69PU4_2_2","cw90",FALSE,
  "OR69PU4_2_4","ccw90",FALSE, "OR69PU4_2_7","ccw90",FALSE)
orient = function(df) df %>% left_join(st, by="section_id") %>%
  mutate(x_rot=case_when(rotate=="none"~x_centroid,rotate=="ccw90"~-y_centroid,rotate=="cw90"~y_centroid),
         y_rot=case_when(rotate=="none"~y_centroid,rotate=="ccw90"~x_centroid,rotate=="cw90"~-x_centroid),
         x_plot=if_else(flip_h,-x_rot,x_rot), y_plot=y_rot) %>%
  select(-rotate,-flip_h,-x_rot,-y_rot)
lge = orient(lge); bg = orient(meta %>% mutate(section_id=factor(section_id,levels=section_order)))
lge$section_id = factor(lge$section_id, levels=section_order)
write_csv(lge %>% select(cell, section_id, first_type, x_plot, y_plot,
                         up_score, down_score, areax_index),
          file.path(out_dir,"areax_index_lge_cells.csv.gz"))

theme_spatial = theme_cowplot() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.line=element_blank(), axis.ticks=element_blank(),
        plot.title=element_text(size=10,hjust=0.5))

## continuous index
panels = map(section_order, function(sid) {
  b = bg %>% filter(section_id==sid); h = lge %>% filter(section_id==sid) %>% arrange(areax_index)
  ggplot() +
    geom_point(data=b, aes(x_plot,y_plot), colour="grey88", size=0.03, alpha=0.4) +
    geom_point(data=h, aes(x_plot,y_plot,colour=areax_index), size=0.35, alpha=0.9) +
    scale_colour_gradient2(low="#2166AC", mid="grey92", high="#B2182B", midpoint=0, name="Area X\nindex") +
    coord_equal() + labs(title=paste0(sid," (LGE n=",nrow(h),")")) + theme_spatial
})
leg = cowplot::get_legend(panels[[1]] + theme(legend.position="right"))
p = cowplot::plot_grid(cowplot::plot_grid(plotlist=map(panels, ~ .x + theme(legend.position="none")), ncol=5),
                       leg, ncol=2, rel_widths=c(1,0.08))
ggsave(file.path(out_dir,"areax_index_continuous.png"), p, width=24, height=10, dpi=300, bg="white")

## top-quantile highlight -- the schematic-facing view
QTOP = 0.80
thr = lge %>% group_by(section_id) %>% summarize(t = quantile(areax_index, QTOP), .groups="drop")
lge = lge %>% left_join(thr, by="section_id") %>% mutate(top = areax_index >= t)
panels2 = map(section_order, function(sid) {
  b = bg %>% filter(section_id==sid); h = lge %>% filter(section_id==sid)
  ggplot() +
    geom_point(data=b, aes(x_plot,y_plot), colour="grey88", size=0.03, alpha=0.4) +
    geom_point(data=h %>% filter(!top), aes(x_plot,y_plot), colour="grey65", size=0.3, alpha=0.6) +
    geom_point(data=h %>% filter(top), aes(x_plot,y_plot), colour="#B2182B", size=0.55, alpha=0.9) +
    coord_equal() + labs(title=paste0(sid," (top ",round((1-QTOP)*100),"% of LGE)")) + theme_spatial
})
ggsave(file.path(out_dir,"areax_index_top_quantile.png"),
       cowplot::plot_grid(plotlist=panels2, ncol=5), width=24, height=10, dpi=300, bg="white")

## How spatially clustered is the top quantile? An sd-based spread is far too
## blunt for this -- a tight focus inside a large field barely moves it. Median
## nearest-neighbour distance among the top cells, against size-matched random
## draws from the same striatal cells, actually responds to a focus.
nn_med = function(x, y) {
  if (length(x) < 5) return(NA_real_)
  m = as.matrix(dist(cbind(x, y))); diag(m) = Inf
  median(apply(m, 1, min))
}
set.seed(1)
conc = map_dfr(section_order, function(sid) {
  h = lge %>% filter(section_id == sid)
  tp = h %>% filter(top)
  if (nrow(tp) < 5) return(tibble(section_id = sid))
  obs = nn_med(tp$x_plot, tp$y_plot)
  nul = replicate(20, { i = sample(nrow(h), nrow(tp)); nn_med(h$x_plot[i], h$y_plot[i]) })
  tibble(section_id = sid, n_top = nrow(tp), nn_obs = obs, nn_null = mean(nul),
         nn_ratio = obs / mean(nul))          # < 1 means clustered
}) %>% arrange(nn_ratio)
cat("\n=== spatial clustering of the top-index striatal cells ===\n")
cat("(nn_ratio = median NN distance of top cells / size-matched random; <1 = clustered)\n")
print(as.data.frame(conc %>% mutate(across(where(is.numeric), ~round(.,2)))), row.names=FALSE)
write_csv(conc, file.path(out_dir,"areax_top_concentration.csv"))
cat("\nDONE\n")
