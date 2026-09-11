## Does subtype composition vary in space more for Glut than Astro than GABA?
##
## Statistic: RAREFIED, PERMUTATION-STANDARDISED SPATIAL COMPOSITION DIVERGENCE.
##
##   - bin each section into BIN-um squares
##   - within a class, draw K subtypes and reweight each to EQUAL total mass
##   - per bin, composition p_b over those K types; global composition q
##   - obs = mean_b JSD(p_b || q) / log2(K)        (normalised, 0..1)
##   - null = same quantity after permuting subtype labels among the SAME cells
##     (positions and bin sizes fixed) -> the divergence expected from finite
##     bin counts alone
##   - report delta = obs - null
##
## Why rarefy/reweight: Glut has 19 subtypes, GABA 13, Astro 3, and abundances
## differ ~100x. Both inflate raw divergence for the richer class, so K is
## equalised by drawing, and abundance by reweighting each drawn type to equal
## total mass (w_t = 1/n_t). Reweighting rather than subsampling matters: a
## first version downsampled to <=400 cells/type, which spread ~1200 cells over
## ~375 bins so almost no bin cleared MIN_BIN and every statistic came back NaN.
## Weighting equalises abundance exactly while keeping all cells for density.
##
## Why a null: with finite cells per bin, even a spatially random class shows
## divergence > 0. delta is the part that is actually spatial.

suppressMessages({ library(tidyverse); library(here) })
source(here::here("config/paths.R"))
select = dplyr::select; set.seed(1)

hpc_dir = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid")
out_dir = here::here("xenium/label_transfer", "spatial_variation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BIN       = 400    # um square bins (200 left GABA/Astro with too few
                   # well-sampled bins: sparser classes were dropping sections,
                   # which would itself bias the between-class comparison)
MIN_BIN   = 20     # min cells of the drawn set per bin
K         = 3      # subtypes per draw (Astro only has 3 -> this is the cap)
B         = 50     # draws
GTE400    = 400

rctd = read_csv(file.path(hpc_dir,"rctd_all.csv.gz"), show_col_types=FALSE) %>%
  select(cell, spot_class, first_type)
md = read_csv(file.path(hpc_dir,"proseg_cell_metadata.csv.gz"), show_col_types=FALSE) %>%
  select(cell, section_id, x_centroid, y_centroid)
nc = read_csv(file.path(hpc_dir,"proseg_ncount.csv.gz"), show_col_types=FALSE)

classify = function(t) case_when(
  str_detect(t, "^Glut-(DACH2|CACNA1H|SATB2)") ~ "Glut",
  str_detect(t, "^GABA-(LGE|MGE|CGE)")         ~ "GABA",
  str_detect(t, "^Astro")                      ~ "Astro",
  str_detect(t, "^Oligo|^OPC")                 ~ "Oligo/OPC",
  TRUE ~ NA_character_)

d = md %>% left_join(rctd, by="cell") %>% left_join(nc, by="cell") %>%
  filter(spot_class %in% c("singlet","doublet_certain"), nCount >= GTE400) %>%
  mutate(class = classify(first_type)) %>% filter(!is.na(class)) %>%
  mutate(bx = floor(x_centroid/BIN), by = floor(y_centroid/BIN))

jsd_stat = function(df) {
  ## mean JSD(bin composition || global), abundance-equalised by weighting each
  ## type to equal total mass, normalised by log2(K). Bin retention uses RAW
  ## counts so a bin still has to be genuinely well sampled.
  tab = df %>% count(bin, first_type) %>%
    pivot_wider(names_from=first_type, values_from=n, values_fill=0)
  m = as.matrix(tab[,-1, drop=FALSE])
  keep = rowSums(m) >= MIN_BIN
  if (sum(keep) < 5) return(NA_real_)
  m = m[keep,,drop=FALSE]
  tot = colSums(m)
  if (any(tot == 0)) return(NA_real_)
  W = sweep(m, 2, 1/tot, "*")          # each type contributes equal total mass
  P = W / rowSums(W)
  q = colSums(W) / sum(W)
  kl = function(a,b){ i = a>0; sum(a[i]*log2(a[i]/b[i])) }
  js = apply(P, 1, function(p){ mm = (p+q)/2; 0.5*kl(p,mm) + 0.5*kl(q,mm) })
  mean(js) / log2(ncol(m))
}

run = function(dat, label) {
  map_dfr(unique(dat$section_id), function(sid) {
    map_dfr(c("Glut","GABA","Astro","Oligo/OPC"), function(cl) {
      sub = dat %>% filter(section_id==sid, class==cl) %>% mutate(bin = paste(bx,by))
      tys = sub %>% count(first_type) %>% filter(n >= 100) %>% pull(first_type)
      if (length(tys) < K) return(tibble())
      res = map_dfr(seq_len(B), function(b) {
        pick = sample(tys, K)
        s = sub %>% filter(first_type %in% pick)
        obs = jsd_stat(s)
        s$first_type = sample(s$first_type)      # positions fixed, labels shuffled
        nul = jsd_stat(s)
        tibble(obs=obs, null=nul)
      })
      tibble(filter_set=label, section_id=sid, class=cl,
             n_subtypes_avail=length(tys),
             draws_usable=sum(is.finite(res$obs)),
             obs=mean(res$obs,na.rm=TRUE), null=mean(res$null,na.rm=TRUE),
             delta=mean(res$obs-res$null,na.rm=TRUE))
    })
  })
}

all_res = bind_rows(
  run(d, "confident_gte400"),
  run(d %>% filter(spot_class=="singlet"), "singlet_only_gte400"))

write_csv(all_res, file.path(out_dir,"spatial_variation_by_class_section.csv"))

summ = all_res %>% group_by(filter_set, class) %>%
  summarize(n_sections=sum(is.finite(delta)),
            delta_mean=mean(delta,na.rm=TRUE), delta_sd=sd(delta,na.rm=TRUE),
            obs_mean=mean(obs,na.rm=TRUE), null_mean=mean(null,na.rm=TRUE),
            .groups="drop") %>%
  arrange(filter_set, desc(delta_mean))
write_csv(summ, file.path(out_dir,"spatial_variation_summary.csv"))
cat("\n=== mean delta (spatial composition divergence above null) ===\n")
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~round(.,4)))), row.names=FALSE)

## paired tests across sections, primary filter only
prim = all_res %>% filter(filter_set=="confident_gte400") %>%
  select(section_id, class, delta) %>% pivot_wider(names_from=class, values_from=delta)
cat("\n=== paired across the 10 sections (primary filter) ===\n")
for (p in list(c("Glut","Astro"), c("Astro","GABA"), c("Glut","GABA"))) {
  a = prim[[p[1]]]; b = prim[[p[2]]]
  ok = is.finite(a) & is.finite(b); a = a[ok]; b = b[ok]
  if (length(a) < 3) { cat(sprintf("%-6s vs %-6s: too few usable sections\n", p[1], p[2])); next }
  tt = t.test(a, b, paired=TRUE)
  cat(sprintf("%-6s vs %-6s: mean diff %+0.4f   t=%6.2f  p=%.2g  (%d/%d sections %s>%s)\n",
      p[1], p[2], mean(a-b), tt$statistic, tt$p.value, sum(a>b), length(a), p[1], p[2]))
}
cat("\nDONE\n")
