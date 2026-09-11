## Per-CELL-TYPE spatial localisation index (companion to spatial_variation.R,
## which works at the class level).
##
## For each type t in each section:
##   p = t's distribution over 400um spatial bins
##   q = ALL profiled cells' distribution over the same bins  (tissue shape and
##       sampling density, so an unrestricted type scores ~0 rather than
##       inheriting the outline of the section)
##   obs   = JSD(p || q)                      -- bounded [0,1] bits, no
##                                               normalisation needed
##   null  = JSD(p_rand || q), where p_rand is the same NUMBER of cells drawn
##           at random from all cells in that section
##   delta = obs - null                        -- the genuinely spatial part
##
## Every type is subsampled to the SAME n (N_FIX) before binning, and the null
## uses that same n, so a rare type is not scored as "localised" purely because
## few cells spread thinly across bins. Types with fewer than N_FIX cells in a
## section are skipped for that section.

suppressMessages({ library(tidyverse); library(here) })
source(here::here("config/paths.R"))
select = dplyr::select; set.seed(1)

hpc_dir = path.expand("~/hdd/rstudio/xenium/260811_brainard_adult-425g/hpc_rctd_proseg_hybrid")
out_dir = here::here("xenium/label_transfer", "spatial_variation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BIN = 400; N_FIX = 100; B = 30; GTE400 = 400
## N_FIX is both the subsample size and the per-section inclusion threshold.
## At 200 several types (e.g. Glut-CACNA1H-RA) cleared it in only one section,
## so their index rested on a single panel; 100 buys much wider section
## coverage. The matched-n null absorbs the extra sampling noise.

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
  t %in% c("Glut-NSC","Glut-NB","Glut-Im","GABA-Im") ~ "Developmental",
  TRUE ~ "Other")

d = md %>% left_join(rctd, by="cell") %>% left_join(nc, by="cell") %>%
  filter(spot_class %in% c("singlet","doublet_certain"), nCount >= GTE400) %>%
  mutate(class = classify(first_type),
         bin = paste(floor(x_centroid/BIN), floor(y_centroid/BIN)))

jsd = function(p, q) {
  m = (p+q)/2
  kl = function(a,b){ i = a>0; sum(a[i]*log2(a[i]/b[i])) }
  0.5*kl(p,m) + 0.5*kl(q,m)
}

res = map_dfr(unique(d$section_id), function(sid) {
  sec = d %>% filter(section_id == sid)
  bins = sort(unique(sec$bin))
  qv = table(factor(sec$bin, levels=bins)); q = as.numeric(qv)/sum(qv)
  pool = sec$bin
  tys = sec %>% count(first_type) %>% filter(n >= N_FIX) %>% pull(first_type)
  map_dfr(tys, function(ty) {
    tb = sec$bin[sec$first_type == ty]
    o = replicate(B, { s = sample(tb, N_FIX)
      jsd(as.numeric(table(factor(s, levels=bins)))/N_FIX, q) })
    nl = replicate(B, { s = sample(pool, N_FIX)
      jsd(as.numeric(table(factor(s, levels=bins)))/N_FIX, q) })
    tibble(section_id=sid, first_type=ty, class=classify(ty),
           n_cells=length(tb), obs=mean(o), null=mean(nl), delta=mean(o)-mean(nl))
  })
})

write_csv(res, file.path(out_dir,"spatial_localisation_by_type_section.csv"))

## NB: compute the sd BEFORE a column called `delta` is redefined -- summarize()
## evaluates in order, so `delta_sd = sd(delta)` after `delta = mean(delta)`
## silently takes the sd of a single scalar and returns NA.
per_type = res %>% group_by(first_type, class) %>%
  summarize(n_sections=n(), delta_sd=sd(delta), delta_mean=mean(delta),
            median_n=median(n_cells), .groups="drop") %>%
  rename(delta = delta_mean) %>% arrange(desc(delta))
write_csv(per_type, file.path(out_dir,"spatial_localisation_by_type.csv"))

cat("\n=== per-type spatial localisation index (delta), ranked ===\n")
print(as.data.frame(per_type %>% mutate(across(where(is.numeric), ~round(.,4)))), row.names=FALSE)

cat("\n=== distribution by class ===\n")
bycl = per_type %>% group_by(class) %>%
  summarize(n_types=n(), mean_delta=mean(delta), median_delta=median(delta),
            mean_sections=mean(n_sections),
            min_delta=min(delta), max_delta=max(delta), .groups="drop") %>%
  arrange(desc(mean_delta))
print(as.data.frame(bycl %>% mutate(across(where(is.numeric), ~round(.,4)))), row.names=FALSE)
write_csv(bycl, file.path(out_dir,"spatial_localisation_by_class_from_types.csv"))

cat("\n=== Wilcoxon on per-type deltas ===\n")
for (p in list(c("Glut","Astro"), c("Astro","GABA"), c("Glut","GABA"))) {
  a = per_type$delta[per_type$class==p[1]]; b = per_type$delta[per_type$class==p[2]]
  if (length(a)>1 && length(b)>1) {
    w = suppressWarnings(wilcox.test(a,b))
    cat(sprintf("%-6s (n=%2d, med %.3f) vs %-6s (n=%2d, med %.3f): p=%.3g\n",
        p[1], length(a), median(a), p[2], length(b), median(b), w$p.value))
  }
}
cat("\nDONE\n")
