## Build data/orthologs/loc_human_aliases.csv: a human gene symbol for each informative LOC model in
## the lonStrDom2 annotation, for use as a display label.
##
##   Rscript scripts/build_loc_aliases.R
##
## Why this is needed. NCBI gives a Gnomon model a symbol only when it can assign orthology
## confidently; otherwise it is left as LOC<GeneID> with just a protein product name. Many of those
## are informative -- LOC110480370 is annotated "multiple epidermal growth factor-like domains protein
## 6" -- but a LOC name says nothing on a figure, so they were excluded from every label.
##
## How the symbol is found. The two kinds of gene carry different naming schemes: named genes have
## HGNC-style products ("multiple EGF like domains 6"), LOC models UniProt-style protein names
## ("multiple epidermal growth factor-like domains protein 6"). So the LOC product is matched against
## the recommended protein name of human Swiss-Prot entries, which is the scheme it was drawn from.
## Match is exact after normalization (case, isoform and transcript-variant suffixes, "homolog",
## "LOW QUALITY PROTEIN"); a trailing "-like" is also tried stripped.
##
## Conservative on purpose: a recommended name shared by several human genes gives no alias, and a LOC
## gene whose products point to different symbols gives no alias. "Uncharacterized protein" models
## have nothing to match. The alias is a label, not an orthology call -- it says only that the model
## is annotated like the named human protein -- which is why it is displayed as "LOC... (SYMBOL-like)".
##
## Inputs: GTF_LONSTR_NCBI (config/paths.R) and a UniProt snapshot of human reviewed entries,
## data/orthologs/uniprot_human_reviewed.tsv, fetched with
##   curl --compressed -o data/orthologs/uniprot_human_reviewed.tsv \
##     "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606%20AND%20reviewed:true&fields=gene_primary,protein_name&format=tsv"
## The snapshot is bundled so the table can be rebuilt without depending on the current release.

suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here::here("config/paths.R"))

uniprot_fname = here::here("data/orthologs/uniprot_human_reviewed.tsv")
out_fname = here::here("data/orthologs/loc_human_aliases.csv")

norm_name = function(x) x %>% tolower() %>%
  str_remove("^low quality protein: ") %>%
  str_remove(", transcript variant .*$") %>%
  str_remove(",? ?isoform x?[0-9a-z]+$") %>%
  str_remove(" homolog$") %>%
  str_squish()

## Gene -> product pairs from the transcript records; one gene can carry several products.
gtf = read_tsv(GTF_LONSTR_NCBI, comment = "#", col_names = FALSE, col_types = cols(.default = "c"),
               progress = FALSE) %>%
  filter(X3 == "transcript", str_detect(X9, 'product "'))
gene_product = tibble(gene = str_match(gtf$X9, 'gene "([^"]+)"')[, 2],
                      product = str_match(gtf$X9, 'product "([^"]+)"')[, 2]) %>%
  filter(str_detect(gene, "^LOC")) %>%
  distinct()

uniprot = read_tsv(uniprot_fname, show_col_types = FALSE) %>%
  setNames(c("symbol", "names")) %>%
  filter(!is.na(symbol)) %>%
  mutate(rec = norm_name(str_remove(names, " \\(.*$"))) %>%
  distinct(symbol, rec) %>%
  group_by(rec) %>% filter(n_distinct(symbol) == 1) %>% ungroup()

aliases = gene_product %>%
  mutate(p = norm_name(product), p_nolike = str_remove(p, "-like( protein)?$")) %>%
  left_join(uniprot %>% rename(sym_exact = symbol), by = c("p" = "rec")) %>%
  left_join(uniprot %>% rename(sym_nolike = symbol), by = c("p_nolike" = "rec")) %>%
  mutate(symbol = coalesce(sym_exact, sym_nolike),
         match = if_else(!is.na(sym_exact), "exact", if_else(!is.na(sym_nolike), "-like stripped", NA))) %>%
  group_by(gene) %>%
  filter(n_distinct(na.omit(symbol)) == 1) %>%
  summarise(product = first(product[!is.na(symbol)]),
            human_symbol = first(na.omit(symbol)),
            match = first(match[!is.na(symbol)]),
            .groups = "drop") %>%
  arrange(gene)

write_csv(aliases, out_fname)
cat(sprintf("LOC genes with a product: %s; aliased: %s -> %s\n",
            n_distinct(gene_product$gene), nrow(aliases), out_fname))
