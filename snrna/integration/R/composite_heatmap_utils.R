## ---------------------------------------------------------------------------
## Shared helpers for the cross-species composite heatmaps
## (snrna/integration/plot_composite_heatmaps_hybrid.R).
##
## R port of finch-integration-toolkit's class_benchmark.py + plot_rank_heatmap.py
## (class inference, column filtering, dot layout). Kept in one place so the
## heatmap driver and any future composite figure classify labels identically.
## ---------------------------------------------------------------------------

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## --- Coarse class inference ---------------------------------------------------
## Finch cluster name -> expected class, from the naming scheme itself.
FINCH_EXPECT <- list(
  c("^PC-",                 "progenitor"),
  c("(-NB(-\\d+)?$)|IP",    "neuroblast"),
  c("^(GABA|Glut)-",        "neuron"),
  c("^Astro(-\\d+)?$",      "astro"),
  c("^(OPC|Oligo)(-\\d+)?$","oligo"),
  c("^Micro$",              "immune"),
  c("^Epen$|^Epen-",        "ependymal"),
  c("^ChP$",                "ependymal"),
  c("^Endo(-\\d+)?$",       "vascular")
)
expected_class <- function(x) {
  vapply(x, function(n) {
    for (p in FINCH_EXPECT) if (grepl(p[1], n, perl = TRUE)) return(p[2])
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

## Reference label -> class. The label's OWN name is tried first (it is specific:
## "Ependymal NN_1", "Astro-OLF NN_1"); the annotation text is consulted only if
## the name alone is uninformative. Order matters: 'immune' before 'ependymal'
## before 'astro' so Yao's combined "Astro-Epen" class string cannot misfire.
REF_PAT <- c(
  progenitor = "Radial glia|\\bRgl|Neural progenitor|\\bNPC",
  neuroblast = "Neuroblast|\\bNbl|IMN\\b",
  neuron     = "Neuron|Gaba|Glut|Chol|Dopa|\\bNeur\\d|IT |ET |CTX|Sst|Pvalb|Vip|Lamp5|MSN|D1|D2",
  astro      = "Astro|Glioblast|\\bGbl",
  oligo      = "Oligo|\\bOPC|COP |MOL|NFOL|Committed oligodendrocyte",
  immune     = "Immune|Microglia|\\bMgl|\\bPvm|BAM |DC NN|Macrophage",
  ependymal  = "Ependymal|Epen\\d|^Epen|Tanycyte|CHOR|Chpl|Choroid|Hypendymal",
  vascular   = "Vascular|Endo|VLMC VLMC|VLMC|Peri|Peric|Pia|Meninges|Arachnoid|Dura|Fibro|Vendo|Angiob|ABC NN|Mesenchyme|SMC|Vsm"
)
REF_ORDER <- c("immune", "ependymal", "astro", "oligo", "vascular", "progenitor",
               "neuroblast", "neuron")
ref_class <- function(labels, annot = NULL) {
  vapply(labels, function(lab) {
    texts <- lab
    if (!is.null(annot) && lab %in% rownames(annot)) {
      texts <- c(lab, paste(lab, paste(as.character(unlist(annot[lab, ])), collapse = " ")))
    }
    for (t in texts) for (cls in REF_ORDER)
      if (grepl(REF_PAT[[cls]], t, perl = TRUE, ignore.case = TRUE)) return(cls)
    "unknown"
  }, character(1), USE.NAMES = FALSE)
}

## 6-class display merge; "neuron" is then split glut/gaba by name pattern.
MERGE <- c(neuron = "neuron", neuroblast = "neuroblast", progenitor = "progenitor",
           astro = "glia", oligo = "glia", ependymal = "glia",
           immune = "immune", vascular = "vascular", unknown = "unknown")
CLASS_ORDER <- c("glut", "gaba", "neuroblast", "progenitor", "glia", "immune", "vascular")
## Validated all-pairs on the light surface (normal-vision floor 15.6, CVD worst 6.9,
## legal with the legend + strip position as secondary encoding).
CLASS_COLORS <- c(glut = "#2a78d6", gaba = "#a34e9e", neuroblast = "#1baf7a",
                  progenitor = "#eda100", glia = "#008300", immune = "#4a3aa7",
                  vascular = "#e34948", unknown = "#c9c9c4", neuron = "#2a78d6")

glut_gaba_split <- function(names, fallback) {
  n <- tolower(names)
  out <- fallback
  is_gaba <- grepl("^gaba|^inh_|^inh-|gaba", n)
  is_glut <- !is_gaba & grepl("^glut|^ex_|^ex-|glut", n)
  out[is_gaba] <- "gaba"; out[is_glut] <- "glut"
  out
}
display_class <- function(names, coarse) {
  cls <- unname(MERGE[ifelse(is.na(coarse), "unknown", coarse)])
  cls[is.na(cls)] <- "unknown"
  out <- ifelse(cls == "neuron", glut_gaba_split(names, cls), cls)
  ## a neuron whose name says neither (e.g. a cholinergic/Sncg label) falls back to
  ## the glut slot, as the Python tool did via its "neuron" colour alias
  out[out == "neuron"] <- "glut"
  out
}

## --- Column filtering ---------------------------------------------------------
## Keep the union of each row's top-k columns, plus each row's peak-agreement
## column (so a column is never dropped while it still holds a row's best
## cross-method consensus); fall back to global strength if over max_cols, never
## trimming an anchor.
filter_columns <- function(M, top_k, max_cols, dots = NULL, dots_min = 2) {
  keep <- unique(unlist(lapply(seq_len(nrow(M)), function(i) {
    r <- M[i, ]; names(sort(r, decreasing = TRUE))[seq_len(min(top_k, length(r)))]
  })))
  anchors <- character(0)
  if (!is.null(dots)) {
    for (i in seq_len(nrow(dots))) {
      r <- dots[i, ]
      if (max(r) >= dots_min) anchors <- c(anchors, names(r)[which.max(r)])
    }
    anchors <- unique(anchors)
  }
  keep <- union(keep, anchors)
  if (length(keep) > max_cols) {
    trimmable <- setdiff(keep, anchors)
    n_keep <- max(max_cols - length(anchors), 0)
    strength <- sort(apply(M[, trimmable, drop = FALSE], 2, max), decreasing = TRUE)
    keep <- union(anchors, names(strength)[seq_len(min(n_keep, length(strength)))])
  }
  colnames(M)[colnames(M) %in% keep]   # preserve original column order
}

## --- Dot layout (Zaremba convention: one dot per agreeing method) -------------
## Offsets in units of one cell ([-0.5, 0.5] per axis).
dot_offsets <- function(n) {
  if (n <= 0) return(matrix(numeric(0), ncol = 2))
  if (n == 1) return(matrix(c(0, 0), ncol = 2))
  if (n == 2) return(cbind(c(-0.16, 0.16), 0))
  if (n == 3) return(cbind(seq(-0.30, 0.30, length.out = 3), 0))
  nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
  g <- expand.grid(x = seq(-0.20, 0.20, length.out = nc), y = seq(-0.20, 0.20, length.out = nr))
  as.matrix(g[seq_len(n), ])
}

read_matrix <- function(path) {
  M <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  M[is.na(M)] <- 0
  M
}
