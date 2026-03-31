## ── scenicplus.R ──────────────────────────────────────────────────────────────
## Helper functions for scenicplus analysis notebooks.

convert_category_to_factor <- function(x) {
  factor(x$codes, labels = x$categories, levels = 0:(length(x$categories) - 1))
}

extract_matrix <- function(obs_name, ereg_file) {
  cur    <- rhdf5::h5read(ereg_file, name = str_interp(obs_name))
  rnames <- cur[["var"]]$`_index`
  cnames <- cur[["obs"]]$Cell
  mat    <- data.frame(t(cur[["X"]]))
  colnames(mat) <- rnames
  rownames(mat) <- cnames
  mat
}
