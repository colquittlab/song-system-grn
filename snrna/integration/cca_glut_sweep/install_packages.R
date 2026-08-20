## Install packages required by the Slurm integration scripts
## (prep.R, transfer.R, aggregate.R, config.R).
## Run once on the cluster before submitting jobs:
##   Rscript snrna/integration/slurm/install_packages.R
##
## If the system R library is read-only, set a personal library first:
##   mkdir -p ~/R/library
##   echo 'R_LIBS_USER=~/R/library' >> ~/.Renviron

cat("R version:", as.character(getRversion()), "\n\n")

## try_install: run expr, then verify pkg_name is loadable.
## install.packages() only warns (never errors) on failure, so we must check.
try_install = function(label, pkg_name, expr) {
  cat(label, "... ")
  withCallingHandlers(
    tryCatch({
      expr
      if (requireNamespace(pkg_name, quietly = TRUE)) {
        cat("OK\n")
        TRUE
      } else {
        cat("FAILED (package not loadable after install)\n")
        FALSE
      }
    }, error = function(e) {
      cat("FAILED\n  ", conditionMessage(e), "\n")
      FALSE
    }),
    warning = function(w) {
      if (grepl("not available", conditionMessage(w)))
        cat("FAILED\n  ", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }
  )
}

## --- pak ---
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://r-lib.github.io/p/pak/stable/")
}

## --- CRAN (via pak) ---
## Seurat 5 requires Matrix >= 1.6.4 which requires R >= 4.4
seurat_pkg = if (getRversion() < "4.4.0") "Seurat@4.4.0" else "Seurat"

cran_pkgs = c(
  "tidyverse",
  seurat_pkg,
  "here",
  "Matrix",
  "RhpcBLASctl",
  "RcppParallel"
)

results = list()
for (pkg in cran_pkgs) {
  pkg_bare = sub("@.*", "", pkg)  # strip version suffix for requireNamespace check
  if (requireNamespace(pkg_bare, quietly = TRUE)) {
    cat(pkg_bare, "... already installed\n")
    results[[pkg_bare]] = TRUE
  } else {
    ok = try_install(pkg, pkg_bare, pak::pkg_install(pkg, ask = FALSE, upgrade = FALSE))
    results[[pkg_bare]] = ok
  }
}

## --- qs: no binary on conda-forge R; install from source ---
if (requireNamespace("qs", quietly = TRUE)) {
  cat("qs ... already installed\n")
  results[["qs"]] = TRUE
} else {
  ok = try_install("qs (pak)", "qs",
                   pak::pkg_install("qs", ask = FALSE, upgrade = FALSE))
  if (!ok) {
    ok = try_install("qs (source)", "qs",
                     install.packages("qs", type = "source",
                                      repos = "https://cloud.r-project.org"))
  }
  results[["qs"]] = ok
}

## --- Bioconductor ---
for (pkg in c("BiocParallel", "glmGamPoi")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "... already installed\n")
    results[[pkg]] = TRUE
  } else {
    ok = try_install(paste0("bioc::", pkg), pkg,
                     pak::pkg_install(paste0("bioc::", pkg), ask = FALSE, upgrade = FALSE))
    results[[pkg]] = ok
  }
}

## --- grr + Matrix.utils (both archived from CRAN) ---
## Always reinstall grr from source to overwrite any stale copy from an old R library.
## Install Matrix.utils via install.packages to bypass pak's dependency solver,
## which chokes on platform-mismatched grr copies in .libPaths().
ok = try_install("grr (CRAN archive)", "grr",
                 install.packages(
                   "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
                   repos = NULL, type = "source"))
results[["grr"]] = ok

ok = try_install("Matrix.utils (GitHub)", "Matrix.utils",
                 install.packages(
                   "https://github.com/cvarrichio/Matrix.utils/archive/refs/heads/master.tar.gz",
                   repos = NULL, type = "source"))
results[["Matrix.utils"]] = ok

## --- Summary ---
cat("\n=== Summary ===\n")
cat(sum(unlist(results)), "/", length(results), "packages installed successfully.\n")
failed = names(results)[!unlist(results)]
if (length(failed) > 0)
  cat("Failed:", paste(failed, collapse = ", "), "\n")
