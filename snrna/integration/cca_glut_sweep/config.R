## Shared parameter configuration for Slurm scripts.
## Must stay in sync with transfer_params in the QMD.

library(tidyverse)

transfer_params = expand_grid(
  k.anchor     = c(5),
  k.filter     = c(50),
  k.score      = c(30),
  max.features = c(100,200,400),
  dims_to_use  = list(1:25),
  reduction    = c("rpca")
)
