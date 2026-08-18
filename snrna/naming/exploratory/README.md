# Exploratory naming work — not production

These notebooks are the record of how the naming scheme was arrived at. They are **not** part of
the analysis chain: nothing outside this directory reads their outputs, and they should not be
re-rendered as part of producing figures or results. The production filtering and naming live in
`../hybrid_division_naming.qmd` (see `../README.md`).

- `zaremba_naming_validation.qmd` — tests whether the Zaremba et al. `Ex_<gene1>[_<gene2>]` naming
  grammar reproduces on this dataset.
- `glut_ex_style_naming.qmd` — applies that grammar to the glutamatergic clusters.
- `apply_composite_labels.qmd` — names every glutamatergic cluster after its specific
  chicken-integration match. Superseded by the hybrid scheme, which deliberately makes the weaker
  claim: the specific per-cluster chicken correspondence is not well enough supported for most
  clusters, and is wrong in principle for the song nuclei.

Each still reads the main clustered object read-only and writes only into its own subdirectory
here, so re-rendering one is harmless — it just is not part of the pipeline.
