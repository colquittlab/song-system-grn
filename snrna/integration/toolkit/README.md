# Vendored finch-integration-toolkit modules

Copies of the two modules `plot_composite_heatmaps_hybrid.py` needs from the lab's
`finch-integration-toolkit` (`/private/groups/colquittlab/finch-integration-toolkit`),
so the composite heatmaps can be regenerated on a machine that has only this repo.

| File | Purpose |
|------|---------|
| `plot_rank_heatmap.py` | Class-annotated clustermap of a finch x reference score matrix (CLI). |
| `class_benchmark.py` | Finch-cluster / reference-label coarse-class inference used for the colour strips. |

Provenance: toolkit commit `174236c` plus its then-uncommitted local edits to
`plot_rank_heatmap.py` (`--transpose`, `--scale`, `--no_dendrograms`, `--label_pt_override`,
`--row_order_csv`), copied 2026-09-07. The toolkit has no public remote yet; when it does,
replace this directory with a pinned dependency and delete these copies. Until then, any fix
made here should be mirrored back to the toolkit (and vice versa).

Python environment: `envs/integration_plots.yaml` at the repo root.
