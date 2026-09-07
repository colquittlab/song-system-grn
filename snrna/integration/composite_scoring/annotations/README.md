# Reference-label annotation tables

Small lookup tables mapping each reference dataset's cluster label to its coarse class,
used by `toolkit/plot_rank_heatmap.py` (colour strips) and `toolkit/class_benchmark.py`
(class-level accuracy). Tracked so the composite heatmaps can be regenerated from the
repo alone.

| File | Reference | Original location |
|------|-----------|-------------------|
| `gg_adult_label_annotation.csv` | Zaremba et al. adult chicken snRNA-seq | `datasets/snrna-bf-adult_snrna-gg-adult/data/` (gitignored) |
| `yao_label_annotation.csv` | Yao et al. 2023 mouse ABC Atlas telencephalon | `/private/groups/colquittlab/saturn/snrna-bf-dev_snrna-yao2023/data/` |
