"""Cluster x k.anchor heatmap of the reciprocity-weighted RPCA sweep for the ORIGINAL
(non-hybrid) glut-only chicken comparison (composite_gg_glutonly) -- same visual
convention as plot_sweep_heatmap_glutonly_hybrid.py (fixed biological row order, square
cells, 6pt labels, anatomical-position row colours, shared vmax across panels), but
reading directly from the *_reciprocal.csv files assemble_glutonly.py's sweep already
computed and saved (sweep_percluster_rpca_kfNA_reciprocal.csv, ..._kf200_reciprocal.csv)
rather than recomputing from raw matrices.

Song/non-song membership and the excluded cluster (Glut-HVC-1a) match
assemble_glutonly.py's own song_group()/EXCLUDED_CLUSTERS exactly (already applied
upstream -- Glut-HVC-1a is simply absent from these CSVs). Anatomical position for row
colour is inferred from each cluster's own name prefix (Glut-HVC-* -> hvc, Glut-RA ->
ra, Glut-LMAN-* -> lman, Glut-NC-* -> nc, Glut-NR-* -> nr, Glut-Arco-* -> arco), reusing
position_palette.py's colour scheme -- note Glut-LMAN-1 is coloured 'lman' by this
naming-based rule but is NON-SONG by assemble_glutonly.py's song_group() (LMAN shell,
reclassified), same split already seen for the hybrid analysis's Glut-DACH2-LMANsh.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
from position_palette import POSITION_COLORS, POSITION_LABELS, POSITION_ORDER

OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly")
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

ROW_ORDER = [
    "Glut-HVC-1", "Glut-HVC-2", "Glut-RA", "Glut-LMAN-2",
    "Glut-LMAN-1",
    "Glut-NC-1", "Glut-NC-2", "Glut-NC-3", "Glut-NC-4",
    "Glut-NR-1", "Glut-NR-2", "Glut-NR-3", "Glut-NR-4", "Glut-NR-5",
    "Glut-Arco-1", "Glut-Arco-2", "Glut-Arco-3", "Glut-Arco-4", "Glut-Arco-5",
]


def cluster_pos(c):
    if c.startswith("Glut-HVC"):
        return "hvc"
    if c.startswith("Glut-RA"):
        return "ra"
    if c.startswith("Glut-LMAN"):
        return "lman"
    if c.startswith("Glut-NC"):
        return "nc"
    if c.startswith("Glut-NR"):
        return "nr"
    if c.startswith("Glut-Arco"):
        return "arco"
    raise ValueError(c)


VARIANTS = [
    ("rpca_kfNA", "RPCA (log-norm, k.filter=NA)"),
    ("rpca_kf200", "RPCA (log-norm, k.filter=200)"),
]

loaded = {}
for tag, title in VARIANTS:
    D = pd.read_csv(OUT / f"sweep_percluster_{tag}_reciprocal.csv", index_col=0)
    ka_cols = [c for c in D.columns if c != "group"]
    M = D[ka_cols].copy()
    M.columns = [int(c) for c in M.columns]
    loaded[tag] = M

vmax = max(M.values.max() for M in loaded.values())
KA = sorted(loaded[VARIANTS[0][0]].columns)

LABEL_PT = 6
n_rows = len(ROW_ORDER)
row_pitch_in = LABEL_PT / 72.0 * 1.15
core_h = row_pitch_in * n_rows
n_cols = len(KA)
core_w = row_pitch_in * n_cols
fig_w = 2 * (core_w + 1.9) + 0.6
fig_h = core_h + 1.55

fig, axes = plt.subplots(1, len(VARIANTS), figsize=(fig_w, fig_h))
for ax, (tag, title) in zip(axes, VARIANTS):
    M = loaded[tag]
    order = [c for c in ROW_ORDER if c in M.index]
    missing = set(M.index) - set(order)
    if missing:
        print(f"WARNING ({tag}): clusters not in ROW_ORDER: {sorted(missing)}")
        order += sorted(missing)
    Mo = M.loc[order, KA]

    sns.heatmap(Mo, ax=ax, cmap="Oranges", vmin=0, vmax=vmax, square=True,
               linewidths=0.4, linecolor="white", cbar_kws={"label": "reciprocity-weighted score"},
               yticklabels=list(Mo.index))
    for tick, c in zip(ax.get_yticklabels(), Mo.index):
        tick.set_color(POSITION_COLORS[cluster_pos(c)])
    ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=LABEL_PT)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=8)
    ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)
    ax.figure.axes[-1].tick_params(labelsize=LABEL_PT)

present = [p for p in POSITION_ORDER if p in {cluster_pos(c) for c in ROW_ORDER}]
handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=len(handles),
          frameon=False, fontsize=8, title="row label colour  (anatomical position)", title_fontsize=8)
obs_min = min(M.values.min() for M in loaded.values())
fig.suptitle("Per-cluster response to k.anchor, reciprocity-weighted RPCA (original, non-hybrid labels)\n"
            "Glut(finch) x Excitatory(chicken)\n"
            f"colour scale fixed 0-{vmax:.2f} across both panels for comparability "
            f"(observed range {obs_min:.2f}-{vmax:.2f})",
            fontsize=8, y=1.20)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"sweep_heatmap_reciprocal.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}/sweep_heatmap_reciprocal.pdf / .png")
