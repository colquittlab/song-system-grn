#!/usr/bin/env python3
"""Figure 2, panel G -- joint UMAP of the finch and chicken glutamatergic cells after RPCA integration (k.anchor = 20).

UMAP of the integrated embedding of the k.anchor = 20 RPCA run of panels E and F (Seurat FindTransferAnchors + IntegrateEmbeddings
on SCTransform data, 40 dimensions, then RunUMAP; finch Glut-SATB2-1 kept, Glut-DACH2-HVCra-Int left out; chicken Ex_Pre_KCNH7
and Ex_Pre_SATB2 removed). Left: finch cells coloured by the dissected region of their type (Glut-SATB2-1 in lilac), song types
labelled at their median position, chicken cells grey. Right: chicken cells coloured by the marker family in their type's name
(SATB2/KIAA1217, DACH2, CACNA1H, other), finch cells grey.
Same drawing as make_finch_chicken_class_figure() for plusSATB2_noPre, k.anchor = 20, in
snrna/integration/gg_glutonly_hybrid/rpca/plot_rpca_umap_satb2_tests.py, which wrote
  snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid/rpca_umap_finchpos_chickenclass_plusSATB2_noPre_ka20_kf200_d40_SCT.{pdf,png}
In the figure the titles and legends were redrawn by hand, the "Meso" label dropped and the song labels moved.

Input: the joint-UMAP coordinates of that run (common.RPCA_RESULTS/gg_glutonly_hybrid_plusSATB2_noPre_rpca_ka20_kf200_d40_SCT_umap.csv).
Output: FIG2_OUT/panel_G_rpca_umap_ka20.{pdf,png} and panel_G_rpca_umap_ka20.csv (plotted cells with their colour group)

  python figures/figure2/panel_G_rpca_umap.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
KA = 20
GREY = "#d9d9d9"
SONG_LABELS = {"Glut-DACH2-HVCra": "HVCra", "Glut-DACH2-HVCx": "HVCx", "Glut-DACH2-LMANco": "LMANco",
               "Glut-DACH2-LMANsh": "LMANsh", "Glut-CACNA1H-RA": "RA"}
# Chicken types grouped by the marker family in their name, coloured to echo the finch palette (a visual pairing only)
CLASS_COLORS = {"SATB2/KIAA1217": c.LILAC, "DACH2": c.POSITION_COLORS["nr"], "CACNA1H": c.POSITION_COLORS["arco"], "other": "#b15928"}
CLASS_ORDER = ["SATB2/KIAA1217", "DACH2", "CACNA1H", "other"]


def chicken_class(label):
    if label.startswith("Ex_SATB2") or label.startswith("Ex_KIAA1217"):
        return "SATB2/KIAA1217"
    if label.startswith("Ex_DACH2"):
        return "DACH2"
    if label.startswith("Ex_CACNA1H"):
        return "CACNA1H"
    return "other"


def text_label(ax, sub, txt, color):
    ax.text(sub["UMAP_1"].median(), sub["UMAP_2"].median(), txt, fontsize=6.5, fontweight="bold", ha="center", va="center",
            color=color, zorder=3, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])


d = pd.read_csv(c.RPCA_RESULTS / f"{c.RUN_TAG}_rpca_ka{KA}_kf{c.KF}_d{c.DIMS}_SCT_umap.csv")
d = d[~d["label"].isin(c.EXCLUDED_CLUSTERS)].copy()
finch = d[d["species"] == "finch"].copy()
chick = d[d["species"] == "chicken"].copy()
assert set(finch["label"]) <= set(c.CLUSTER_POSITION) | {c.SATB2_CLUSTER}

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6))

# left: finch by dissected region, chicken grey
ax = axes[0]
ax.scatter(chick["UMAP_1"], chick["UMAP_2"], s=1.2, c=GREY, alpha=0.5, linewidths=0, rasterized=True, zorder=1)
present_positions = []
for cl in finch["label"].unique():   # drawn in order of first appearance in the table, as in the original
    if cl == c.SATB2_CLUSTER:
        continue
    pos = c.CLUSTER_POSITION[cl]
    present_positions.append(pos)
    sub = finch[finch["label"] == cl]
    ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=c.POSITION_COLORS[pos], alpha=0.6, linewidths=0, rasterized=True, zorder=2)
satb2 = finch[finch["label"] == c.SATB2_CLUSTER]
ax.scatter(satb2["UMAP_1"], satb2["UMAP_2"], s=1.6, c=c.LILAC, alpha=0.7, linewidths=0, rasterized=True, zorder=2)
ax.set_title("finch only, by anatomical position\n(chicken cells shown in grey)", fontsize=9)
handles = [Patch(facecolor=c.POSITION_COLORS[p], edgecolor="none", label=c.POSITION_LABELS[p])
           for p in c.POSITION_ORDER if p in set(present_positions)]
handles.append(Patch(facecolor=c.LILAC, edgecolor="none", label="Meso (=old NR-5)"))
ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)
for cl, short in SONG_LABELS.items():
    text_label(ax, finch[finch["label"] == cl], short, c.POSITION_COLORS[c.CLUSTER_POSITION[cl]])
text_label(ax, satb2, "Meso", c.LILAC)

# right: chicken by marker family, finch grey
ax = axes[1]
ax.scatter(finch["UMAP_1"], finch["UMAP_2"], s=1.0, c=GREY, alpha=0.4, linewidths=0, rasterized=True, zorder=1)
cls_of = chick["label"].map(chicken_class)
present = [k for k in CLASS_ORDER if k in set(cls_of)]
for k in present:
    sub = chick[cls_of == k]
    ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=CLASS_COLORS[k], alpha=0.6, linewidths=0, rasterized=True, zorder=2)
ax.set_title("chicken only, by excitatory class\n(finch cells shown in grey)", fontsize=9)
ax.legend(handles=[Patch(facecolor=CLASS_COLORS[k], edgecolor="none", label=k) for k in present],
          fontsize=6.5, frameon=False, loc="upper right")

for ax in axes:
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
fig.suptitle("Joint RPCA embedding (SCTransform), Glut(finch) x Excitatory(chicken)\n"
             "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed\n"
             f"k.anchor={KA}, k.filter={c.KF} (nominal only, forced NA under SCT), dims={c.DIMS} -- "
             "finch position vs. chicken class",
             fontsize=9, y=1.08)
fig.tight_layout()
out = c.FIG2_OUT / f"panel_G_rpca_umap_ka{KA}"
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

group = pd.Series("", index=d.index)
group[finch.index] = [c.SATB2_CLUSTER if l == c.SATB2_CLUSTER else c.CLUSTER_POSITION[l] for l in finch["label"]]
group[chick.index] = chick["label"].map(chicken_class)
d.assign(colour_group=group).to_csv(f"{out}.csv", index=False)
print(f"wrote {out}.pdf / .png / .csv")
