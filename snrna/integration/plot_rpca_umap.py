"""UMAP(s) of the joint RPCA embedding (gg_glutonly_hybrid), produced by
R/cca_transfer_umap.R -- the first genuine shared low-dimensional embedding built for
this comparison (cca_transfer.R itself only ever computes label-transfer anchors/scores,
never a joint embedding to visualize).

Four panels per (k.anchor, dims) configuration:
  1. species overlay -- both species alpha-blended together, to see overall mixing.
  2. chicken only, highlighted against a grey backdrop of ALL cells -- shows exactly the
     territory chicken cells occupy without finch points competing for the same pixels.
  3. finch only, same treatment -- the territory finch cells occupy.
  4. finch only, coloured by anatomical position (position_palette.py's
     CLUSTER_POSITION/POSITION_COLORS, the same scheme used throughout this analysis's
     heatmaps), chicken cells shown in grey -- where each finch population specifically
     lands, not just "finch" as a whole.
Panels 2+3 are the direct answer to "split by species to better see the space occupied by
finch/chicken" -- panel 1 alone can't show this because dense overlap hides which species
is filling a given region once points are drawn on top of each other.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              SONG_CLUSTERS, EXCLUDED_CLUSTERS, song_group)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

# The 5 song-nucleus clusters labelled on-plot in the 4-panel figure's rightmost panel --
# same set used throughout this analysis's heatmaps/stats (Glut-DACH2-HVCra-Pre excluded
# as an edge case there, per earlier explicit request). The dedicated song/non-song split
# figure below labels all 6 SONG_CLUSTERS members instead, since HVCra-Pre genuinely
# belongs to the song group being displayed there and excluding it would leave an
# unexplained coloured-but-unlabelled population.
SONG_LABELS = {
    "Glut-DACH2-HVCra": "HVCra",
    "Glut-DACH2-HVCx": "HVCx",
    "Glut-DACH2-LMANco": "LMANco",
    "Glut-DACH2-LMANsh": "LMANsh",
    "Glut-CACNA1H-RA": "RA",
}

R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly_hybrid")
TAG = "gg_glutonly_hybrid"
GREY = "#d9d9d9"
SPECIES_COLORS = {"chicken": "#b15928", "finch": "#2a78d6"}


def make_figure(ka, kf, dims, norm="log"):
    norm_suffix = "_SCT" if norm == "SCT" else ""
    in_path = R / f"{TAG}_rpca_ka{ka}_kf{kf}_d{dims}{norm_suffix}_umap.csv"
    if not in_path.exists():
        print(f"SKIP: {in_path} not found")
        return
    d = pd.read_csv(in_path)
    d = d[~d["label"].isin(EXCLUDED_CLUSTERS)].copy()
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()
    unmatched = sorted(set(finch["label"]) - set(CLUSTER_POSITION))
    if unmatched:
        print(f"WARNING (ka={ka}): {len(unmatched)} finch labels not in CLUSTER_POSITION: {unmatched}")

    fig, axes = plt.subplots(1, 4, figsize=(17.5, 4.3))

    ax = axes[0]
    for sp, color in SPECIES_COLORS.items():
        sub = d[d["species"] == sp]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.2, c=color, alpha=0.35, linewidths=0, label=sp)
    ax.set_title("coloured by species\n(overlaid)", fontsize=9)
    ax.legend(markerscale=8, fontsize=7.5, frameon=False, loc="upper right")

    for ax, sp, other, sub in [
        (axes[1], "chicken", "finch", chick),
        (axes[2], "finch", "chicken", finch),
    ]:
        ax.scatter(d["UMAP_1"], d["UMAP_2"], s=1.0, c=GREY, alpha=0.4, linewidths=0, zorder=1)
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.4, c=SPECIES_COLORS[sp], alpha=0.6,
                  linewidths=0, zorder=2)
        ax.set_title(f"{sp} only\n(all cells shown in grey)", fontsize=9)

    ax = axes[3]
    ax.scatter(chick["UMAP_1"], chick["UMAP_2"], s=1.2, c=GREY, alpha=0.5, linewidths=0, zorder=1)
    present_positions = []
    for c in finch["label"].unique():
        if c not in CLUSTER_POSITION:
            continue
        pos = CLUSTER_POSITION[c]
        present_positions.append(pos)
        sub = finch[finch["label"] == c]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=POSITION_COLORS[pos], alpha=0.6,
                  linewidths=0, zorder=2)
    ax.set_title("finch only, by anatomical position\n(chicken cells shown in grey)", fontsize=9)
    present = [p for p in POSITION_ORDER if p in set(present_positions)]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
    ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)

    # On-plot text labels for the 5 song-nucleus clusters, placed at each cluster's own
    # median UMAP position -- a white-outlined stroke keeps them legible over the coloured
    # scatter without needing a separate callout/legend.
    for c, short in SONG_LABELS.items():
        sub = finch[finch["label"] == c]
        if sub.empty:
            continue
        x, y = sub["UMAP_1"].median(), sub["UMAP_2"].median()
        ax.text(x, y, short, fontsize=6.5, fontweight="bold", ha="center", va="center",
               color=POSITION_COLORS[CLUSTER_POSITION[c]], zorder=3,
               path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])

    for ax in axes:
        ax.set_xlabel("UMAP 1", fontsize=8)
        ax.set_ylabel("UMAP 2", fontsize=8)
        ax.tick_params(labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_aspect("equal", adjustable="datalim")

    norm_label = "SCTransform (SEPARATE from log-norm)" if norm == "SCT" else "log-norm"
    fig.suptitle("Joint RPCA embedding (hybrid labels), Glut(finch) x Excitatory(chicken)\n"
                f"k.anchor={ka}, k.filter={kf}, {norm_label}, dims={dims} -- IntegrateEmbeddings + UMAP",
                fontsize=9.5, y=1.08)
    fig.tight_layout()
    stub = f"rpca_umap_ka{ka}_kf{kf}_d{dims}{norm_suffix}"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


def make_song_split_figure(ka, kf, dims, norm="log"):
    """Variation on panel 4 of make_figure(): song and non-song finch clusters drawn in
    SEPARATE panels rather than sharing one, so the two groups' territories can be
    compared without their points competing for the same pixels."""
    norm_suffix = "_SCT" if norm == "SCT" else ""
    in_path = R / f"{TAG}_rpca_ka{ka}_kf{kf}_d{dims}{norm_suffix}_umap.csv"
    if not in_path.exists():
        print(f"SKIP: {in_path} not found")
        return
    d = pd.read_csv(in_path)
    d = d[~d["label"].isin(EXCLUDED_CLUSTERS)].copy()
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
    for ax, grp, labels in [(axes[0], "song", SONG_LABELS), (axes[1], "non-song", {})]:
        ax.scatter(chick["UMAP_1"], chick["UMAP_2"], s=1.2, c=GREY, alpha=0.5, linewidths=0, zorder=1)
        sub_grp = finch[finch["label"].map(song_group) == grp]
        present_positions = []
        for c in sub_grp["label"].unique():
            pos = CLUSTER_POSITION[c]
            present_positions.append(pos)
            sub = sub_grp[sub_grp["label"] == c]
            ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=POSITION_COLORS[pos], alpha=0.6,
                      linewidths=0, zorder=2)
        for c, short in labels.items():
            sub = finch[finch["label"] == c]
            if sub.empty:
                continue
            x, y = sub["UMAP_1"].median(), sub["UMAP_2"].median()
            ax.text(x, y, short, fontsize=6.5, fontweight="bold", ha="center", va="center",
                   color=POSITION_COLORS[CLUSTER_POSITION[c]], zorder=3,
                   path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
        ax.set_title(f"{grp} finch clusters only\n(chicken cells shown in grey)", fontsize=9)
        present = [p for p in POSITION_ORDER if p in set(present_positions)]
        handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
        ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)
        ax.set_xlabel("UMAP 1", fontsize=8)
        ax.set_ylabel("UMAP 2", fontsize=8)
        ax.tick_params(labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_aspect("equal", adjustable="datalim")

    norm_label = "SCTransform (SEPARATE from log-norm)" if norm == "SCT" else "log-norm"
    fig.suptitle("Joint RPCA embedding (hybrid labels), Glut(finch) x Excitatory(chicken)\n"
                f"k.anchor={ka}, k.filter={kf}, {norm_label}, dims={dims} -- song vs. non-song split",
                fontsize=9.5, y=1.08)
    fig.tight_layout()
    stub = f"rpca_umap_songsplit_ka{ka}_kf{kf}_d{dims}{norm_suffix}"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


if __name__ == "__main__":
    make_figure(ka=50, kf=200, dims=30)
    make_song_split_figure(ka=50, kf=200, dims=30)
    for ka in [20, 30, 50, 75]:
        make_figure(ka=ka, kf=200, dims=40)
        make_song_split_figure(ka=ka, kf=200, dims=40)
    # SCTransform variant (separate analysis, kept out of the log-norm loops above so a
    # re-run of this script never silently regenerates SCT figures unless explicitly
    # requested). Two k.filter values computed: 50 (this project's original "manuscript
    # config" for SCT) and 200 (added per explicit request, to match the k.filter used
    # throughout the log-norm analysis and make the two directly comparable).
    for kf in [50, 200]:
        for ka in [20, 30, 50, 75]:
            make_figure(ka=ka, kf=kf, dims=40, norm="SCT")
            make_song_split_figure(ka=ka, kf=kf, dims=40, norm="SCT")
