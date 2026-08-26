"""UMAPs of the joint RPCA embedding for the two SATB2/mesopallial follow-up subsets
(plusSATB2, noMeso) -- same conventions as plot_rpca_umap.py (4-panel + song/non-song
split figures), extended with a 5th panel that highlights the Ex_SATB2*/Ex_KIAA1217*
("mesopallial") chicken clusters specifically (same prefix definition used by
prepare_satb2_test_subsets.py's own noMeso filter).

Glut-SATB2-1 (plusSATB2's extra finch cluster, =old scheme's Glut-NR-5, barcode-confirmed
100% cell-identical) and the chicken mesopallial highlight are both drawn in a light
purple/lilac (#cab2d6) -- common_aesthetics.R position_colors' 'fieldl' slot, the light
half of the same Brewer-Paired purple pair whose dark half (#6a3d9a, 'ov') was already
used for this cluster/cell-type family in plot_sweep_heatmap_satb2_tests.py's legend, so
this stays within the analysis's existing colour scheme rather than introducing an
unrelated hue.

Also generates a standalone chicken-by-excitatory-class figure (make_chicken_class_figure)
-- chicken cells coloured by the transcription-factor family each cluster is named for
(SATB2/KIAA1217, DACH2, CACNA1H, or "other"), finch cells shown in grey -- the mirror
image of make_figure()'s panel 4 (finch by anatomical position, chicken in grey).

Only ka in {20, 30, 50, 75} have UMAP coordinates computed (matching plot_rpca_umap.py's
SCT loop) -- kf=200, dims=40, SCT-only (no log-norm equivalent exists for these two
follow-up subsets). k.filter=200 is nominal-only here: Seurat forces k.filter to NA
internally whenever normalization.method="SCT", regardless of what's passed (see
R/cca_transfer.R) -- kept in the filename purely for naming parallelism with the main
analysis's log-norm runs.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import pandas as pd

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from position_palette import (POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION,
                              EXCLUDED_CLUSTERS, song_group)

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

LILAC = "#cab2d6"
SATB2_CLUSTER = "Glut-SATB2-1"
MESO_PREFIXES = ("Ex_SATB2", "Ex_KIAA1217")

# Broad chicken excitatory "class" groupings by cluster-name prefix (the transcription
# factor family each cluster is named for) -- coloured to MATCH the finch anatomical-
# position palette (POSITION_COLORS) rather than assigning fresh hues, per explicit
# request: SATB2/KIAA1217 reuses the mesopallial lilac (same as finch's Meso/SATB2-1),
# DACH2 takes NR's light blue, CACNA1H takes Arco's pink. This is a purely visual
# pairing across the two species' plots, not a claim of biological correspondence
# between e.g. chicken CACNA1H and finch Arco. "other" (ungrouped: Ex_Pre_*,
# Ex_TSHZ2_*, standalone Ex_BCL6) keeps its own colour, unused elsewhere in either palette.
CLASS_COLORS = {
    "SATB2/KIAA1217": LILAC,
    "DACH2": POSITION_COLORS["nr"],
    "CACNA1H": POSITION_COLORS["arco"],
    "other": "#b15928",
}
CLASS_ORDER = ["SATB2/KIAA1217", "DACH2", "CACNA1H", "other"]


def chicken_class(label):
    if label.startswith("Ex_SATB2") or label.startswith("Ex_KIAA1217"):
        return "SATB2/KIAA1217"
    if label.startswith("Ex_DACH2"):
        return "DACH2"
    if label.startswith("Ex_CACNA1H"):
        return "CACNA1H"
    return "other"

SONG_LABELS = {
    "Glut-DACH2-HVCra": "HVCra",
    "Glut-DACH2-HVCx": "HVCx",
    "Glut-DACH2-LMANco": "LMANco",
    "Glut-DACH2-LMANsh": "LMANsh",
    "Glut-CACNA1H-RA": "RA",
}

R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid")
GREY = "#d9d9d9"
SPECIES_COLORS = {"chicken": "#b15928", "finch": "#2a78d6"}
KF, DIMS = 200, 40


def is_meso(label):
    return any(label.startswith(p) for p in MESO_PREFIXES)


def load(tag, ka):
    in_path = R / f"{tag}_rpca_ka{ka}_kf{KF}_d{DIMS}_SCT_umap.csv"
    if not in_path.exists():
        print(f"SKIP: {in_path} not found")
        return None
    d = pd.read_csv(in_path)
    return d[~d["label"].isin(EXCLUDED_CLUSTERS)].copy()


def draw_finch_position_panel(ax, finch, chick):
    """Finch cells coloured by anatomical position (Glut-SATB2-1 in lilac if present),
    chicken cells shown in grey -- shared by make_figure()'s panel 4 and
    make_finch_chicken_class_figure()."""
    has_satb2 = SATB2_CLUSTER in finch["label"].values
    ax.scatter(chick["UMAP_1"], chick["UMAP_2"], s=1.2, c=GREY, alpha=0.5, linewidths=0, rasterized=True, zorder=1)
    present_positions = []
    for c in finch["label"].unique():
        if c == SATB2_CLUSTER or c not in CLUSTER_POSITION:
            continue
        pos = CLUSTER_POSITION[c]
        present_positions.append(pos)
        sub = finch[finch["label"] == c]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=POSITION_COLORS[pos], alpha=0.6,
                  linewidths=0, rasterized=True, zorder=2)
    if has_satb2:
        sub = finch[finch["label"] == SATB2_CLUSTER]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=LILAC, alpha=0.7, linewidths=0, rasterized=True, zorder=2)
    ax.set_title("finch only, by anatomical position\n(chicken cells shown in grey)", fontsize=9)
    present = [p for p in POSITION_ORDER if p in set(present_positions)]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
    if has_satb2:
        handles.append(Patch(facecolor=LILAC, edgecolor="none", label="Meso (=old NR-5)"))
    ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)

    for c, short in SONG_LABELS.items():
        sub = finch[finch["label"] == c]
        if sub.empty:
            continue
        x, y = sub["UMAP_1"].median(), sub["UMAP_2"].median()
        ax.text(x, y, short, fontsize=6.5, fontweight="bold", ha="center", va="center",
               color=POSITION_COLORS[CLUSTER_POSITION[c]], zorder=3,
               path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
    if has_satb2:
        sub = finch[finch["label"] == SATB2_CLUSTER]
        x, y = sub["UMAP_1"].median(), sub["UMAP_2"].median()
        ax.text(x, y, "Meso", fontsize=6.5, fontweight="bold", ha="center", va="center",
               color=LILAC, zorder=3, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])


def draw_chicken_class_panel(ax, finch, chick):
    """Chicken cells coloured by broad excitatory class, finch cells shown in grey --
    shared by make_chicken_class_figure() and make_finch_chicken_class_figure()."""
    ax.scatter(finch["UMAP_1"], finch["UMAP_2"], s=1.0, c=GREY, alpha=0.4, linewidths=0, rasterized=True, zorder=1)
    cls_of = chick["label"].map(chicken_class)
    present = [c for c in CLASS_ORDER if c in set(cls_of)]
    for cls in present:
        sub = chick[cls_of == cls]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=CLASS_COLORS[cls], alpha=0.6,
                  linewidths=0, rasterized=True, zorder=2)
    ax.set_title("chicken only, by excitatory class\n(finch cells shown in grey)", fontsize=9)
    handles = [Patch(facecolor=CLASS_COLORS[c], edgecolor="none", label=c) for c in present]
    ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right")


def make_figure(tag, ka, title_suffix):
    d = load(tag, ka)
    if d is None:
        return
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()

    fig, axes = plt.subplots(1, 5, figsize=(21.5, 4.3))

    ax = axes[0]
    for sp, color in SPECIES_COLORS.items():
        sub = d[d["species"] == sp]
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.2, c=color, alpha=0.35, linewidths=0, rasterized=True, label=sp)
    ax.set_title("coloured by species\n(overlaid)", fontsize=9)
    ax.legend(markerscale=8, fontsize=7.5, frameon=False, loc="upper right")

    for ax, sp, sub in [(axes[1], "chicken", chick), (axes[2], "finch", finch)]:
        ax.scatter(d["UMAP_1"], d["UMAP_2"], s=1.0, c=GREY, alpha=0.4, linewidths=0, rasterized=True, zorder=1)
        ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.4, c=SPECIES_COLORS[sp], alpha=0.6,
                  linewidths=0, rasterized=True, zorder=2)
        ax.set_title(f"{sp} only\n(all cells shown in grey)", fontsize=9)

    draw_finch_position_panel(axes[3], finch, chick)

    ax = axes[4]
    ax.scatter(finch["UMAP_1"], finch["UMAP_2"], s=1.0, c=GREY, alpha=0.4, linewidths=0, rasterized=True, zorder=1)
    meso_mask = chick["label"].map(is_meso)
    n_meso = int(meso_mask.sum())
    ax.scatter(chick.loc[~meso_mask, "UMAP_1"], chick.loc[~meso_mask, "UMAP_2"], s=1.2,
              c=SPECIES_COLORS["chicken"], alpha=0.35, linewidths=0, rasterized=True, zorder=2, label="other chicken Ex")
    if n_meso:
        ax.scatter(chick.loc[meso_mask, "UMAP_1"], chick.loc[meso_mask, "UMAP_2"], s=1.8,
                  c=LILAC, alpha=0.85, linewidths=0, rasterized=True, zorder=3,
                  label="mesopallial\n(Ex_SATB2*/Ex_KIAA1217*)")
    ax.set_title(f"chicken mesopallial highlight\n(finch cells shown in grey, n={n_meso})", fontsize=9)
    ax.legend(markerscale=6, fontsize=6.5, frameon=False, loc="upper right")

    for ax in axes:
        ax.set_aspect("equal", adjustable="datalim")
        ax.axis("off")

    fig.suptitle("Joint RPCA embedding (SCTransform), Glut(finch) x Excitatory(chicken)\n"
                f"{title_suffix}\n"
                f"k.anchor={ka}, k.filter={KF} (nominal only, forced NA under SCT), dims={DIMS}",
                fontsize=9.5, y=1.10)
    fig.tight_layout()
    stub = f"rpca_umap_{tag.replace('gg_glutonly_hybrid_', '')}_ka{ka}_kf{KF}_d{DIMS}_SCT"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


def make_chicken_class_figure(tag, ka, title_suffix):
    """Single-panel companion to make_figure()'s panel 4 (finch by anatomical position),
    mirrored for chicken: chicken cells coloured by broad excitatory class (the
    transcription-factor family each cluster is named for), finch cells shown in grey."""
    d = load(tag, ka)
    if d is None:
        return
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    draw_chicken_class_panel(ax, finch, chick)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")

    fig.suptitle("Joint RPCA embedding (SCTransform), Glut(finch) x Excitatory(chicken)\n"
                f"{title_suffix}\n"
                f"k.anchor={ka}, k.filter={KF} (nominal only, forced NA under SCT), dims={DIMS}",
                fontsize=8.5, y=1.06)
    fig.tight_layout()
    stub = f"rpca_umap_chickenclass_{tag.replace('gg_glutonly_hybrid_', '')}_ka{ka}_kf{KF}_d{DIMS}_SCT"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


def make_finch_chicken_class_figure(tag, ka, title_suffix):
    """Two-panel figure pairing make_figure()'s finch-by-anatomical-position panel with
    make_chicken_class_figure()'s chicken-by-excitatory-class panel, side by side, so
    each species' own sub-population layout can be compared directly in one image."""
    d = load(tag, ka)
    if d is None:
        return
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6))
    draw_finch_position_panel(axes[0], finch, chick)
    draw_chicken_class_panel(axes[1], finch, chick)
    for ax in axes:
        ax.set_aspect("equal", adjustable="datalim")
        ax.axis("off")

    fig.suptitle("Joint RPCA embedding (SCTransform), Glut(finch) x Excitatory(chicken)\n"
                f"{title_suffix}\n"
                f"k.anchor={ka}, k.filter={KF} (nominal only, forced NA under SCT), dims={DIMS} -- "
                "finch position vs. chicken class",
                fontsize=9, y=1.08)
    fig.tight_layout()
    stub = f"rpca_umap_finchpos_chickenclass_{tag.replace('gg_glutonly_hybrid_', '')}_ka{ka}_kf{KF}_d{DIMS}_SCT"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


def make_song_split_figure(tag, ka, title_suffix):
    d = load(tag, ka)
    if d is None:
        return
    finch = d[d["species"] == "finch"].copy()
    chick = d[d["species"] == "chicken"].copy()
    has_satb2 = SATB2_CLUSTER in finch["label"].values

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
    for ax, grp, labels in [(axes[0], "song", SONG_LABELS), (axes[1], "non-song", {})]:
        ax.scatter(chick["UMAP_1"], chick["UMAP_2"], s=1.2, c=GREY, alpha=0.5, linewidths=0, rasterized=True, zorder=1)
        sub_grp = finch[(finch["label"] != SATB2_CLUSTER) & (finch["label"].map(song_group) == grp)]
        present_positions = []
        for c in sub_grp["label"].unique():
            pos = CLUSTER_POSITION[c]
            present_positions.append(pos)
            sub = sub_grp[sub_grp["label"] == c]
            ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=POSITION_COLORS[pos], alpha=0.6,
                      linewidths=0, rasterized=True, zorder=2)
        if grp == "non-song" and has_satb2:
            sub = finch[finch["label"] == SATB2_CLUSTER]
            ax.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.6, c=LILAC, alpha=0.7, linewidths=0, rasterized=True, zorder=2)
            x, y = sub["UMAP_1"].median(), sub["UMAP_2"].median()
            ax.text(x, y, "Meso", fontsize=6.5, fontweight="bold", ha="center", va="center",
                   color=LILAC, zorder=3, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
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
        if grp == "non-song" and has_satb2:
            handles.append(Patch(facecolor=LILAC, edgecolor="none", label="Meso (=old NR-5)"))
        ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)
        ax.set_aspect("equal", adjustable="datalim")
        ax.axis("off")

    fig.suptitle("Joint RPCA embedding (SCTransform), Glut(finch) x Excitatory(chicken)\n"
                f"{title_suffix}\n"
                f"k.anchor={ka}, k.filter={KF} (nominal only), dims={DIMS} -- song vs. non-song split",
                fontsize=9.5, y=1.08)
    fig.tight_layout()
    stub = f"rpca_umap_songsplit_{tag.replace('gg_glutonly_hybrid_', '')}_ka{ka}_kf{KF}_d{DIMS}_SCT"
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stub}.pdf / .png")


TESTS = [
    ("gg_glutonly_hybrid_plusSATB2", "Glut-SATB2-1 (=old Glut-NR-5) added back to the finch subset"),
    ("gg_glutonly_hybrid_noMeso", "chicken Ex_SATB2*/Ex_KIAA1217* mesopallial types removed"),
    ("gg_glutonly_hybrid_plusSATB2_noPre",
     "Glut-SATB2-1 added back to the finch subset; chicken Ex_Pre_KCNH7/Ex_Pre_SATB2 removed"),
]

if __name__ == "__main__":
    for tag, title_suffix in TESTS:
        for ka in [20, 30, 50, 75]:
            make_figure(tag, ka, title_suffix)
            make_song_split_figure(tag, ka, title_suffix)
            make_chicken_class_figure(tag, ka, title_suffix)
            make_finch_chicken_class_figure(tag, ka, title_suffix)
