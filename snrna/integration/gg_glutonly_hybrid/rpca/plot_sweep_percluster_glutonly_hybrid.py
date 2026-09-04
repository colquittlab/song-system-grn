"""Per-cluster RPCA k.anchor titration curves for the HYBRID-labeled glut-only finch x
excitatory chicken comparison -- repeat of plot_sweep_percluster.py's reciprocity-weighted
scoring and plotting, using the gg_glutonly_hybrid RPCA sweep (log-norm, k.filter in
{NA, 200}, dims=30) in place of the original gg_glutonly tag. No CCA (reduction=cca) or
RPCA_SCT variant was run for the hybrid tag, so only the two plain-RPCA panels are drawn.

song_group() is redefined here (not imported from assemble_glutonly.py) because that
module's song/non-song membership is keyed to the OLD flat 'cluster' scheme (Glut-HVC-*,
Glut-RA-*, Glut-LMAN-*), which does not exist under the hybrid celltype_hybrid labels.
Song membership here matches the set used throughout this hybrid-label analysis: the six
position-restricted song-nucleus populations identified in bf_adult_hybrid.h5ad's own
labelling (see prepare_finch_adult_hybrid.py) -- Glut-DACH2-HVCra, -HVCra-Pre, -HVCx,
Glut-CACNA1H-RA, Glut-DACH2-LMANco, and Glut-DACH2-LMANsh (LMAN shell, reclassified song
per explicit user request -- overriding this project's earlier non-song call on LMAN
shell, e.g. the old scheme's Glut-LMAN-1).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from reciprocal_score import best_per_cluster

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
GROUP_COLORS = {"song": "#e34948", "non-song": "#2a78d6"}
KA = [2, 5, 10, 20, 30, 50, 75, 100]
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid")
OUT.mkdir(parents=True, exist_ok=True)

SONG_CLUSTERS = {"Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Pre", "Glut-DACH2-HVCx",
                 "Glut-CACNA1H-RA", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh"}
# HVCra-Pre is a putative HVCra precursor population (per bf_adult_hybrid.h5ad's own
# naming) -- an edge case excluded from the song-vs-non-song comparison at the user's
# request, same mechanism as assemble_glutonly.py's EXCLUDED_CLUSTERS.
EXCLUDED_CLUSTERS = {"Glut-DACH2-HVCra-Pre"}


def song_group(c):
    return "song" if c in SONG_CLUSTERS else "non-song"


def load_variant(fwd_pattern, rev_pattern):
    rows = {}
    for ka in KA:
        f = R / fwd_pattern.format(ka=ka)
        rv = R / rev_pattern.format(ka=ka)
        if not (f.exists() and rv.exists()):
            continue
        best_score, _, _ = best_per_cluster(f, rv)
        rows[ka] = best_score
    if not rows:
        return None
    D = pd.DataFrame(rows)
    D.columns.name = "k_anchor"
    D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
    return D


VARIANTS = [
    ("rpca_kfNA", "RPCA (hybrid labels)\n(log-norm, k.filter=NA)",
     "gg_glutonly_hybrid_rpca_log_finch_from_mouse_ka{ka}_kfNA_d30_matrix.csv",
     "gg_glutonly_hybrid_rpca_log_mouse_from_finch_ka{ka}_kfNA_d30_matrix.csv"),
    ("rpca_kf200", "RPCA (hybrid labels)\n(log-norm, k.filter=200)",
     "gg_glutonly_hybrid_rpca_log_finch_from_mouse_ka{ka}_kf200_d30_matrix.csv",
     "gg_glutonly_hybrid_rpca_log_mouse_from_finch_ka{ka}_kf200_d30_matrix.csv"),
]

datasets = [(tag, title, load_variant(fp, rp)) for tag, title, fp, rp in VARIANTS]
datasets = [(tag, t, d) for tag, t, d in datasets if d is not None]

fig, axes = plt.subplots(1, len(datasets), figsize=(4.2 * len(datasets), 4.2), sharey=True)
if len(datasets) == 1:
    axes = [axes]

stats_rows = []
for ax, (tag, title, D) in zip(axes, datasets):
    present_ka = list(D.columns)
    for cluster in D.index:
        grp = song_group(cluster)
        ax.plot(present_ka, D.loc[cluster, present_ka], color=GROUP_COLORS[grp],
                alpha=0.35, linewidth=0.9, marker="o", markersize=2.2, zorder=2)
    for grp, color in GROUP_COLORS.items():
        rows = [c for c in D.index if song_group(c) == grp]
        mean_curve = D.loc[rows, present_ka].mean(axis=0)
        ax.plot(present_ka, mean_curve, color=color, linewidth=3.0, zorder=5,
                label=f"{grp} mean (n={len(rows)})", solid_capstyle="round")
    ax.set_xscale("log")
    ax.set_xticks(present_ka); ax.set_xticklabels([str(k) for k in present_ka])
    ax.set_xlabel("k.anchor  (weak -> strong integration)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)

    for ka in present_ka:
        song_vals = D.loc[[c for c in D.index if song_group(c) == "song"], ka]
        nonsong_vals = D.loc[[c for c in D.index if song_group(c) == "non-song"], ka]
        u, p = mannwhitneyu(song_vals, nonsong_vals, alternative="less")
        stats_rows.append({"variant": tag, "k_anchor": ka,
                            "song_mean": song_vals.mean(), "nonsong_mean": nonsong_vals.mean(),
                            "U": u, "p_one_sided": p})

axes[0].set_ylabel("reciprocity-weighted best score per cluster", fontsize=8)
axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
fig.suptitle("Per-cluster integration-strength titration, reciprocity-weighted (HYBRID labels)\n"
            "Glut(finch) x Excitatory(chicken)", fontsize=10, y=1.06)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"sweep_percluster_reciprocal.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}/sweep_percluster_reciprocal.pdf / .png  ({len(datasets)} panels)")

for tag, title, D in datasets:
    Dout = D.copy()
    Dout["group"] = [song_group(c) for c in Dout.index]
    Dout = Dout.sort_values(["group"] + list(D.columns), ascending=[True] + [False] * len(D.columns))
    Dout.to_csv(OUT / f"sweep_percluster_{tag}_reciprocal.csv")
    print(f"  wrote {OUT}/sweep_percluster_{tag}_reciprocal.csv")

stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "sweep_percluster_reciprocal_stats.csv", index=False)
pd.set_option("display.width", 120)
print()
print(stats.to_string(index=False))
