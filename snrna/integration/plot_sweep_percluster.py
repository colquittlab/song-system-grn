"""Per-cluster titration curves across a k.anchor sweep, one panel per sweep variant.

Each finch cluster's own best score is tracked across k.anchor, plotted as an individual
line (thin, low-alpha) so within-group heterogeneity is visible rather than hidden behind
a group mean. Group means are overlaid as thick lines.

Score = reciprocity-weighted: sqrt(fwd[i,j] * rev[j,i]), fwd = chicken-reference/finch-query
(the manuscript's own scoring direction), rev = finch-reference/chicken-query. A pair only
scores high if BOTH directions independently nominate it. This matters because the raw
one-directional fwd score can be artifactually high for a genuinely isolated query cluster
(no reject/novel option in TransferData's vote) -- observed for Glut-RA under the manuscript-
exact RPCA/SCT config: fwd score 0.968 (highest of any cluster) despite Glut-RA being, per
the source UMAP, the most non-overlapping cluster in the whole integration; the reverse
direction shows 374/377 (99%) of its forward-assigned partner's own cells vote to Arco
clusters, NONE to Glut-RA -- i.e. the fwd score alone was a one-directional vote artifact,
not evidence of a real match. Reciprocity weighting sends Glut-RA's score to 0.0.

Colour: validated song(red)/non-song(blue) pair (CVD 21.6, normal-vision 32.3, all-pairs,
light surface) -- reused from plot_rank_heatmap.py's GROUP_COLORS.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from assemble_glutonly import song_group, EXCLUDED_CLUSTERS
from reciprocal_score import best_per_cluster

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
GROUP_COLORS = {"song": "#e34948", "non-song": "#2a78d6"}
KA = [2, 5, 10, 20, 30, 50]
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly")


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
    ("cca_kfNA", "CCA\n(log-norm, default k.filter=NA)",
     "gg_glutonly_cca_matrix_finch_from_mouse_ka{ka}.csv",
     "gg_glutonly_cca_matrix_mouse_from_finch_ka{ka}.csv"),
    ("rpca_kfNA", "RPCA\n(log-norm, k.filter=NA)",
     "gg_glutonly_rpca_matrix_finch_from_mouse_ka{ka}.csv",
     "gg_glutonly_rpca_matrix_mouse_from_finch_ka{ka}.csv"),
    ("rpca_kf200", "RPCA\n(log-norm, k.filter=200)",
     "gg_glutonly_rpca_matrix_finch_from_mouse_ka{ka}_kf200.csv",
     "gg_glutonly_rpca_matrix_mouse_from_finch_ka{ka}_kf200.csv"),
    ("rpca_SCT_kf50", "RPCA\n(SCTransform, k.filter=50, dims=40 -- manuscript config)",
     "gg_glutonly_rpca_SCT_finch_from_mouse_ka{ka}_kf50_d40_matrix.csv",
     "gg_glutonly_rpca_SCT_mouse_from_finch_ka{ka}_kf50_d40_matrix.csv"),
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
    ax.set_xlabel("k.anchor  (weak → strong integration)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)

    for ka in present_ka:
        song = D.loc[[c for c in D.index if song_group(c) == "song"], ka]
        nonsong = D.loc[[c for c in D.index if song_group(c) == "non-song"], ka]
        u, p = mannwhitneyu(song, nonsong, alternative="less")
        stats_rows.append({"variant": tag, "k_anchor": ka,
                            "song_mean": song.mean(), "nonsong_mean": nonsong.mean(),
                            "U": u, "p_one_sided": p})

axes[0].set_ylabel("reciprocity-weighted best score per cluster", fontsize=8)
axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
fig.suptitle("Per-cluster integration-strength titration, reciprocity-weighted\n"
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
