"""SAMap crossK sweep, rescored as a vote SHARE rather than SAMap's native mean-edge-weight.

SAMap's own alignment score is sum(cross-species edge weight to candidate) / crossK -- an
absolute mean over the crossK-neighbor window, which mechanically shrinks as crossK grows
regardless of whether one candidate still dominates (see plot_samap_crossk_sweep.py). To ask
the question that's actually analogous to Seurat's k.anchor sweep -- does widening the
cross-species neighborhood erode SPECIFIC discrimination, independent of the window-size
normalization artifact -- each finch cluster's row is renormalized to sum to 1 (a vote share,
exactly Seurat's prediction.score semantics) using the full alignment matrix already computed
per crossK run (samap_mapping_table.csv), no new SAMap runs needed.

Excludes Glut-HVC-1a, matching the rest of this analysis (see assemble_glutonly.py).
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

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
GROUP_COLORS = {"song": "#e34948", "non-song": "#2a78d6"}
CK = [5, 10, 20, 30, 50, 100]
SW = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly/results_sweep")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly")

rows = {}
for ck in CK:
    f = SW / f"ck{ck}" / "samap_mapping_table.csv"
    if not f.exists():
        print(f"missing: {f}")
        continue
    MT = pd.read_csv(f, index_col=0)
    bf_rows = [i for i in MT.index if str(i).startswith("bf_")]
    gg_cols = [c for c in MT.columns if str(c).startswith("gg_")]
    sub = MT.loc[bf_rows, gg_cols]
    sub.index = [i[3:] for i in sub.index]
    sub.columns = [c[3:] for c in sub.columns]
    share = sub.div(sub.sum(axis=1), axis=0)  # renormalize each finch cluster's row to sum to 1
    rows[ck] = share.max(axis=1)

D = pd.DataFrame(rows)
D.columns.name = "crossK"
D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
D["group"] = [song_group(c) for c in D.index]
D.to_csv(OUT / "samap_crossk_voteshare_percluster.csv")

present_ck = [n for n in CK if n in D.columns]

fig, ax = plt.subplots(figsize=(5.2, 4.4))
for cluster in D.index:
    grp = D.loc[cluster, "group"]
    ax.plot(present_ck, D.loc[cluster, present_ck], color=GROUP_COLORS[grp],
            alpha=0.35, linewidth=0.9, marker="o", markersize=2.5, zorder=2)
stats_rows = []
for grp, color in GROUP_COLORS.items():
    rows_ = D.index[D["group"] == grp]
    mean_curve = D.loc[rows_, present_ck].mean(axis=0)
    ax.plot(present_ck, mean_curve, color=color, linewidth=3.0, zorder=5,
            label=f"{grp} mean (n={len(rows_)})", solid_capstyle="round")
for ck in present_ck:
    song = D.loc[D.group == "song", ck]
    nonsong = D.loc[D.group == "non-song", ck]
    u, p = mannwhitneyu(song, nonsong, alternative="less")
    stats_rows.append({"crossK": ck, "song_mean": song.mean(), "nonsong_mean": nonsong.mean(),
                        "U": u, "p_one_sided": p})

ax.set_xscale("log")
ax.set_xlabel("crossK  (weak → strong integration)", fontsize=9)
ax.set_ylabel("vote share of best match per cluster (row-renormalized)", fontsize=9)
ax.set_title("SAMap crossK titration, vote-share rescoring\nGlut(finch) x Excitatory(chicken)", fontsize=10)
ax.tick_params(labelsize=8)
ax.set_xticks(present_ck); ax.set_xticklabels([str(c) for c in present_ck])
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)
ax.legend(fontsize=8, frameon=False, loc="best")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"samap_crossk_voteshare_percluster.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "samap_crossk_voteshare_stats.csv", index=False)
pd.set_option("display.width", 120)
print(D.sort_values(["group"] + present_ck, ascending=[True] + [False]*len(present_ck)))
print()
print(stats.to_string(index=False))
print(f"\nwrote {OUT}/samap_crossk_voteshare_percluster.pdf/.png, .csv, samap_crossk_voteshare_stats.csv")
