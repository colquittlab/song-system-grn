"""SAMap integration-strength titration: crossK (cross-species edges per cell) sweep.

crossK is the number of cross-species edges identified per cell when building the initial
mutual-nearest-neighbour graph -- the direct analogue of Seurat's k.anchor (small crossK =
few, tight cross-species candidate edges = weak/conservative; large crossK = many candidate
edges = aggressive/strong). Contrast with NUMITERS, which only refines an already-fixed
graph and (per plot_samap_sweep.py) barely moves the result -- crossK controls what edges
exist in the first place.

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
    f = SW / f"ck{ck}" / "samap_finch_calls.csv"
    if not f.exists():
        print(f"missing: {f}")
        continue
    d = pd.read_csv(f, index_col="finch_cluster")
    rows[ck] = d["samap_score"]

D = pd.DataFrame(rows)
D.columns.name = "crossK"
D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
D["group"] = [song_group(c) for c in D.index]
D.to_csv(OUT / "samap_crossk_sweep_percluster.csv")

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
ax.set_ylabel("SAMap alignment score (best match per cluster)", fontsize=9)
ax.set_title("SAMap crossK titration\nGlut(finch) x Excitatory(chicken)", fontsize=10)
ax.tick_params(labelsize=8)
ax.set_xticks(present_ck); ax.set_xticklabels([str(c) for c in present_ck])
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)
ax.legend(fontsize=8, frameon=False, loc="best")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"samap_crossk_sweep_percluster.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "samap_crossk_sweep_stats.csv", index=False)
pd.set_option("display.width", 120)
print(D.sort_values(["group"] + present_ck, ascending=[True] + [False]*len(present_ck)))
print()
print(stats.to_string(index=False))
print(f"\nwrote {OUT}/samap_crossk_sweep_percluster.pdf/.png, .csv, samap_crossk_sweep_stats.csv")
