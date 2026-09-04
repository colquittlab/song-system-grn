"""SAMap integration-strength titration: NUMITERS (graph-reweighting iterations) sweep.

Analogous to the Seurat k.anchor sweep -- NUMITERS controls how many rounds SAMap spends
reweighting the cross-species gene-gene homology graph using expression correlation in
already-mapped cells, so low NUMITERS = weak/conservative integration, high NUMITERS =
more aggressively refined integration. Unlike Seurat's TransferData vote, SAMap's alignment
score (mutual-nearest-neighbour edge fraction, via get_mapping_scores) is already inherently
bidirectional, so no separate reciprocity correction is layered on here.

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
NI = [1, 2, 3, 5, 8, 10]
SW = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-adult_gg-adult-glutonly/results_sweep")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep/results/gg_glutonly")

rows = {}
for ni in NI:
    f = SW / f"ni{ni}" / "samap_finch_calls.csv"
    if not f.exists():
        print(f"missing: {f}")
        continue
    d = pd.read_csv(f, index_col="finch_cluster")
    rows[ni] = d["samap_score"]

D = pd.DataFrame(rows)
D.columns.name = "numiters"
D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
D["group"] = [song_group(c) for c in D.index]
D.to_csv(OUT / "samap_sweep_percluster.csv")

present_ni = [n for n in NI if n in D.columns]

fig, ax = plt.subplots(figsize=(5.2, 4.4))
for cluster in D.index:
    grp = D.loc[cluster, "group"]
    ax.plot(present_ni, D.loc[cluster, present_ni], color=GROUP_COLORS[grp],
            alpha=0.35, linewidth=0.9, marker="o", markersize=2.5, zorder=2)
stats_rows = []
for grp, color in GROUP_COLORS.items():
    rows_ = D.index[D["group"] == grp]
    mean_curve = D.loc[rows_, present_ni].mean(axis=0)
    ax.plot(present_ni, mean_curve, color=color, linewidth=3.0, zorder=5,
            label=f"{grp} mean (n={len(rows_)})", solid_capstyle="round")
for ni in present_ni:
    song = D.loc[D.group == "song", ni]
    nonsong = D.loc[D.group == "non-song", ni]
    u, p = mannwhitneyu(song, nonsong, alternative="less")
    stats_rows.append({"numiters": ni, "song_mean": song.mean(), "nonsong_mean": nonsong.mean(),
                        "U": u, "p_one_sided": p})

ax.set_xlabel("NUMITERS  (weak → strong integration)", fontsize=9)
ax.set_ylabel("SAMap alignment score (best match per cluster)", fontsize=9)
ax.set_title("SAMap integration-strength titration\nGlut(finch) x Excitatory(chicken)", fontsize=10)
ax.tick_params(labelsize=8)
ax.set_xticks(present_ni)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)
ax.legend(fontsize=8, frameon=False, loc="best")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"samap_sweep_percluster.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "samap_sweep_stats.csv", index=False)
pd.set_option("display.width", 120)
print(D.sort_values(["group"] + present_ni, ascending=[True] + [False]*len(present_ni)))
print()
print(stats.to_string(index=False))
print(f"\nwrote {OUT}/samap_sweep_percluster.pdf/.png, samap_sweep_percluster.csv, samap_sweep_stats.csv")
