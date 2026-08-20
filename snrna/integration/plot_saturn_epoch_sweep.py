"""SATURN integration-strength titration: training-epoch checkpoint sweep.

Analogous to Seurat's k.anchor and SAMap's crossK -- more metric-learning training epochs
means the triplet loss has had more opportunity to pull matching cross-species cell types
together (and push non-matches apart), so few epochs = weak/undertrained integration, many
epochs = strong/fully-trained integration. Unlike the other two methods this is measured
as checkpoints along a SINGLE training run (seed=0, N_MACRO=2000, N_HV=8000, matching the
manuscript composite's own SATURN config) rather than independent re-runs, since epoch count
is inherently a trajectory rather than a free top-level hyperparameter.

Score = top1_p from analyze_integration.py's KNN-weighted label transfer (same metric used
throughout this project for the SATURN method, e.g. assemble_glutonly.py's "saturn" matrix).

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
EP = [25, 50, 75, 100, 125, 150]
GGA = Path("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult-glutonly/analysis")
OUT = Path("/private/groups/colquittlab/saturn/zaremba_composite/results/composite_gg_glutonly")

rows = {}
for ep in EP:
    f = GGA / f"epoch_sweep_ep{ep}" / "finch_cluster_calls.csv"
    if not f.exists():
        print(f"missing: {f}")
        continue
    d = pd.read_csv(f, index_col="finch_cluster")
    rows[ep] = d["top1_p"]

D = pd.DataFrame(rows)
D.columns.name = "epoch"
D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
D["group"] = [song_group(c) for c in D.index]
D.to_csv(OUT / "saturn_epoch_sweep_percluster.csv")

present_ep = [e for e in EP if e in D.columns]

fig, ax = plt.subplots(figsize=(5.2, 4.4))
for cluster in D.index:
    grp = D.loc[cluster, "group"]
    ax.plot(present_ep, D.loc[cluster, present_ep], color=GROUP_COLORS[grp],
            alpha=0.35, linewidth=0.9, marker="o", markersize=2.5, zorder=2)
stats_rows = []
for grp, color in GROUP_COLORS.items():
    rows_ = D.index[D["group"] == grp]
    mean_curve = D.loc[rows_, present_ep].mean(axis=0)
    ax.plot(present_ep, mean_curve, color=color, linewidth=3.0, zorder=5,
            label=f"{grp} mean (n={len(rows_)})", solid_capstyle="round")
for ep in present_ep:
    song = D.loc[D.group == "song", ep]
    nonsong = D.loc[D.group == "non-song", ep]
    u, p = mannwhitneyu(song, nonsong, alternative="less")
    stats_rows.append({"epoch": ep, "song_mean": song.mean(), "nonsong_mean": nonsong.mean(),
                        "U": u, "p_one_sided": p})

ax.set_xlabel("training epoch  (weak → strong integration)", fontsize=9)
ax.set_ylabel("SATURN KNN-transfer top1 probability (best match per cluster)", fontsize=9)
ax.set_title("SATURN training-epoch titration\nGlut(finch) x Excitatory(chicken)", fontsize=10)
ax.tick_params(labelsize=8)
ax.set_xticks(present_ep)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.grid(axis="y", color="#e8e7e1", linewidth=0.6, zorder=0)
ax.legend(fontsize=8, frameon=False, loc="best")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"saturn_epoch_sweep_percluster.{ext}", dpi=220, bbox_inches="tight")
plt.close(fig)

stats = pd.DataFrame(stats_rows)
stats.to_csv(OUT / "saturn_epoch_sweep_stats.csv", index=False)
pd.set_option("display.width", 120)
print(D.sort_values(["group"] + present_ep, ascending=[True] + [False]*len(present_ep)))
print()
print(stats.to_string(index=False))
print(f"\nwrote {OUT}/saturn_epoch_sweep_percluster.pdf/.png, .csv, saturn_epoch_sweep_stats.csv")
