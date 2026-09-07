"""Cluster x {max.features, k.score} heatmaps for the hybrid-labeled glut-only finch x
excitatory chicken RPCA sweep -- companions to plot_sweep_heatmap_glutonly_hybrid.py's
k.anchor view, exploring the two other axes exposed in cca_transfer.R (see that file's
docstrings for what each parameter does). Both axes are swept ONE-AT-A-TIME from the
Seurat-default baseline (k.anchor=30, max.features=200, k.score=30) -- per explicit user
request, not a full factorial with k.anchor or with each other.

max.features={100,200,400,800}, k.score={10,30,50,100}, each x k.filter in {NA, 200}. The
baseline value on each axis (max.features=200, k.score=30) reuses the SAME file as the
k.anchor=30 point in the k.anchor sweep -- these are the identical Seurat run, not a
separate one, so no extra compute was spent re-deriving that point.

Same conventions as plot_sweep_heatmap_glutonly_hybrid.py: fixed biological row order
(song nuclei HVCra/HVCx/LMANco/LMANsh/RA, then non-song DACH2, then non-song CACNA1H),
Glut-DACH2-HVCra-Int excluded, square cells, 6pt labels, reciprocity-weighted scores
(reciprocal_score.best_per_cluster), vmax shared across the two k.filter panels WITHIN
each axis's own figure (max.features and k.score are plotted as separate figures, each
with its own scale -- not cross-comparable to each other or to the k.anchor figure without
checking both scales).
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

import sys
sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from reciprocal_score import best_per_cluster
from position_palette import POSITION_COLORS, POSITION_LABELS, POSITION_ORDER, CLUSTER_POSITION

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
R = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/cca")
OUT = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid/rpca/results/gg_glutonly_hybrid")
TAG = "gg_glutonly_hybrid"
BASE_KA = 30

SONG_CLUSTERS = {"Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Int", "Glut-DACH2-HVCx",
                 "Glut-CACNA1H-RA", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh"}
EXCLUDED_CLUSTERS = {"Glut-DACH2-HVCra-Int"}
ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]


def song_group(c):
    return "song" if c in SONG_CLUSTERS else "non-song"


def file_suffix(mf, ks):
    """mf/ks suffixes appended ONLY when non-default (matches cca_transfer.R's own
    filename convention) -- the baseline point (mf=200, ks=30) is the plain, unsuffixed
    ka30 file shared by both axes."""
    s = ""
    if mf != 200:
        s += f"_mf{mf}"
    if ks != 30:
        s += f"_ks{ks}"
    return s


def load_axis(axis_values, mf_for, ks_for, kf):
    """axis_values: values to sweep. mf_for/ks_for: functions(value) -> (mf, ks) for that point."""
    rows = {}
    for v in axis_values:
        mf, ks = mf_for(v), ks_for(v)
        suf = file_suffix(mf, ks)
        f = R / f"{TAG}_rpca_log_finch_from_mouse_ka{BASE_KA}_kf{kf}_d30{suf}_matrix.csv"
        rv = R / f"{TAG}_rpca_log_mouse_from_finch_ka{BASE_KA}_kf{kf}_d30{suf}_matrix.csv"
        if not (f.exists() and rv.exists()):
            print(f"WARNING: missing {f.name} or {rv.name}, skipping {v}")
            continue
        best_score, _, _ = best_per_cluster(f, rv)
        rows[v] = best_score
    D = pd.DataFrame(rows)
    D.columns.name = None
    D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])
    return D


def make_figure(axis_name, axis_values, mf_for, ks_for, xlabel, out_stub):
    loaded = {kf: load_axis(axis_values, mf_for, ks_for, kf) for kf in ["NA", "200"]}
    vmax = max(D.values.max() for D in loaded.values())
    obs_min = min(D.values.min() for D in loaded.values())

    LABEL_PT = 6
    n_rows = len(ROW_ORDER)
    row_pitch_in = LABEL_PT / 72.0 * 1.15
    core_h = row_pitch_in * n_rows
    fig_h = core_h + 1.55
    n_cols = len(axis_values)
    core_w = row_pitch_in * n_cols
    fig_w = 2 * (core_w + 1.9) + 0.6

    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    stats_rows = []
    for ax, kf in zip(axes, ["NA", "200"]):
        D = loaded[kf]
        order = [c for c in ROW_ORDER if c in D.index]
        missing = set(D.index) - set(order)
        if missing:
            print(f"WARNING ({axis_name}, kf={kf}): clusters not in ROW_ORDER: {sorted(missing)}")
            order += sorted(missing)
        Do = D.loc[order]

        sns.heatmap(Do, ax=ax, cmap="Oranges", vmin=0, vmax=vmax, square=True,
                   linewidths=0.4, linecolor="white",
                   cbar_kws={"label": "reciprocity-weighted score"},
                   yticklabels=list(Do.index))
        for tick, c in zip(ax.get_yticklabels(), Do.index):
            tick.set_color(POSITION_COLORS[CLUSTER_POSITION[c]])
        ax.set_xlabel(xlabel, fontsize=LABEL_PT)
        ax.set_ylabel("")
        ax.set_title(f"RPCA (hybrid labels), log-norm, k.filter={kf}", fontsize=8)
        ax.tick_params(axis="both", labelsize=LABEL_PT, length=2, pad=1.5)
        ax.figure.axes[-1].tick_params(labelsize=LABEL_PT)

        for v in axis_values:
            if v not in D.columns:
                continue
            song_vals = D.loc[[c for c in D.index if song_group(c) == "song"], v]
            nonsong_vals = D.loc[[c for c in D.index if song_group(c) == "non-song"], v]
            u, p = mannwhitneyu(song_vals, nonsong_vals, alternative="less")
            stats_rows.append({"axis": axis_name, "k_filter": kf, "value": v,
                               "song_mean": song_vals.mean(), "nonsong_mean": nonsong_vals.mean(),
                               "U": u, "p_one_sided": p})

    from matplotlib.patches import Patch
    present = [p for p in POSITION_ORDER if p in {CLUSTER_POSITION[c] for c in ROW_ORDER}]
    handles = [Patch(facecolor=POSITION_COLORS[p], edgecolor="none", label=POSITION_LABELS[p]) for p in present]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=len(handles),
              frameon=False, fontsize=8, title="row label colour  (anatomical position)", title_fontsize=8)
    fig.suptitle(f"Per-cluster response to {axis_name}, reciprocity-weighted RPCA (hybrid labels)\n"
                f"Glut(finch) x Excitatory(chicken)  --  k.anchor fixed at baseline ({BASE_KA})  --  "
                "rows in fixed order: song nuclei (HVCra, HVCx, LMANco, LMANsh, RA), "
                "then non-song DACH2, then non-song CACNA1H\n"
                f"colour scale fixed 0-{vmax:.2f} across both panels for comparability "
                f"(observed range {obs_min:.2f}-{vmax:.2f})",
                fontsize=8, y=1.22)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.9)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{out_stub}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{out_stub}.pdf / .png  ({fig_w:.1f}x{fig_h:.1f} in)")

    stats = pd.DataFrame(stats_rows)
    stats.to_csv(OUT / f"{out_stub}_stats.csv", index=False)
    pd.set_option("display.width", 120)
    print(stats.to_string(index=False))
    for kf, D in loaded.items():
        Dout = D.copy()
        Dout["group"] = [song_group(c) for c in Dout.index]
        Dout.to_csv(OUT / f"{out_stub}_kf{kf}_reciprocal.csv")


# --- max.features sweep: mf in {100,200,400,800}, ks fixed at baseline (30) ---
make_figure(
    axis_name="max.features", axis_values=[100, 200, 400, 800],
    mf_for=lambda v: v, ks_for=lambda v: 30,
    xlabel="max.features  (anchor-filtering search-space cap)",
    out_stub="sweep_heatmap_maxfeatures_reciprocal",
)

# --- k.score sweep: ks in {10,30,50,100}, mf fixed at baseline (200) ---
make_figure(
    axis_name="k.score", axis_values=[10, 30, 50, 100],
    mf_for=lambda v: 200, ks_for=lambda v: v,
    xlabel="k.score  (anchor-scoring neighbourhood size)",
    out_stub="sweep_heatmap_kscore_reciprocal",
)
