#!/usr/bin/env python3
"""Inputs to MAFB (the TFs and regions that regulate it) in GABAergic neurons versus Glut-CACNA1H-RA / Glut-CACNA1H-1.

"Inputs to MAFB" = the eRegulon rows whose target gene is MAFB: for each row, the regulating TF, the region it acts through, and the
TF->gene / region->gene evidence. Read from the small outputs of the three subset runs (results/config1 and config2, direct + extended):

    GABA     ra-arco-hvc-nc_gaba_hybrid                (3,275 cells)
    RA_run2  glut-ra_glut-c1h-1_hybrid                 (745 cells; run3 samples excluded, as in the earlier run)
    RA_all   glut-ra_glut-c1h-1_allruns_hybrid         (1,966 cells; run3 included)

MAFB listed as a target of MAFB (TF2G correlation exactly 1.000) is the TF->gene self-link, not autoregulation, and is dropped (the perturbation notebook
drops it too). Outputs: two small CSVs next to this script, and mafb_inputs_gaba_vs_ra.{pdf,png} in OUT_DIR.

    python mafb_inputs_gaba_vs_ra.py [OUT_DIR]
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("ignore")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

HERE = os.path.dirname(os.path.abspath(__file__))
HB = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/ssd/rstudio/multiome/motor-pathway/scenicplus/mafb_inputs")
os.makedirs(OUT_DIR, exist_ok=True)
DATASETS = {"GABA": "ra-arco-hvc-nc_gaba_hybrid", "RA_run2": "glut-ra_glut-c1h-1_hybrid", "RA_all": "glut-ra_glut-c1h-1_allruns_hybrid"}
LABEL = {"GABA": "GABA", "RA_run2": "RA vs C1H-1\n(run3 excl.)", "RA_all": "RA vs C1H-1\n(all runs)"}
# repo figure colours (config/figure_theme.R)
BLUE, RED, MID, GRID, INK2, MUTED = "#2a78d6", "#e34948", "#f0efec", "#e1e0d9", "#52514e", "#898781"

# ---- MAFB input rows
cols = ["TF", "Gene", "Region", "regulation", "importance_TF2G", "rho_TF2G", "triplet_rank", "eRegulon_name", "importance_R2G", "rho_R2G"]
rows = []
for ds, d in DATASETS.items():
    for n in (1, 2):
        for kind, f in (("direct", "eRegulon_direct.tsv"), ("extended", "eRegulons_extended.tsv")):
            t = pd.read_table(f"{HB}{d}/results/config{n}/outs/{f}", usecols=cols)
            m = t[(t.Gene == "MAFB") & (t.TF != "MAFB")].copy()
            m["ds"], m["config"], m["kind"] = ds, n, kind
            rows.append(m)
A = pd.concat(rows, ignore_index=True)

# ---- TF table
tf_tab = (A.groupby(["ds", "TF"]).agg(configs=("config", lambda s: "".join(map(str, sorted(set(s))))), kinds=("kind", lambda s: "+".join(sorted(set(s)))),
          n_regions=("Region", "nunique"), median_rho_TF2G=("rho_TF2G", "median"), max_importance_TF2G=("importance_TF2G", "max"),
          best_triplet_rank=("triplet_rank", "min")).reset_index().sort_values(["ds", "best_triplet_rank"]))
tf_tab.round(3).to_csv(f"{HERE}/mafb_inputs_tf_table.csv", index=False)

# ---- region table (distance to the MAFB TSS; MAFB is on the minus strand: TSS = End; negative = upstream)
ga = pd.read_table(f"{HB}{DATASETS['GABA']}/scenicplus/genome_annotation.tsv")
mg = ga[ga.Gene == "MAFB"].iloc[0]
tss = int(mg.End) if mg.Strand == "-" else int(mg.Start)
sign = 1 if mg.Strand == "-" else -1
gene_end = sign * (tss - (int(mg.Start) if mg.Strand == "-" else int(mg.End)))   # gene body: 0 .. gene_end (kb offsets, downstream positive)
reg_rows = []
for ds in DATASETS:
    for r, sub in A[A.ds == ds].groupby("Region"):
        c, se = r.split(":"); s, e = map(int, se.split("-")); mid = (s + e) // 2
        reg_rows.append(dict(ds=ds, region=r, offset_from_TSS_bp=(tss - mid) if mg.Strand == "-" else (mid - tss), n_TF=sub.TF.nunique(),
                             TFs=",".join(sorted(set(sub.TF))), configs="".join(map(str, sorted(set(sub.config)))), median_rho_R2G=sub.rho_R2G.median()))
R = pd.DataFrame(reg_rows)
R["shared_with"] = R.apply(lambda x: ",".join(sorted(set(R[(R.region == x.region) & (R.ds != x.ds)].ds))), axis=1)
R.sort_values(["ds", "offset_from_TSS_bp"]).round(3).to_csv(f"{HERE}/mafb_inputs_regions.csv", index=False)

# ---- figure
fig = plt.figure(figsize=(6.4, 6.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.9], hspace=0.38)
ax = fig.add_subplot(gs[0])
ax.set_xscale("symlog", linthresh=15, linscale=1.0)
for y, ds in enumerate(DATASETS):
    s = R[R.ds == ds]
    shared = s.shared_with != ""
    for flag, face, edge in ((False, "white", MUTED), (True, BLUE, BLUE)):
        q = s[shared == flag]
        ax.scatter(q.offset_from_TSS_bp / 1000, [y] * len(q), s=18 + 22 * q.n_TF, facecolor=face, edgecolor=edge, linewidth=0.9, zorder=3,
                   label=("region used in >= 2 datasets" if flag else "region found in this dataset only") if y == 0 else None)
    ax.axhline(y, color=GRID, linewidth=0.6, zorder=1)
ax.axvspan(0, gene_end / 1000, color="#c9c8c2", alpha=0.7, zorder=0)
ax.axvline(0, color=INK2, linewidth=0.8, zorder=2)
ax.text(0, -0.78, "MAFB\n(TSS, gene body)", ha="center", va="top", fontsize=5, color=INK2)
ax.set_yticks(range(len(DATASETS)))
ax.set_yticklabels([LABEL[d] for d in DATASETS], fontsize=6)
ax.set_ylim(-1.1, len(DATASETS) - 0.4)
ax.invert_yaxis()
ticks = [-100, -50, -15, 0, 15, 50, 100]
ax.set_xticks(ticks); ax.set_xticklabels([f"{t:+d}" if t else "0" for t in ticks], fontsize=5)
ax.set_xlabel("position relative to the MAFB TSS, kb (negative = upstream; log scale beyond +/-15 kb)", fontsize=6)
ax.tick_params(axis="x", labelsize=5, length=2)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(fontsize=5, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.28), ncol=2, handletextpad=0.2, columnspacing=0.8)
ax.set_title("Regions linked to MAFB (marker size = number of regulating TFs)", fontsize=7, fontweight="bold", loc="left", pad=14)

# TF presence map: rows = TFs, columns = dataset x config; colour = sign of the TF-MAFB expression correlation (rho_TF2G)
ax2 = fig.add_subplot(gs[1])
colspec = [(ds, n) for ds in DATASETS for n in (1, 2)]
order = []
for ds in DATASETS:
    present = A[A.ds == ds].groupby("TF").size().index
    others = set(A[A.ds != ds].TF)
    order += [t for t in tf_tab[tf_tab.ds == ds].sort_values("best_triplet_rank").TF if t not in order]
M = pd.DataFrame(np.nan, index=order, columns=pd.MultiIndex.from_tuples(colspec))
for (ds, n), sub in A.groupby(["ds", "config"]):
    for tf, v in sub.groupby("TF").rho_TF2G.median().items():
        M.loc[tf, (ds, n)] = v
cmap = LinearSegmentedColormap.from_list("div", [BLUE, MID, RED])
im = ax2.imshow(M.to_numpy(), cmap=cmap, vmin=-0.6, vmax=0.6, aspect="auto")
ax2.set_xticks(range(len(colspec)))
ax2.set_xticklabels([f"config{n}" for _, n in colspec], fontsize=5, rotation=0)
ax2.set_yticks(range(len(order))); ax2.set_yticklabels(order, fontsize=5)
for k, ds in enumerate(DATASETS):
    ax2.text(2 * k + 0.5, -1.6, LABEL[ds].replace("\n", " "), ha="center", fontsize=6, fontweight="bold", color=INK2)
    if k:
        ax2.axvline(2 * k - 0.5, color="white", linewidth=2)
ax2.set_xticks(np.arange(-0.5, len(colspec), 1), minor=True); ax2.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
ax2.grid(which="minor", color="white", linewidth=0.6); ax2.tick_params(which="both", length=0)
for sp in ax2.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax2, fraction=0.025, pad=0.02)
cb.set_label("TF-MAFB expression correlation (rho)", fontsize=5); cb.ax.tick_params(labelsize=5, length=2)
ax2.set_title("TFs regulating MAFB (blank = not an input in that config)", fontsize=7, fontweight="bold", loc="left", pad=16)
fig.savefig(f"{OUT_DIR}/mafb_inputs_gaba_vs_ra.pdf", bbox_inches="tight")
fig.savefig(f"{OUT_DIR}/mafb_inputs_gaba_vs_ra.png", dpi=300, bbox_inches="tight")
print("wrote", OUT_DIR, "| TFs:", {ds: tf_tab[tf_tab.ds == ds].TF.nunique() for ds in DATASETS}, "| regions:", R.groupby("ds").size().to_dict())
