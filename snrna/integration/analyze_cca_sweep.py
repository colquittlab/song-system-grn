"""Weak->strong k.anchor titration: does confidence rise earlier for non-song than song?

For each k.anchor, load both CCA direction matrices (finch<-mouse gives, per finch
cluster, frac*mean_score against every chicken label), take each finch cluster's own
best (top-1) score as its "confidence at this integration strength", and track that
across the sweep. The prediction: non-song clusters should reach a given confidence
threshold at LOWER k.anchor (weaker integration) than song clusters, i.e. their score
curves should sit above the song curves throughout, not just at one endpoint.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from assemble_glutonly import song_group  # reuse the exact, already-corrected grouping

import sys as _sys
R = Path("/private/groups/colquittlab/saturn/zaremba_composite/results/cca")
KA = [2, 5, 10, 20, 30, 50]
REDUCTION = _sys.argv[1] if len(_sys.argv) > 1 else "cca"

rows = []
for ka in KA:
    f = R / f"gg_glutonly_{REDUCTION}_matrix_finch_from_mouse_ka{ka}.csv"
    if not f.exists():
        print(f"missing: {f}")
        continue
    M = pd.read_csv(f, index_col=0)
    top = M.max(axis=1)
    for c, v in top.items():
        rows.append(dict(k_anchor=ka, cluster=c, group=song_group(c), best_score=v))

D = pd.DataFrame(rows)
if D.empty:
    print("no data yet")
    sys.exit(0)

pd.set_option("display.width", 200)
piv = D.pivot(index="cluster", columns="k_anchor", values="best_score")
piv["group"] = [song_group(c) for c in piv.index]
piv = piv.sort_values("group")
print("=== per-cluster best CCA score across k.anchor ===")
print(piv.round(3).to_string())

print("\n=== group mean best-score by k.anchor (the titration curve) ===")
curve = D.groupby(["k_anchor", "group"])["best_score"].mean().unstack()
print(curve.round(3).to_string())

print("\n=== group median best-score by k.anchor ===")
curve_med = D.groupby(["k_anchor", "group"])["best_score"].median().unstack()
print(curve_med.round(3).to_string())

# Threshold-crossing: at what k.anchor does each cluster first exceed 0.5?
THRESH = 0.5
cross = {}
for c in piv.index:
    row = piv.loc[c, KA]
    above = row[row >= THRESH]
    cross[c] = above.index.min() if len(above) else np.nan
X = pd.DataFrame({"cluster": list(cross.keys()), "group": [song_group(c) for c in cross],
                  "k_anchor_first_above_0.5": list(cross.values())}).set_index("cluster")
print(f"\n=== first k.anchor at which best-score >= {THRESH} ===")
print(X.sort_values(["group", "k_anchor_first_above_0.5"]).to_string())
print("\nmean crossing point by group (NaN = never crossed, excluded from mean):")
print(X.groupby("group")["k_anchor_first_above_0.5"].agg(["mean", "count", lambda s: s.isna().sum()]))

out = Path("/private/groups/colquittlab/saturn/zaremba_composite/results/composite_gg_glutonly")
piv.to_csv(out / f"{REDUCTION}_kanchor_sweep_per_cluster.csv")
curve.to_csv(out / f"{REDUCTION}_kanchor_sweep_group_mean.csv")
print(f"\nwrote {out}/{REDUCTION}_kanchor_sweep_*.csv")
