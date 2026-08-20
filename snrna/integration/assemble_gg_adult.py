"""Composite for ADULT finch (song system) x ADULT chicken (Zaremba et al. 2025).

Both datasets are adult, so this is the closest match in the whole project to
Zaremba's own protocol: all four methods included (GSI + SAMap + CCA + SATURN), same
as the finch-mature x Yao comparison, unlike the two dev-vs-dev pairs which drop GSI
per Zaremba's explicit rule.

SATURN weight is fixed at 0.5 (only method consuming the query labels; worst-agreeing
with label-free evidence throughout this project) WITHOUT a per-cluster seed-stability
discount, because only a single seed was trained for this pair -- no replicate seeds
exist to estimate stability from. This is a real gap relative to the Yao/La Manno
composites and should be treated as a caveat on any single-cluster SATURN-driven claim
here specifically.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from composite_score import aggregate, summarise
from class_benchmark import score as class_score

BASE = Path("/private/groups/colquittlab/saturn/zaremba_composite")
GGA = Path("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult/analysis/macro2000_hv8000_seed0")
SAMAP_DIR = Path("/private/groups/colquittlab/saturn/samap_bf-adult_gg-adult/results")
ANN = pd.read_csv("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult/data/gg_adult_label_annotation.csv",
                  index_col=0)

gsi = pd.read_csv(BASE / "results" / "gsi_corr_gg_adult.csv", index_col=0)
sam = pd.read_csv(SAMAP_DIR / "samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")],
              [c for c in sam.columns if str(c).startswith("gg_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]
sat = pd.read_csv(GGA / "transfer_matrix.csv", index_col=0)
cca = pd.read_csv(BASE / "results" / "cca" / "gg_adult_cca_matrix_finch_from_mouse.csv", index_col=0)

mats = {"gsi": gsi, "samap": sam, "saturn": sat, "cca": cca}
idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
print(f"{len(idx)} finch clusters x {len(cols)} chicken labels\n")

weights = {"gsi": 1.0, "samap": 1.0, "cca": 1.0, "saturn": 0.5}
print(f"weights: {weights}  (SATURN static 0.5, NO seed-stability discount -- single seed only)")

rs, bo, zs = aggregate(mats, weights=weights)
outdir = BASE / "results" / "composite_gg_adult"
D = summarise(mats, rs, bo, zs, outdir, n_top=5)

a, tab = class_score(D.top_rank_agg, ANN)
print(f"\nclass-benchmark accuracy: {a:.1%} (n={len(tab)})")

pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 34)
foc = [c for c in ["Astro-1","Astro-2","Astro-3","OPC","Oligo-1","Oligo-2","Oligo-3","Micro","Endo","Epen",
                   "GABA-3","GABA-4-1","GABA-4-2","GABA-8","GABA-Pre","Glut-HVC-1","Glut-RA","Glut-LMAN-1",
                   "Glut-Arco-1","Glut-Arco-2","Glut-Nido-1"] if c in D.index]
print("\n=== focus clusters ===")
cols_show = ["top_rank_agg", "rank_agg_score", "confidence", "confidence_tier",
             "reciprocal_support_for_rank_agg", "n_clusters_sharing_winner"]
print(D.loc[foc, cols_show].to_string())
D.to_csv(outdir / "composite_calls_summary.csv")
print(f"\nwrote {outdir}")
