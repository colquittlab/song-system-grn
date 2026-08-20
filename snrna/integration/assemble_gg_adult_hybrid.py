"""Composite for ADULT finch (song system, HYBRID relabeling) x ADULT chicken (Zaremba et
al. 2025) -- repeat of assemble_gg_adult.py's full-suite (GSI + SAMap + CCA + SATURN)
composite, using the updated 'celltype_hybrid' finch cell-type labels (47 types) in place
of the original flat 'cluster' scheme (49 types). Same 34,295 finch cells (659 dropped for
missing celltype_hybrid), same chicken side, unchanged.

celltype_hybrid's naming already encodes song/non-song structure directly: position-
restricted song-nucleus Glut populations get their own dedicated clusters
(Glut-DACH2-HVCra, -HVCra-Pre, -HVCx, Glut-CACNA1H-RA, Glut-DACH2-LMANco), while
Glut-DACH2-LMANsh (LMAN shell, non-song per this project's prior determination) and the
widely shared non-song populations (Glut-DACH2-1..8, Glut-CACNA1H-1..4, spanning
NC/NR/Arco/HVC diffusely) are not further subdivided by position. New
developmental/ambiguous populations (Glut-Im, -NB, -NSC, -GABA, Glut-SATB2-1) have no
clear analogue in the original scheme and are left ungrouped here.

SATURN weight fixed at 0.5, NO seed-stability discount (single seed trained), matching
assemble_gg_adult.py's caveat.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from composite_score import aggregate, summarise
from class_benchmark import score as class_score

BASE = Path("/private/groups/colquittlab/saturn/zaremba_composite")
GGA = Path("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult/analysis/hybrid_seed0")
SAMAP_DIR = Path("/private/groups/colquittlab/saturn/samap_bf-adult_gg-adult/results_hybrid")
ANN = pd.read_csv("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult/data/gg_adult_label_annotation.csv",
                  index_col=0)

gsi = pd.read_csv(BASE / "results" / "gsi_corr_gg_adult_hybrid.csv", index_col=0)
sam = pd.read_csv(SAMAP_DIR / "samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")],
              [c for c in sam.columns if str(c).startswith("gg_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]
sat = pd.read_csv(GGA / "transfer_matrix.csv", index_col=0)
cca = pd.read_csv(BASE / "results" / "cca" / "gg_adult_hybrid_cca_log_finch_from_mouse_ka5_kfNA_d30_matrix.csv", index_col=0)

mats = {"gsi": gsi, "samap": sam, "saturn": sat, "cca": cca}
idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
print(f"{len(idx)} finch clusters (celltype_hybrid) x {len(cols)} chicken labels\n")

outdir_early = BASE / "results" / "composite_gg_adult_hybrid"
outdir_early.mkdir(parents=True, exist_ok=True)
for k, m in mats.items():
    m.to_csv(outdir_early / f"method_{k}_matrix.csv")
print("wrote per-method matrices (method_<name>_matrix.csv) for individual heatmaps")

weights = {"gsi": 1.0, "samap": 1.0, "cca": 1.0, "saturn": 0.5}
print(f"weights: {weights}  (SATURN static 0.5, NO seed-stability discount -- single seed only)")

rs, bo, zs = aggregate(mats, weights=weights)
outdir = BASE / "results" / "composite_gg_adult_hybrid"
D = summarise(mats, rs, bo, zs, outdir, n_top=5)

a, tab = class_score(D.top_rank_agg, ANN)
print(f"\nclass-benchmark accuracy: {a:.1%} (n={len(tab)})")

pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 34)
foc = [c for c in [
    "Astro-1", "Astro-2", "Astro-3", "OPC", "Oligo-1", "Oligo-2", "Oligo-3", "Micro", "Endo", "Epen",
    "GABA-3", "GABA-4-1", "GABA-4-2", "GABA-8", "GABA-Im",
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Pre", "Glut-DACH2-HVCx", "Glut-CACNA1H-RA",
    "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh",
    "Glut-DACH2-1", "Glut-CACNA1H-1", "Glut-Im", "Glut-NB", "Glut-NSC",
] if c in D.index]
print("\n=== focus clusters ===")
cols_show = ["top_rank_agg", "rank_agg_score", "confidence", "confidence_tier",
             "reciprocal_support_for_rank_agg", "n_clusters_sharing_winner"]
print(D.loc[foc, cols_show].to_string())
D.to_csv(outdir / "composite_calls_summary.csv")
print(f"\nwrote {outdir}")
