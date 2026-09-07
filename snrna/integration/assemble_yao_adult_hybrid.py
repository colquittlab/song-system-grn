"""Composite for ADULT finch (song system, HYBRID celltype_hybrid labels) x Yao et al.
2023 mouse telencephalon (ABC Atlas) -- full-suite (GSI + SAMap + CCA + SATURN), same
methodology as assemble_gg_adult_hybrid.py (chicken pairing), swapping in the mouse Yao
reference. Finch side is the identical bf_adult_hybrid.h5ad object used throughout the
hybrid-label analyses.

SATURN weight fixed at 0.5, NO seed-stability discount (single seed trained).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from composite_score import aggregate, summarise
from class_benchmark import score as class_score

COMPOSITE = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/composite_scoring")
RPCA = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/rpca_sweep")
YAOA = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/snrna-bf-adult_snrna-yao/analysis/macro2000_hv8000_seed0_ep50")
SAMAP_DIR = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/datasets/samap_bf-mature_yao/results_adult_hybrid")
ANN = pd.read_csv("/private/groups/colquittlab/saturn/snrna-bf-dev_snrna-yao2023/data/yao_label_annotation.csv",
                  index_col=0)

gsi = pd.read_csv(COMPOSITE / "results" / "gsi_corr_yao_adult_hybrid.csv", index_col=0)
sam = pd.read_csv(SAMAP_DIR / "samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")],
              [c for c in sam.columns if str(c).startswith("mm_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]
sat = pd.read_csv(YAOA / "transfer_matrix.csv", index_col=0)
cca = pd.read_csv(RPCA / "results" / "cca" / "yao_adult_hybrid_cca_log_finch_from_mouse_ka5_kfNA_d30_matrix.csv", index_col=0)

mats = {"gsi": gsi, "samap": sam, "saturn": sat, "cca": cca}
idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
print(f"{len(idx)} finch clusters (celltype_hybrid) x {len(cols)} Yao mouse labels\n")

outdir_early = COMPOSITE / "results" / "yao_adult_hybrid"
outdir_early.mkdir(parents=True, exist_ok=True)
for k, m in mats.items():
    m.to_csv(outdir_early / f"method_{k}_matrix.csv")
print("wrote per-method matrices (method_<name>_matrix.csv) for individual heatmaps")

weights = {"gsi": 1.0, "samap": 1.0, "cca": 1.0, "saturn": 0.5}
print(f"weights: {weights}  (SATURN static 0.5, NO seed-stability discount -- single seed only)")

rs, bo, zs = aggregate(mats, weights=weights)
outdir = COMPOSITE / "results" / "yao_adult_hybrid"
D = summarise(mats, rs, bo, zs, outdir, n_top=5)

a, tab = class_score(D.top_rank_agg, ANN)
print(f"\nclass-benchmark accuracy: {a:.1%} (n={len(tab)})")

pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 34)
foc = [c for c in [
    "Astro-1", "Astro-2", "Astro-3", "OPC", "Oligo-1", "Oligo-2", "Oligo-3", "Micro", "Endo", "Epen",
    "GABA-3", "GABA-4-1", "GABA-4-2", "GABA-8", "GABA-Im",
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCra-Int", "Glut-DACH2-HVCx", "Glut-CACNA1H-RA",
    "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh",
    "Glut-DACH2-1", "Glut-CACNA1H-1", "Glut-Im", "Glut-NB", "Glut-NSC",
] if c in D.index]
print("\n=== focus clusters ===")
cols_show = ["top_rank_agg", "rank_agg_score", "confidence", "confidence_tier",
             "reciprocal_support_for_rank_agg", "n_clusters_sharing_winner"]
print(D.loc[foc, cols_show].to_string())
D.to_csv(outdir / "composite_calls_summary.csv")
print(f"\nwrote {outdir}")
