"""Composite for the Glut(finch)-only x Excitatory(chicken)-only adult subset.

Manuscript context: an RPCA analysis found that song-nucleus Glut types (HVC, RA, LMAN,
Area X/Arco) do not map onto any single well-defined chicken excitatory type, while
non-song Glut types (NC, NR) do -- interpreted as the song nuclei containing derived,
novel cell types rather than conserved ones. This assembles the same 4-method ensemble
(GSI + SAMap + CCA + SATURN) used throughout this project, restricted to the identical
Glut-vs-Excitatory subset, as an independent replication for reviewers.

Orthologs for GSI/CCA are OrthoFinder 2.5.4 (not RBH) for this pair specifically, per
explicit request -- matches Zaremba's own cited method with no substitution.

SATURN weight fixed at 0.5, NO seed-stability discount (single seed trained; same caveat
as the full adult composite).

Song vs non-song grouping is read directly off the cluster name (HVC-/RA/LMAN- =
song nucleus; NC-/NR-/Arco- = non-song), matching the manuscript's own framing -- Arco
(arcopallium generically) is distinct from Area X/song-specific striatal nuclei here.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from composite_score import aggregate, summarise
from class_benchmark import score as class_score

BASE = Path("/private/groups/colquittlab/saturn/zaremba_composite")
GGA = Path("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult-glutonly/analysis/macro2000_hv8000_seed0")
SAMAP_DIR = Path("/private/groups/colquittlab/saturn/samap_bf-adult_gg-adult-glutonly/results")
ANN = pd.read_csv("/private/groups/colquittlab/saturn/snrna-bf-adult_snrna-gg-adult-glutonly/data/gg_glutonly_label_annotation.csv",
                  index_col=0)

gsi = pd.read_csv(BASE / "results" / "gsi_corr_gg_glutonly.csv", index_col=0)
sam = pd.read_csv(SAMAP_DIR / "samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")],
              [c for c in sam.columns if str(c).startswith("gg_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]
sat = pd.read_csv(GGA / "transfer_matrix.csv", index_col=0)
cca = pd.read_csv(BASE / "results" / "cca" / "gg_glutonly_cca_matrix_finch_from_mouse.csv", index_col=0)

mats = {"gsi": gsi, "samap": sam, "saturn": sat, "cca": cca}
idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
print(f"{len(idx)} finch Glut clusters x {len(cols)} chicken Excitatory labels\n")

outdir_early = BASE / "results" / "composite_gg_glutonly"
outdir_early.mkdir(parents=True, exist_ok=True)
for k, m in mats.items():
    m.to_csv(outdir_early / f"method_{k}_matrix.csv")
print(f"wrote per-method matrices (method_<name>_matrix.csv) for individual heatmaps")

weights = {"gsi": 1.0, "samap": 1.0, "cca": 1.0, "saturn": 0.5}
print(f"weights: {weights}  (SATURN static 0.5, NO seed-stability discount -- single seed only)")

rs, bo, zs = aggregate(mats, weights=weights)
outdir = BASE / "results" / "composite_gg_glutonly"
D = summarise(mats, rs, bo, zs, outdir, n_top=5)

# Excluded from the song-vs-non-song comparison entirely (not relabeled as either group).
EXCLUDED_CLUSTERS = {"Glut-HVC-1a"}

def song_group(c):
    # Arco is NOT a song nucleus in this dissection scheme -- only HVC/RA/LMAN are.
    # Glut-LMAN-1 is reclassified non-song: complementary Xenium data places it in LMAN
    # shell (not LMAN core), which is not part of the classic song circuit. Glut-LMAN-2
    # is retained as song (LMAN core).
    if c == "Glut-LMAN-1":
        return "non-song"
    return "song" if any(c.startswith(p) for p in ["Glut-HVC", "Glut-RA", "Glut-LMAN"]) else "non-song"
D["group"] = [song_group(c) for c in D.index]
D.to_csv(outdir / "composite_calls_summary.csv")
D = D.drop(index=[c for c in EXCLUDED_CLUSTERS if c in D.index])

pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 34)
cols_show = ["group", "top_rank_agg", "rank_agg_score", "confidence", "confidence_tier",
             "reciprocal_support_for_rank_agg", "n_clusters_sharing_winner"]
print("\n=== ALL clusters, sorted by confidence within group ===")
print(D.sort_values(["group", "confidence"], ascending=[True, False])[cols_show].to_string())

print("\n=== SONG vs NON-SONG: confidence summary ===")
g = D.groupby("group")["confidence"].agg(["count", "mean", "median", "std", "min", "max"])
print(g.round(3).to_string())

from scipy.stats import mannwhitneyu
song = D.loc[D.group == "song", "confidence"]
nonsong = D.loc[D.group == "non-song", "confidence"]
u, p = mannwhitneyu(song, nonsong, alternative="less")
print(f"\nMann-Whitney U (song < non-song), one-sided: U={u:.1f}  p={p:.4f}")

print(f"\nwrote {outdir}")
