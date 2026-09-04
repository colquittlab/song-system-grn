"""3-method composite (GSI + SAMap + RPCA/CCA) for the plusSATB2_noPre finch/chicken
subset -- Glut-SATB2-1 added back to the finch side, chicken Ex_Pre_KCNH7/Ex_Pre_SATB2
removed from the reference -- the current/most up-to-date test set in this analysis.

SATURN is DROPPED entirely for this composite (per explicit request, after confirming no
SATURN embedding has ever been trained on this exact subset -- only on the original base
glutonly-hybrid dataset; reusing that mismatched embedding as a proxy was explicitly
declined in favour of a clean 3-method comparison). Matches assemble_lamanno.py's
established precedent of falling back to a 3-method composite when one method's input
isn't available for a given comparison, rather than reusing a mismatched substitute.

"CCA" here is this project's RPCA reciprocity-weighted score (reciprocal_score.symmetrize
-- sqrt(fwd*rev), the fix for TransferData's one-directional vote artifact), NOT Seurat's
separate CCA-reduction integration method -- "RPCA" and "CCA" have been used
interchangeably for this method slot throughout this analysis (matching
assemble_glutonly.py's own mats dict key "cca" for what is actually an RPCA-reduction
run); per explicit request, ka=20 (not the ka=50 used as the headline example
elsewhere in this analysis) is the fixed k.anchor for this composite's CCA input.
k.filter=200 is nominal only under SCT (Seurat forces it to NA internally).

Orthologs for GSI are OrthoFinder (matching assemble_glutonly.py's convention); SAMap
uses its own BLAST-graph homology, independent of any ortholog table.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, "/private/groups/colquittlab/song-system-grn/snrna/integration")
sys.path.insert(0, "/private/groups/colquittlab/finch-integration-toolkit")
from composite_score import aggregate, summarise
from reciprocal_score import symmetrize
from position_palette import EXCLUDED_CLUSTERS, song_group

GG_HYBRID = Path("/private/groups/colquittlab/song-system-grn/snrna/integration/gg_glutonly_hybrid")
GSI_DIR = GG_HYBRID / "gsi" / "results"
RPCA_DIR = GG_HYBRID / "rpca" / "results" / "cca"
SAMAP_DIR = GG_HYBRID / "samap" / "results" / "plusSATB2_noPre"
COMPOSITE_OUT = GG_HYBRID / "composite" / "results"

TAG = "gg_glutonly_hybrid_plusSATB2_noPre"
KA, KF, DIMS = 20, 200, 40

gsi = pd.read_csv(GSI_DIR / "gsi_corr_gg_glutonly_hybrid_plusSATB2_noPre_SCT.csv",
                  index_col=0)

sam = pd.read_csv(SAMAP_DIR / "samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")],
              [c for c in sam.columns if str(c).startswith("gg_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]

cca = symmetrize(
    RPCA_DIR / f"{TAG}_rpca_SCT_finch_from_mouse_ka{KA}_kf{KF}_d{DIMS}_matrix.csv",
    RPCA_DIR / f"{TAG}_rpca_SCT_mouse_from_finch_ka{KA}_kf{KF}_d{DIMS}_matrix.csv",
)

mats = {"gsi": gsi, "samap": sam, "cca": cca}
idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
print(f"{len(idx)} finch Glut clusters x {len(cols)} chicken Excitatory labels (3-method: {list(mats)})\n")

outdir = COMPOSITE_OUT / "composite_plusSATB2_noPre"
outdir.mkdir(parents=True, exist_ok=True)
for k, m in mats.items():
    m.to_csv(outdir / f"method_{k}_matrix.csv")
print(f"wrote per-method matrices (method_<name>_matrix.csv) for individual heatmaps")

weights = {"gsi": 1.0, "samap": 1.0, "cca": 1.0}
print(f"weights: {weights}  (no SATURN, no seed-stability discount needed)")

rs, bo, zs = aggregate(mats, weights=weights)
D = summarise(mats, rs, bo, zs, outdir, n_top=5)

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

# Per-method top call, side by side -- the direct answer to "how consistent are the top
# calls across methods" (each method's own argmax chicken cluster per finch cluster).
print("\n=== per-method top call (argmax), side by side ===")
per_method_top = pd.DataFrame({k: m.idxmax(axis=1) for k, m in mats.items()})
per_method_top["n_methods_agreeing"] = per_method_top.apply(lambda r: r.value_counts().iloc[0], axis=1)
per_method_top.to_csv(outdir / "per_method_top_calls.csv")
print(per_method_top.to_string())
print(f"\nmethods fully agreeing (all 3 same top call): {(per_method_top['n_methods_agreeing'] == 3).sum()}/{len(per_method_top)}")
print(f"exactly 2/3 agreeing: {(per_method_top['n_methods_agreeing'] == 2).sum()}/{len(per_method_top)}")
print(f"all 3 disagreeing: {(per_method_top['n_methods_agreeing'] == 1).sum()}/{len(per_method_top)}")

from scipy.stats import mannwhitneyu
song = D.loc[D.group == "song", "confidence"]
nonsong = D.loc[D.group == "non-song", "confidence"]
u, p = mannwhitneyu(song, nonsong, alternative="less")
print(f"\nMann-Whitney U (song < non-song), one-sided: U={u:.1f}  p={p:.4f}")

print(f"\nwrote {outdir}")
