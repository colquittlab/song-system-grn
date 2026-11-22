"""Compare the RPCA sweep on OrthoFinder orthologs (this folder) with the original sweep on
reciprocal-best-BLASTP-hit orthologs (../rpca/results/cca), and rebuild the 3-method
composite (GSI + SAMap + RPCA at k.anchor 20) with the OrthoFinder RPCA.

Nothing in ../rpca or ../composite is read-modified; all outputs go to results/ here.

Per k.anchor:
  * reciprocal score S = sqrt(fwd * rev) (reciprocal_score.symmetrize, as in the original)
  * per finch type: best score over chicken types and the chicken type giving it
  * song vs non-song: two-sided Mann-Whitney U and Cliff's delta on the best scores
    (song_nonsong_helpers.R: 5 song vs 13 non-song types, Glut-DACH2-HVCra-Int excluded)
  * RBH vs OrthoFinder agreement: Spearman / Pearson over all finch x chicken cells and
    over per-type best scores; share of finch types with the same best chicken match

Composite: assemble_glutonly_plusSATB2_noPre.py with the RPCA input swapped; song vs
non-song on each type's row maximum of composite_confidence_matrix.csv, as
plot_song_vs_nonsong_plusSATB2_noPre.R does.

    python snrna/integration/gg_glutonly_hybrid/rpca_orthofinder/compare_rbh_orthofinder.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
GG = HERE.parent
sys.path.insert(0, str(GG.parent))  # snrna/integration -> integration_paths puts the toolkit on sys.path
import integration_paths  # noqa: F401,E402
from reciprocal_score import symmetrize  # noqa: E402
from composite_score import aggregate, reciprocal_top, summarise  # noqa: E402

KA = [2, 5, 10, 20, 30, 50, 75, 100]
OLD = ("gg_glutonly_hybrid_plusSATB2_noPre", GG / "rpca/results/cca")
NEW = ("gg_glutonly_hybrid_plusSATB2_noPre_orthofinder", HERE / "results/cca")
SONG = ["Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA"]
EXCLUDED = ["Glut-DACH2-HVCra-Int"]
OUT = HERE / "results"


def recip(tag, d, ka):
    stub = f"{tag}_rpca_SCT_%s_ka{ka}_kf200_d40_matrix.csv"
    S = symmetrize(d / (stub % "finch_from_mouse"), d / (stub % "mouse_from_finch"))
    return S.drop(index=[c for c in EXCLUDED if c in S.index])


def cliffs(x, y):
    return float(np.mean(np.sign(np.subtract.outer(np.asarray(x), np.asarray(y)))))


def song_stats(best):
    s, n = best[best.index.isin(SONG)], best[~best.index.isin(SONG)]
    u, p = mannwhitneyu(s, n, alternative="two-sided")
    return dict(n_song=len(s), n_nonsong=len(n), song_mean=s.mean(), nonsong_mean=n.mean(),
                cliffs_delta=cliffs(s, n), U=u, p_mannwhitney=p)


rows, per_type = [], []
for ka in KA:
    So, Sn = recip(*OLD, ka), recip(*NEW, ka)
    idx = So.index.intersection(Sn.index)
    cols = So.columns.union(Sn.columns)
    So, Sn = So.reindex(index=idx, columns=cols, fill_value=0), Sn.reindex(index=idx, columns=cols, fill_value=0)
    bo, bn = So.max(axis=1), Sn.max(axis=1)
    to, tn = So.idxmax(axis=1), Sn.idxmax(axis=1)
    for name, best in (("RBH", bo), ("OrthoFinder", bn)):
        rows.append(dict(k_anchor=ka, orthologs=name, **song_stats(best)))
    rows[-1].update(
        matrix_spearman=spearmanr(So.values.ravel(), Sn.values.ravel())[0],
        matrix_pearson=pearsonr(So.values.ravel(), Sn.values.ravel())[0],
        best_spearman=spearmanr(bo, bn)[0], best_pearson=pearsonr(bo, bn)[0],
        same_top_match=f"{int((to == tn).sum())}/{len(idx)}")
    per_type.append(pd.DataFrame({"k_anchor": ka, "finch_type": idx,
                                  "group": ["song" if c in SONG else "non-song" for c in idx],
                                  "best_rbh": bo.values, "top_rbh": to.values,
                                  "best_orthofinder": bn.values, "top_orthofinder": tn.values}))
    Sn.to_csv(OUT / f"reciprocal_orthofinder_ka{ka}_matrix.csv")

stats = pd.DataFrame(rows)
stats.to_csv(OUT / "song_vs_nonsong_rbh_vs_orthofinder.csv", index=False)
pd.concat(per_type).to_csv(OUT / "best_match_by_type_rbh_vs_orthofinder.csv", index=False)

# ---- composite with the OrthoFinder RPCA at k.anchor 20 -------------------------------------
gsi = pd.read_csv(GG / "gsi/results/gsi_corr_gg_glutonly_hybrid_plusSATB2_noPre_SCT.csv", index_col=0)
sam = pd.read_csv(GG / "samap/results/SCT_plusSATB2_noPre/samap_mapping_table.csv", index_col=0)
sam = sam.loc[[i for i in sam.index if str(i).startswith("bf_")], [c for c in sam.columns if str(c).startswith("gg_")]]
sam.index = [i[3:] for i in sam.index]; sam.columns = [c[3:] for c in sam.columns]
def build_composite(cca, cdir):
    """assemble_glutonly_plusSATB2_noPre.py, with the RPCA matrix passed in."""
    mats = {"gsi": gsi, "samap": sam, "cca": cca}
    mats = {k: m.drop(index=[c for c in EXCLUDED if c in m.index]) for k, m in mats.items()}
    idx = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
    cols = sorted(set.intersection(*[set(m.columns) for m in mats.values()]))
    mats = {k: m.reindex(index=idx, columns=cols).fillna(0.0) for k, m in mats.items()}
    cdir.mkdir(parents=True, exist_ok=True)
    for k, m in mats.items():
        m.to_csv(cdir / f"method_{k}_matrix.csv")
    rs, bo_, zs = aggregate(mats, weights={"gsi": 1.0, "samap": 1.0, "cca": 1.0})
    summarise(mats, rs, bo_, zs, cdir, n_top=5)
    for n in (1, 2):
        sum(reciprocal_top(m, n).astype(int) for m in mats.values()).to_csv(cdir / f"composite_agreement_count_top{n}_matrix.csv")
    return (pd.read_csv(cdir / "composite_confidence_matrix.csv", index_col=0),
            pd.read_csv(cdir / "composite_agreement_count_top2_matrix.csv", index_col=0))


# RBH composite rebuilt here (the original's top-2 agreement matrix is not tracked) and
# checked against the original confidence matrix. NB: the version committed as of this
# script predates the current run (19 rows: Glut-DACH2-HVCra-Int still scored, older SAMap
# input); the current version (18 rows) is the one behind the manuscript's numbers and was
# uncommitted in the main checkout, hence the override.
old_c, old_a = build_composite(recip(*OLD, 20), OUT / "composite_plusSATB2_noPre_rbh_rebuilt")
ref = Path(os.environ.get("COMPOSITE_REFERENCE",
                          GG / "composite/results/composite_plusSATB2_noPre/composite_confidence_matrix.csv"))
tracked = pd.read_csv(ref, index_col=0)
same_rows = sorted(tracked.index) == sorted(old_c.index)
diff = (old_c - tracked.reindex(index=old_c.index, columns=old_c.columns)).abs().max().max()
print(f"rebuilt RBH composite vs {ref}: same finch types {same_rows}, max abs diff {diff:.2e}")
if not (same_rows and diff < 1e-3):
    print("WARNING: rebuilt RBH composite does not match the reference composite")
new_c, new_a = build_composite(recip(*NEW, 20), OUT / "composite_plusSATB2_noPre_orthofinder")
crow = []
for name, C, A in (("RBH", old_c, old_a), ("OrthoFinder", new_c, new_a)):
    C = C.drop(index=[c for c in EXCLUDED if c in C.index])
    A = A.reindex(index=C.index)
    best = C.max(axis=1)
    # per finch type: methods agreeing (reciprocal top-2) on its best-agreed chicken type
    agree = A.max(axis=1)
    song = best.index.isin(SONG)
    crow.append(dict(orthologs=name, **song_stats(best),
                     all3_agree_song=f"{int((agree[song] == 3).sum())}/{song.sum()}",
                     all3_agree_nonsong=f"{int((agree[~song] == 3).sum())}/{(~song).sum()}",
                     ge2_agree_song=f"{int((agree[song] >= 2).sum())}/{song.sum()}",
                     ge2_agree_nonsong=f"{int((agree[~song] >= 2).sum())}/{(~song).sum()}"))
cstats = pd.DataFrame(crow)
cstats.to_csv(OUT / "composite_song_vs_nonsong_rbh_vs_orthofinder.csv", index=False)

pd.set_option("display.width", 250)
print(stats.round(4).to_string(index=False))
print()
print(cstats.round(4).to_string(index=False))
