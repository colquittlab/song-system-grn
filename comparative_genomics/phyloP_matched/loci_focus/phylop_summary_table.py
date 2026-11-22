#!/usr/bin/env python
"""Per-class summary of constraint at high-accessibility-difference bases, for plotting.

From the per-element tables written by phylop_density.py (unsmoothed runs), for each span
(domain, +/-1.5 Mb, genome-wide), track, window (2 kb / 500 bp / 200 bp) and element set
(all; Log2FC-matched) and for each class separately -- Glut-RA activity bases (Glut-RA minus
Glut-Arco-1 fragment difference > 5) and GABA-4 activity bases (GABA-4-1 minus GABA-1-1) --
computes, pooled over that class's high-difference bases:
  pct   percentage with phyloP > Q   (Q < 0: percentage with phyloP < Q -- acceleration)
  mean  mean phyloP of those bases
each with a 95% peak-cluster bootstrap interval (elements resampled with all their bases).
The G4 - RA difference and its p come from the same run's summary table.

Writes results/compare_mafb_spans/phylop_summary_q<Q>.tsv; plot_phylop_summary.R draws it.
"""
import os, sys
import numpy as np, pandas as pd

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import result

Q = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
NBOOT = 2000
SPANS = {"mafb_domain": "DHX35–TOP1 domain", "mafb_1p5mb": "MAFB ±1.5 Mb", "genome_wide": "genome-wide"}
PAIRS = {"RA": "Glut-RA − Glut-Arco-1", "G4": "GABA-4-1 − GABA-1-1"}
TRACKS = {"": ["avian all", "vert. all", "oscine (Anc218)"],
          "_perbase": ["whole tree (matched)", "oscine lineage, birdAnc222 (matched)"]}
WINDOWS = {"": "2 kb", "_w500": "500 bp", "_w200": "200 bp"}
SETS = {"": "all", "_matchlfc": "Log2FC-matched"}
CENTERS = {"": "element centre", "_snapcenter": "SnapATAC midpoint", "_summit": "SnapATAC summit"}


def main():
    q = f"{abs(Q):g}"; rows = []; rng = np.random.default_rng(5)
    op = "<-" if Q < 0 else ">"                       # column operator in the per-element tables
    for tr, tracks in TRACKS.items():
        for st, stlab in SETS.items():
            for w, wlab in WINDOWS.items():
              for ce, celab in CENTERS.items():
                for sp, splab in SPANS.items():
                    stem = f"{tr}{st}{w}{ce}_nosmooth"
                    pk, sm = (f"{F}/results/{sp}/phylop_density_{k}{stem}.tsv" for k in ("peaks", "summary"))
                    if not (os.path.exists(pk) and os.path.exists(sm)): continue
                    d = pd.read_csv(pk, sep="\t"); S = pd.read_csv(sm, sep="\t").set_index("track")
                    for t in tracks:
                        if f"tn{op}{q}|{PAIRS['RA']}|{t}" not in d: continue
                        for cl, pr in PAIRS.items():
                            n = d[f"nbins|{pr}"].values; tn = d[f"tn{op}{q}|{pr}|{t}"].values; ts = d[f"ts{op}{q}|{pr}|{t}"].values
                            keep = n > 0; n, tn, ts = n[keep], tn[keep], ts[keep]
                            pct = 100 * tn.sum() / n.sum() if n.sum() else np.nan
                            mean = ts.sum() / tn.sum() if tn.sum() else np.nan
                            bp, bm = [], []
                            for _ in range(NBOOT):
                                ix = rng.integers(0, len(n), len(n))
                                if n[ix].sum(): bp.append(100 * tn[ix].sum() / n[ix].sum())
                                if tn[ix].sum(): bm.append(ts[ix].sum() / tn[ix].sum())
                            s = S.loc[t]
                            rows.append(dict(tracks="per-base matched" if tr else "3-model", track=t, set=stlab, window=wlab, centring=celab,
                                             span=splab, cls=cl, n_elements=int(keep.sum()), n_bases=int(n.sum()),
                                             pct=pct, pct_lo=np.percentile(bp, 2.5) if bp else np.nan,
                                             pct_hi=np.percentile(bp, 97.5) if bp else np.nan,
                                             mean=mean, mean_lo=np.percentile(bm, 2.5) if bm else np.nan,
                                             mean_hi=np.percentile(bm, 97.5) if bm else np.nan,
                                             pct_diff_p=s.get(f"base_frac>{q}_p" if Q > 0 else f"base_acc<-{q}_p", np.nan),
                                             mean_diff_p=s.get(f"tailmean>{q}_p" if Q > 0 else f"accmean<-{q}_p", np.nan)))
    D = pd.DataFrame(rows)
    out = result(f"phylop_summary_q{Q:g}.tsv", "compare_mafb_spans")
    D.to_csv(out, sep="\t", index=False, float_format="%.4g")
    print(f"  wrote {out}: {len(D)} rows")


if __name__ == "__main__":
    main()
