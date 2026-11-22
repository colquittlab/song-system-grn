#!/usr/bin/env python
"""Collate the RA vs G4 phyloP comparison with and without an alignment-depth filter.

Reads phylop_density_summary_perbase[_depth<T>][_matchlfc]_w<W>_summit_nosmooth.tsv (phylop_density.py
--tracks perbase --center summit --smooth 1 [--depth T] [--match lfc] --window W) for each span and
lists, per track, G4 minus RA with 95% peak-cluster bootstrap intervals and p: pooled mean phyloP of
high-difference bases, share > 1 and > 2, and the elements kept.

    python depth_filter_compare.py
Writes results/compare_mafb_spans/depth_filter_compare.tsv.
"""
import os
import pandas as pd

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results"
SPANS = {"domain": "mafb_domain", "1p5mb": "mafb_1p5mb", "gw": "genome_wide"}
TR = {"whole tree (matched)": "whole", "oscine lineage, birdAnc222 (matched)": "oscine"}


def main():
    rows = []
    for sp, d in SPANS.items():
        for W in (500, 200):
            for m in ("", "_matchlfc"):
                for dep in ("", "_depth0.8", "_depth0.9"):
                    f = f"{R}/{d}/phylop_density_summary_perbase{dep}{m}_w{W}_summit_nosmooth.tsv"
                    if not os.path.exists(f): continue
                    for r in pd.read_csv(f, sep="\t").itertuples(index=False):
                        g = r._asdict() if hasattr(r, "_asdict") else dict(r)
                        x = pd.read_csv(f, sep="\t"); x = x[x.track == r.track].iloc[0]
                        rows.append(dict(span=sp, window=W, matched="lfc" if m else "", depth=dep.replace("_depth", "") or "none",
                                         track=TR.get(r.track, r.track), n_RA=x.n_peaks_RA, n_G4=x.n_peaks_G4,
                                         **{f"{k}": f"{100 * x[f'{c}_diff'] if 'frac' in c else x[f'{c}_diff']:+.2f} "
                                                     f"[{100 * x[f'{c}_lo'] if 'frac' in c else x[f'{c}_lo']:+.2f}, "
                                                     f"{100 * x[f'{c}_hi'] if 'frac' in c else x[f'{c}_hi']:+.2f}]"
                                            for k, c in (("mean G4-RA", "base_mean"), ("%>1 G4-RA", "base_frac>1"), ("%>2 G4-RA", "base_frac>2"))},
                                         p_mean=x.base_mean_p, p_frac1=x["base_frac>1_p"], p_frac2=x["base_frac>2_p"]))
    T = pd.DataFrame(rows).drop_duplicates()
    T.to_csv(f"{R}/compare_mafb_spans/depth_filter_compare.tsv", sep="\t", index=False, float_format="%.3g")
    pd.set_option("display.width", 250)
    for (W, m), t in T.groupby(["window", "matched"], sort=False):
        print(f"\n== window {W} bp{', Log2FC-matched' if m else ''}")
        print(t.drop(columns=["window", "matched"]).to_string(index=False, float_format=lambda v: f"{v:.2g}"))


if __name__ == "__main__":
    main()
