#!/usr/bin/env python
"""Reverse-ascertainment calibration on the corrected presence call.

gain_calib_radius.py / gain_calib_sweep.py rerun with presence.project() in
place of row concatenation. Everything else is unchanged: the domain is tiled
in 76 bp windows at a 38 bp step on lonStrDom2 coordinates; each of fourteen
genomes in turn plays the role lonStrDom2 played, and a window counts for that
genome when its own sequence carries a site; near/far sets are matched on
patristic distance (radius 0.11, >= 8 near and >= 60 far); gain = fraction of
near genomes with a site minus fraction of far genomes.

Coverage here is per window: a genome counts only if >= 80% of the window's
reference positions carry a base, so a genome with no orthologous sequence is
missing rather than absent -- in the ascertainment genome itself, and in both
comparison arms.

All four calling thresholds are scored in one pass. Writes data/gain_sweep.json
(the schema the report reads; the original kept as gain_sweep_concat.json) and
data/gain_calib_radius.json for the 99.99th-percentile threshold.
"""
import collections, json, os, shutil, sys
import numpy as np

sys.path.insert(0, ".")
import presence as PR
from gain_calib_radius import (LO, dist, RADIUS, MINNEAR, ASC_OSC, ASC_OUT, OSC, OUT,
                               WIN, STEP)

QS = ["99.99", "99.9", "99.5", "99.0"]
RC = {k: v[::-1, ::-1] for k, v in LO.items()}
CODE = np.full(256, 4, np.int8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i; CODE[ord(c.lower())] = i


def thresh(lo, q, n=200000):
    rng = np.random.default_rng(0); w = len(lo)
    s = rng.integers(0, 4, (n, w))
    return float(np.percentile(np.array([lo[k][s[:, k]] for k in range(w)]).sum(0), q))


TH = {(mo, q): thresh(lo, float(q)) for mo, lo in LO.items() for q in QS}


def best(seq, mo):
    lo = LO[mo]; w = len(lo)
    a = CODE[np.frombuffer(seq.encode(), np.uint8)]
    if a.size < w: return None
    sl = np.lib.stride_tricks.sliding_window_view(a, w)
    ok = (sl < 4).all(1)
    if not ok.any(): return None
    sl = sl[ok]; idx = np.arange(w)
    return float(max(lo[idx, sl].sum(1).max(), RC[mo][idx, sl].sum(1).max()))


def main():
    for f in ("data/gain_sweep.json", "data/gain_calib_radius.json"):
        b = f.replace(".json", "_concat.json")
        if os.path.exists(f) and not os.path.exists(b): shutil.copy(f, b)
    B = PR.load("data/domain_gap.maf")
    blocks, starts = B
    ext = sorted(set(OSC) | set(OUT))
    wins = set()
    for bs, be, _, _ in blocks:
        for w in range((bs // STEP) * STEP, be, STEP): wins.add(w)
    wins = sorted(wins)
    print(f"  blocks {len(blocks):,}   windows {len(wins):,}", flush=True)

    # one projection per window; scores per genome per matrix
    score = {}
    for n, w in enumerate(wins):
        bl = PR.overlapping(blocks, starts, w, w + WIN)
        if not bl: continue
        pr = PR.project(bl, w, w + WIN, w, w + WIN, genomes=ext)
        cov = {g: v for g, v in pr.items() if PR.covered(v[1], w, w + WIN)}
        if len(cov) < 200: continue
        score[w] = {mo: {g: best(v[0], mo) for g, v in cov.items()} for mo in LO}
        if (n + 1) % 400 == 0: print(f"    {n + 1}/{len(wins)} windows", flush=True)
    print(f"  windows with >=200 covered genomes: {len(score):,}", flush=True)

    sweep = []; radius = {}
    for q in QS:
        for mo in LO:
            th = TH[(mo, q)]
            per = {}
            for G in ASC_OSC + ASC_OUT:
                gains = []
                for w, sc in score.items():
                    d = sc[mo]
                    if G not in d or d[G] is None or d[G] < th: continue
                    dd = [(dist(G, g), g) for g in d if g != G]
                    dd = [(x, g) for x, g in dd if x is not None]
                    if len(dd) < 80: continue
                    near = [g for x, g in dd if x <= RADIUS]
                    far = [g for x, g in dd if x > RADIUS]
                    if len(near) < MINNEAR or len(far) < 60: continue
                    pres = lambda g: d[g] is not None and d[g] >= th
                    gains.append(sum(map(pres, near)) / len(near) - sum(map(pres, far)) / len(far))
                if len(gains) >= 8: per[G] = gains
            o = [np.mean(v) for G, v in per.items() if G in OSC]
            u = [np.mean(v) for G, v in per.items() if G in OUT]
            if not o or not u: continue
            row = dict(q=q, motif=mo, osc=round(float(np.mean(o)), 3),
                       out=round(float(np.mean(u)), 3),
                       diff=round(float(np.mean(o) - np.mean(u)), 3),
                       win=int(round(np.mean([len(v) for v in per.values()]))),
                       n_osc=len(o), n_out=len(u))
            sweep.append(row)
            if q == "99.99":
                radius.update({f"{mo}|{G}": float(np.mean(v)) for G, v in per.items()})
            print(f"  {q:>6}th  {mo:<7} oscine-ascertained {row['osc']:+.3f} ({len(o)} genomes)"
                  f"   outgroup-ascertained {row['out']:+.3f} ({len(u)})"
                  f"   difference {row['diff']:+.3f}   ~{row['win']} windows/genome", flush=True)
    json.dump(sweep, open("data/gain_sweep.json", "w"), indent=1)
    json.dump(radius, open("data/gain_calib_radius.json", "w"), indent=1)
    print("  wrote data/gain_sweep.json and data/gain_calib_radius.json")


if __name__ == "__main__":
    main()
