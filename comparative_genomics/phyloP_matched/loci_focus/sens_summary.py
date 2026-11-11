#!/usr/bin/env python
"""Summarise the motif-definition sensitivity analysis and draw it.

For each regulator and statistic: the medoid's value, the spread across the
regulator's defensible matrices, sign agreement with the medoid, the fraction
reaching p < 0.05, and a call:

  stable null      medoid p >= 0.05, and <= 20% of matrices reach p < 0.05
  robust effect    medoid p < 0.05, >= 80% share its sign, >= 50% reach p < 0.05
  matrix-dependent anything else -- the conclusion changes with the matrix
  untestable       fewer than 3 matrices give a value

The call uses SCENIC+ cluster members and HOMER known motifs. HOMER de novo
motifs were discovered on Glut-RA DARs, which biases a Glut-RA vs GABA-4 contrast;
they are drawn, and summarised separately, but not used for the call.

AR under every definition: for each AR matrix, AR's rank among the regulators
(each other regulator at its medoid), per statistic.

Writes results/regulator_motifs/{sens_summary.json, fig_motif_sensitivity.png/pdf}.
"""
import collections, json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"], "pdf.fonttype": 42})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

F = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0, F)
from span import result

PH = json.load(open(result("sens_phylop.json", "regulator_motifs")))
GA = {r["id"]: r for r in json.load(open(result("sens_gain.json", "regulator_motifs")))}
CHOICE = json.load(open(result("tf_scenic_choice.json", "regulator_motifs")))
ORDER = ["AR", "KLF12", "NFIA", "LIN28B", "STAT1", "FOSL2", "MAF", "SOX6", "TCF12", "ARX", "SATB1", "MEIS2"]
LABEL = {"KLF12": "KLF", "AR": "AR (= NR3C2)"}
REGUL = {t: CHOICE[t]["regulon"] for t in ORDER}
INK, INK3, RULE = "#171B22", "#6E7887", "#DDE1E7"
RA, G4 = "#C8161A", "#2C6E7F"

# statistic: (label, getter -> (value, p) or None)
def ph(key):
    def g(r):
        v = r.get(key)
        return (v["dod"], v["p"]) if v else None
    return g
STATS = [
    ("oscine DoD, genome-wide", ph("subtree222|gw")),
    ("whole-tree DoD, genome-wide", ph("wholetree|gw")),
    ("outlier z, genome-wide", lambda r: (r["outlier_gw"]["z"], r["outlier_gw"]["p"]) if r.get("outlier_gw") else None),
    ("oscine DoD, DHX35–TOP1 domain", ph("subtree222|domain")),
    ("oscine DoD, MAFB ±1.5 Mb", ph("subtree222|MAFB")),
    ("mean site gain, domain", lambda r: ((GA[r["id"]]["domain"]["mean_gain"], None)
                                          if GA.get(r["id"]) and GA[r["id"]].get("domain") else None)),
]


def call(med, vals):
    """vals: list of (value, p) over cluster + HOMER known matrices."""
    if len(vals) < 3 or med is None: return "untestable"
    v = np.array([x[0] for x in vals]); sgn = np.mean(np.sign(v) == np.sign(med[0]))
    if med[1] is None:      # no p (site gain): sign agreement only
        return "sign-stable" if sgn >= .8 else "matrix-dependent"
    ps = np.array([x[1] for x in vals]); fs = np.mean(ps < .05)
    if med[1] >= .05 and fs <= .2: return "stable null"
    if med[1] < .05 and sgn >= .8 and fs >= .5: return "robust effect"
    return "matrix-dependent"


def main():
    by = collections.defaultdict(list)
    for r in PH: by[r["tf"]].append(r)
    out = {}; lines = []
    for lab, get in STATS:
        out[lab] = {}
        for tf in ORDER:
            R = by[tf]
            med = next((get(r) for r in R if r["medoid"]), None)
            main = [get(r) for r in R if r["source"] != "denovo" and get(r)]
            dn = [get(r) for r in R if r["source"] == "denovo" and get(r)]
            v = np.array([x[0] for x in main]) if main else np.array([])
            ps = [x[1] for x in main if x[1] is not None]
            rec = dict(n=len(main), medoid=med[0] if med else None, medoid_p=med[1] if med else None,
                       median=float(np.median(v)) if len(v) else None,
                       lo=float(np.percentile(v, 10)) if len(v) else None,
                       hi=float(np.percentile(v, 90)) if len(v) else None,
                       same_sign=float(np.mean(np.sign(v) == np.sign(med[0]))) if (len(v) and med) else None,
                       frac_p05=float(np.mean(np.array(ps) < .05)) if ps else None,
                       denovo=[x[0] for x in dn], call=call(med, main))
            out[lab][tf] = rec
    # AR under every definition: rank among regulators at their medoids
    rank = {}
    for lab, get in STATS:
        others = {tf: next((get(r) for r in by[tf] if r["medoid"]), None) for tf in ORDER if tf != "AR"}
        others = {k: v[0] for k, v in others.items() if v}
        rs = []
        for r in by["AR"]:
            v = get(r)
            if not v: continue
            rs.append(dict(matrix=r["matrix"], source=r["source"], value=v[0],
                           rank=1 + sum(1 for x in others.values() if x > v[0]), of=len(others) + 1))
        rank[lab] = rs
    json.dump(dict(stats=out, ar_rank=rank), open(result("sens_summary.json", "regulator_motifs"), "w"),
              indent=1, default=float)

    # ---------------- table
    for lab, _ in STATS:
        print(f"\n  === {lab}")
        print(f"    {'TF':14s}{'n':>4}{'medoid':>9}{'p':>8}{'median':>9}{'10-90%':>18}{'same sign':>10}{'p<.05':>7}  call")
        for tf in ORDER:
            r = out[lab][tf]
            if r["medoid"] is None:
                print(f"    {LABEL.get(tf, tf):14s}{r['n']:>4}   -- no medoid value --"); continue
            pp = f"{r['medoid_p']:.3f}" if r["medoid_p"] is not None else "   -"
            rng_ = f"[{r['lo']:+.3f},{r['hi']:+.3f}]" if r["lo"] is not None else ""
            ss = f"{100 * r['same_sign']:.0f}%" if r["same_sign"] is not None else ""
            fp = f"{100 * r['frac_p05']:.0f}%" if r["frac_p05"] is not None else ""
            print(f"    {LABEL.get(tf, tf):14s}{r['n']:>4}{r['medoid']:>+9.3f}{pp:>8}"
                  f"{(r['median'] if r['median'] is not None else float('nan')):>+9.3f}{rng_:>18}{ss:>10}{fp:>7}  {r['call']}")
        rk = rank[lab]
        if rk:
            rr = [x["rank"] for x in rk]
            print(f"    AR rank among {rk[0]['of']}: best {min(rr)}, median {int(np.median(rr))}, worst {max(rr)}"
                  f" over {len(rk)} AR matrices")

    # ---------------- figure
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.6))
    rngj = np.random.default_rng(3)
    for ax, (lab, get) in zip(axes.flat, STATS):
        for i, tf in enumerate(ORDER):
            y = len(ORDER) - 1 - i
            col = RA if REGUL[tf] == "Glut-RA" else G4
            if tf == "AR": ax.axhspan(y - .45, y + .45, color=RA, alpha=.06, lw=0, zorder=0)
            for r in by[tf]:
                v = get(r)
                if not v: continue
                sig = v[1] is not None and v[1] < .05
                jy = y + rngj.uniform(-.22, .22)
                if r["medoid"]:
                    ax.plot([v[0]], [y], "D", ms=7, color=INK, mec="white", mew=.8, zorder=5)
                elif r["source"] == "homer":
                    ax.plot([v[0]], [jy], "^", ms=5.5, color="#3B6FB6", mfc="#3B6FB6" if sig else "white", mew=1, zorder=4)
                elif r["source"] == "denovo":
                    ax.plot([v[0]], [jy], "s", ms=5, color="#D07A12", mfc="#D07A12" if sig else "white", mew=1, zorder=4)
                else:
                    ax.plot([v[0]], [jy], "o", ms=3.2, color=col, alpha=.55 if sig else .25, mew=0, zorder=3)
        ax.axvline(0, color=INK3, lw=.8, zorder=1)
        ax.set_yticks(range(len(ORDER)))
        ax.set_yticklabels([LABEL.get(t, t) for t in ORDER[::-1]], fontsize=8.5)
        for t, tf in zip(ax.get_yticklabels(), ORDER[::-1]):
            t.set_color(RA if REGUL[tf] == "Glut-RA" else G4)
            if tf == "AR": t.set_fontweight("bold")
        ax.set_title(lab, fontsize=10, color=INK, loc="left", fontweight="bold")
        ax.tick_params(labelsize=8)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        for s in ("left", "bottom"): ax.spines[s].set_color(RULE)
        ax.grid(axis="x", color="#EEF1F3", lw=.7); ax.set_axisbelow(True)
    fig.legend(handles=[Line2D([], [], marker="D", ls="", color=INK, ms=7, label="matrix used in the main analyses (cluster medoid)"),
                        Line2D([], [], marker="o", ls="", color=INK3, ms=5, alpha=.6, label="other members of the SCENIC+ cluster"),
                        Line2D([], [], marker="^", ls="", color="#3B6FB6", ms=6, label="HOMER known"),
                        Line2D([], [], marker="s", ls="", color="#D07A12", ms=6, label="HOMER de novo (found on Glut-RA DARs)")],
               loc="lower center", ncol=4, frameon=False, fontsize=8.8, bbox_to_anchor=(.5, .0))
    fig.suptitle("How much do the regulator comparisons depend on which matrix represents each TF?",
                 fontsize=12.5, fontweight="bold", color=INK, y=.985)
    fig.text(.5, .948, "one point per matrix; filled = p < 0.05 (DoD panels: Glut-RA − GABA-4 permutation test); "
                       "Glut-RA regulators red, GABA-4 teal", ha="center", fontsize=9, color=INK3)
    fig.tight_layout(rect=(0, .04, 1, .935))
    for ext in ("png", "pdf"):
        fig.savefig(result(f"fig_motif_sensitivity.{ext}", "regulator_motifs"), dpi=190)
    print("\n  wrote results/regulator_motifs/fig_motif_sensitivity.png / .pdf and sens_summary.json")


if __name__ == "__main__":
    main()
