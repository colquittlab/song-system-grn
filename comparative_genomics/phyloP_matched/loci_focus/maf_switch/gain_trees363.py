#!/usr/bin/env python
"""Per-site gain/loss trees for the regulator (TF) sites on the avian363 tree --
the avian363 counterpart of results/per_site/gain_trees_tf (74-way, rebuild_gain.py).

Presence at every node, tips and cactus's reconstructed ancestors, is the corrected
call (presence.py projection, missing below 80% motif coverage), the same call
tiers363.py places gains with; no event is drawn on a branch touching a missing node.

Which sites get a tree. The 74-way rule -- oscine-enriched (>= 60% of oscines,
<= 35% of non-oscines) with at most two gains -- does not transfer: avian363 has
five times the branches, and homoplasy is reported there, not used to discard
sites. A tree is drawn for every site that passes the same frequency rule, and for
every site whose principal oscine gain is at or near the oscine stem (tiers363).

Uniform clades -- every node, tips and ancestors, in one state -- are collapsed,
oscine clades included, and labelled with their genus range.

Writes results/per_site/gain_trees_tf_363/<span>/<motif>_<site>.pdf/.png and
trees363/manifest_tf.json.
    python gain_trees363.py [AR|TF]      (default TF)
"""
import collections, json, sys
sys.path.insert(0, "."); sys.path.insert(0, "..")
import tiers363 as T
import presence as PR
from gain_trees import draw
from span import per_site

FAM = (sys.argv[1:] or ["TF"])[0]
MAF = {"AR": "data/ar_windows_gap.maf", "TF": "data/tfsc_windows_gap_anc.maf"}[FAM]
COLL = {"AR": "gain_trees_ar_363", "TF": "gain_trees_tf_363"}[FAM]
LAB = {"full": "ARE-full", "jaspar": "JASPAR AR"}
TREE = open("data/tree363.nw").read().strip()


def main():
    cnt = {(c["motif"], c["name"]): c for c in json.load(open(f"trees363/counts_{FAM.lower()}.json"))}
    B = PR.load(MAF)
    man = []
    for st in T.sites(FAM):
        c = cnt.get((st["motif"], st["name"]))
        if c is None: continue                       # below the coverage floor
        pf = c["f_osc"] >= .6 and c["f_out"] <= .35
        if not (pf or c["tier"] in ("stem", "near-stem")): continue
        pres = T.call(B, st)
        o = [g for g in T.OSC if pres.get(g) is not None]; u = [g for g in T.OUT if pres.get(g) is not None]
        where = ("stem" if c["tier"] == "stem" else
                 f'near-stem ({c["depth"]} node{"s" if c["depth"] > 1 else ""} in)' if c["tier"] == "near-stem"
                 else c["tier"])
        nm = f'{st["motif"]}_{st["name"]}'
        N, g, l = draw(pres, TREE, f"{per_site(COLL, nm)}/{nm}", T.OSC,
                       f'{st["name"]}  ({LAB.get(st["motif"], st["motif"])})  chr20:{st["s"] + 1:,}-{st["e"]:,}'
                       f'   oscine {sum(pres[x] for x in o)}/{len(o)}  outgroup {sum(pres[x] for x in u)}/{len(u)}'
                       f'   gain: {where}; {c["nonosc_gains"]} outside oscines',
                       collapse_all=True)
        assert (g, l) == (c["gains"], c["losses"]), (nm, g, l, c["gains"], c["losses"])
        man.append(dict(file=nm + ".pdf", motif=st["motif"], name=st["name"], cls=st["cls"], elem=st["elem"],
                        dom=c["dom"], f_osc=c["f_osc"], f_out=c["f_out"], pass_freq=pf, tier=c["tier"],
                        depth=c["depth"], nonosc_gains=c["nonosc_gains"], tips=N, gains=g, losses=l))
        print(f"    {nm:<46} {c['tier']:9s} osc {c['f_osc']:.2f} out {c['f_out']:.2f}  {g} gains {l} losses  {N} rows",
              flush=True)
    json.dump(man, open(f"trees363/manifest_{FAM.lower()}.json", "w"), indent=1)
    t = collections.Counter(m["tier"] for m in man)
    print(f"  {FAM}: {len(man)} trees   " + ", ".join(f"{k} {v}" for k, v in sorted(t.items())))


if __name__ == "__main__":
    main()
