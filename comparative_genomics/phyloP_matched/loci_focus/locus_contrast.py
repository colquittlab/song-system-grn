#!/usr/bin/env python
"""phyloP matched-set contrast restricted to the MAFB and EMX2 neighbourhoods.

Element sets come from the bLonStr1 liftover (all group-assigned DARs within
1 Mb of the lifted locus extent), which is what makes these analysable: the
curated locus sets held 56/51 at MAFB and 37/12 at EMX2, and at EMX2 22% of the
curated set was not at the locus at all -- lonStrDom2 carries a ~5.5 Mb
inversion whose breakpoint falls just distal to the EMX2 core, confirmed against
bTaeGut7 and bPasDom1.

delta = log2(rate_oscine) - log2(rate_matched_background); POSITIVE means FASTER
in oscines, i.e. less constrained. interaction = mean(delta|A) - mean(delta|B).
"""
import numpy as np
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"; P=f"{A}/phyloP_matched"
rng=np.random.default_rng(0)
d={};cls={}
for l in open(f"{P}/final/DA_primary_per_element.tsv"):
    f=l.rstrip("\n").split("\t")
    if f[0]=="element_id": continue
    d[f[0]]=float(f[3]); cls[f[0]]=f[1]
SETS={"genome-wide":None}
for nm in ("mafb","emx2"):
    SETS[nm.upper()+" +/-1Mb"]={l.strip() for l in open(f"{A}/vgp577/lifted/{nm}_1mb.txt")}
A_,B_="RA_over_Arco","GABA4_over_GABA1"
print(f"  {'set':<18}{'nRA':>6}{'nG4':>6}{'mean d RA':>12}{'mean d G4':>12}"
      f"{'interaction':>13}{'95% CI':>22}{'p_perm':>9}")
out={}
for nm,keep in SETS.items():
    ea=[e for e in d if cls[e]==A_ and (keep is None or e in keep)]
    eb=[e for e in d if cls[e]==B_ and (keep is None or e in keep)]
    if len(ea)<10 or len(eb)<10: continue
    va=np.array([d[e] for e in ea]); vb=np.array([d[e] for e in eb])
    obs=va.mean()-vb.mean()
    bs=[np.mean(rng.choice(va,len(va)))-np.mean(rng.choice(vb,len(vb))) for _ in range(4000)]
    pool=np.concatenate([va,vb]); n=len(va)
    perm=[]
    for _ in range(4000):
        p=rng.permutation(pool); perm.append(p[:n].mean()-p[n:].mean())
    perm=np.array(perm); pv=(np.abs(perm)>=abs(obs)).mean()
    lo,hi=np.percentile(bs,[2.5,97.5])
    out[nm]=(len(ea),len(eb),va.mean(),vb.mean(),obs,lo,hi,pv)
    print(f"  {nm:<18}{len(ea):>6}{len(eb):>6}{va.mean():>+12.5f}{vb.mean():>+12.5f}"
          f"{obs:>+13.5f}   [{lo:+.5f},{hi:+.5f}]{max(pv,1/4000):>9.4f}")
import json; json.dump(out,open(f"{P}/loci_focus/results/compare_gw_mafb_emx2/locus_contrast.json","w"),indent=1)
