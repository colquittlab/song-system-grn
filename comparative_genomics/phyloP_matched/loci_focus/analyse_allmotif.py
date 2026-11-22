#!/usr/bin/env python
"""Paired motif-vs-own-enhancer constraint for 192 motif archetypes.

Same statistic as the AR analysis: for each motif instance set within a DAR,
the motif's phyloP minus the phyloP of the remainder of that same element. The
remainder is obtained by subtraction (DAR total minus the motif's own merged
hits), so every motif is compared against an identically-defined local
background and the values are directly rankable.

The reference still contains other motifs -- the masking curve showed that is
unavoidable -- so these are relative contrasts, not constraint-versus-neutral.
"""
import numpy as np, collections, sys
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
TREE=sys.argv[1] if len(sys.argv)>1 else "subtree222"
rng=np.random.default_rng(0)
cls={}
for l in open(f"{F}/allelem_cls.bed"):
    f=l.rstrip("\n").split("\t"); cls[f[3]]=f[4]
dar={}
for l in open(f"{F}/{TREE}_am_dars.tab"):
    p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))   # size, covered, sum
agg=collections.defaultdict(lambda:[0,0,0.0])                      # (motif,elem) -> size,cov,sum
for l in open(f"{F}/{TREE}_am_hits.tab"):
    p=l.split("\t"); m,e,_=p[0].split("@")
    a=agg[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
print(f"  tree={TREE}  DARs scored={len(dar):,}  (motif,element) pairs={len(agg):,}",flush=True)
per=collections.defaultdict(lambda: collections.defaultdict(list))
for (m,e),(msz,mcov,msum) in agg.items():
    if e not in dar or e not in cls: continue
    dsz,dcov,dsum=dar[e]
    rcov=dcov-mcov; rsum=dsum-msum
    if mcov<0.5*msz or rcov<0.5*(dsz-msz) or rcov<50: continue
    per[m][cls[e]].append(msum/mcov - rsum/rcov)
def boot_p(d,n=10000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
def bh(ps):
    o=np.argsort(ps); r=np.empty(len(ps)); s=np.array(ps)[o]
    q=s*len(ps)/np.arange(1,len(ps)+1); q=np.minimum.accumulate(q[::-1])[::-1]
    r[o]=np.clip(q,0,1); return r
rows=[]
for m,d in per.items():
    ra=np.array(d.get("RA_over_Arco",[])); g4=np.array(d.get("GABA4_over_GABA1",[]))
    if len(ra)<40 or len(g4)<40: continue
    both=np.concatenate([ra,g4])
    rows.append(dict(motif=m,n_ra=len(ra),n_g4=len(g4),
                     d_both=both.mean(),p_both=boot_p(both),
                     d_ra=ra.mean(),p_ra=boot_p(ra),d_g4=g4.mean(),
                     d_diff=ra.mean()-g4.mean()))
for k in ("p_both","p_ra"):
    q=bh([r[k] for r in rows])
    for r,v in zip(rows,q): r["q_"+k.split("_")[1]]=v
rows.sort(key=lambda r:-r["d_both"])
print(f"  motifs with >=40 elements in both classes: {len(rows)}\n")
def short(m):
    return m.split(".",2)[-1].replace("jaspar__","").replace("hocomoco__","")\
            .replace("swissregulon__mm__","").replace("swissregulon__hs__","")\
            .replace("cisbp__","").replace("transfac_public__","")[:26]
print(f"  {'rank':>5}  {'motif':<28}{'n RA':>7}{'n G4':>7}{'both':>9}{'q':>9}"
      f"{'RA':>9}{'GABA-4':>9}{'RA-G4':>9}")
AR=[i for i,r in enumerate(rows) if "MA0007.3" in r["motif"]]
for i,r in enumerate(rows):
    tag="  <- AR" if i in AR else ""
    if i<12 or i>=len(rows)-6 or i in AR:
        print(f"  {i+1:>5}  {short(r['motif']):<28}{r['n_ra']:>7}{r['n_g4']:>7}"
              f"{r['d_both']:>+9.3f}{r['q_both']:>9.4f}{r['d_ra']:>+9.3f}{r['d_g4']:>+9.3f}"
              f"{r['d_diff']:>+9.3f}{tag}")
    elif i==12: print(f"  {'...':>5}")
sig=[r for r in rows if r["q_both"]<0.05]
print(f"\n  significant at FDR 5% (both classes pooled): {len(sig)} / {len(rows)}")
print(f"    positive (motif > its enhancer): {sum(1 for r in sig if r['d_both']>0)}")
print(f"    negative: {sum(1 for r in sig if r['d_both']<0)}")
if AR:
    r=rows[AR[0]]
    pct=100*(1-AR[0]/len(rows))
    print(f"\n  AR (MA0007.3): rank {AR[0]+1}/{len(rows)} ({pct:.0f}th percentile)"
          f"  d={r['d_both']:+.3f}  q={r['q_both']:.4f}")
    dd=np.array([x["d_diff"] for x in rows])
    rk=int((dd>r["d_diff"]).sum())+1
    print(f"  AR RA-minus-GABA4 contrast {r['d_diff']:+.3f}: rank {rk}/{len(rows)}"
          f"  (distribution median {np.median(dd):+.3f}, sd {dd.std(ddof=1):.3f})")
import json
json.dump(rows,open(f"{F}/allmotif_{TREE}.json","w"),indent=1)
