#!/usr/bin/env python
"""Motif scoring per TF, using the motifs SCENIC+ itself assigns.

Motif sets come from the v10nr_clust annotation used by config4: Direct_annot
plus Orthology_annot, matching that run's `annotations_to_use`. Each assigned
motif is scored independently with the symmetric-ascertainment statistic
(forward = oscine-ascertained windows, reverse = outgroup-ascertained, midpoint
= winner's-curse-corrected), and the TF's value is the mean across its motifs.

Consensus sequences are built once per element and scanned vectorised, so the
~3,000 matrices are tractable.
"""
import sys, csv, collections, json, glob, os, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, OSC, OUT, cls, dom, MINO, MING
V="/private/groups/colquittlab/scenicplus/v10nr_clust_public"
TBL=f"{V}/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"
PCT=99.5; PC=0.001
REG={"KLF12":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA","LIN28B":"Glut-RA",
     "STAT1":"Glut-RA","FOSL2":"Glut-RA","AR":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
want=collections.defaultdict(set)
for r in csv.DictReader(open(TBL),delimiter="\t"):
    g=r["gene_name"]
    if g in REG: want[g].add(r["#motif_id"])        # direct + orthology, as config4 used
print("  motifs assigned per TF:",{k:len(v) for k,v in sorted(want.items())})
need={m for s in want.values() for m in s}
def read_cb(p):
    rows=[]
    for l in open(p):
        if l.startswith(">"): continue
        q=l.split()
        if len(q)==4:
            try: rows.append([float(x) for x in q])
            except ValueError: pass
    if not rows: return None
    a=np.array(rows); a=a/a.sum(1,keepdims=True)
    return np.log2((a+PC)/(1+4*PC)/.25)
PWM={}
for m in need:
    p=f"{V}/singletons/{m}.cb"
    if os.path.exists(p):
        w=read_cb(p)
        if w is not None and 4<=len(w)<=30: PWM[m]=w
print(f"  matrices loaded: {len(PWM):,} of {len(need):,} assigned "
      f"(missing or out of width range: {len(need)-len(PWM):,})")
# ---- consensus sequence per element, as integer arrays
col=load_cols("data/domain_363.maf")
IDX={"A":0,"C":1,"G":2,"T":3}
seqs=[]
for e in dom:
    ch,es,ee,cl=cls[e]
    o=np.full(ee-es,-1,np.int8); g=np.full(ee-es,-1,np.int8)
    for i,rp in enumerate(range(es,ee)):
        d=col.get(rp)
        if not d: continue
        co=collections.Counter(b for s,b in d.items() if s in OSC)
        cg=collections.Counter(b for s,b in d.items() if s in OUT)
        if sum(co.values())>=MINO: o[i]=IDX[co.most_common(1)[0][0]]
        if sum(cg.values())>=MING: g[i]=IDX[cg.most_common(1)[0][0]]
    seqs.append((e,"Glut-RA" if cl=="RA_over_Arco" else "GABA-4",o,g))
print(f"  consensus built for {len(seqs)} elements")
RC=np.array([3,2,1,0],np.int8)
def scan(seq,lo):
    """max of both strands at every valid offset; nan where any base is missing"""
    W=len(lo); n=len(seq)-W+1
    if n<=0: return np.array([])
    idx=np.lib.stride_tricks.sliding_window_view(seq,W)
    ok=(idx>=0).all(1)
    safe=np.where(idx>=0,idx,0)
    f=np.zeros(len(idx)); r=np.zeros(len(idx))
    for k in range(W):
        f+=lo[k][safe[:,k]]
        r+=lo[W-1-k][RC[safe[:,k]]]
    out=np.maximum(f,r); out[~ok]=np.nan
    return out
res={}
for ti,(tf,ms) in enumerate(sorted(want.items())):
    vals=[]
    for m in sorted(ms):
        if m not in PWM: continue
        lo=PWM[m]; rec={"Glut-RA":[],"GABA-4":[]}
        for e,kls,o,g in seqs:
            a=scan(o,lo); b=scan(g,lo)
            if a.size==0: continue
            ok=~(np.isnan(a)|np.isnan(b))
            if ok.sum(): rec[kls].append(np.c_[a[ok],b[ok]])
        allv=np.concatenate([x for k in rec for x in rec[k]]) if any(rec.values()) else None
        if allv is None or len(allv)<100: continue
        th=np.percentile(allv.ravel(),PCT)
        mids={}
        for kls in ("Glut-RA","GABA-4"):
            if not rec[kls]: continue
            A=np.concatenate(rec[kls]); d=A[:,0]-A[:,1]
            f=d[A[:,0]>=th]; v=d[A[:,1]>=th]
            if len(f)>=4 and len(v)>=4: mids[kls]=(f.mean()+v.mean())/2
        if len(mids)==2: vals.append((mids["Glut-RA"],mids["GABA-4"]))
    if len(vals)<5: print(f"  {tf:<8} too few usable matrices ({len(vals)})"); continue
    a=np.array(vals)
    res[tf]=dict(reg=REG[tf],n_motif=len(vals),n_assigned=len(ms),
                 ra=float(a[:,0].mean()),g4=float(a[:,1].mean()),
                 diff=float((a[:,0]-a[:,1]).mean()),
                 sd=float((a[:,0]-a[:,1]).std()),
                 se=float((a[:,0]-a[:,1]).std()/np.sqrt(len(a))))
    print(f"    {tf} done ({len(vals)} matrices)",flush=True)
print(f"\n  {'TF':<9}{'regulon':<9}{'motifs':>8}{'RA mean':>10}{'G4 mean':>10}"
      f"{'RA-G4':>9}{'sd':>8}{'se':>8}")
for tf,d in sorted(res.items(),key=lambda kv:-kv[1]["diff"]):
    print(f"  {tf:<9}{d['reg']:<9}{d['n_motif']:>4}/{d['n_assigned']:<3}{d['ra']:>+10.2f}"
          f"{d['g4']:>+10.2f}{d['diff']:>+9.2f}{d['sd']:>8.2f}{d['se']:>8.2f}")
json.dump(res,open("data/ascertain_tf_scenic.json","w"),indent=1)
