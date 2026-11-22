#!/usr/bin/env python
"""Per-TF motif scoring on the motif sets SCENIC+ actually searched.

Motif set per TF = its annotation in motifs-v10-nr.hgnc (Direct_annot plus
Orthology_annot at >=70% orthologous identity) intersected with the 10,249
motifs in the v10nr_clust collection, which is exactly what the cistarget
database contains. Those are the matrices that could have produced the
eRegulon calls.

Statistic per motif: symmetric ascertainment inside the DHX35-TOP1 domain --
forward (oscine-ascertained windows), reverse (outgroup-ascertained), midpoint
= winner's-curse-corrected. The TF value is the mean across its matrices, with
a bootstrap over matrices for the spread.
"""
import sys, csv, collections, json, os, glob, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from ascertain_sym import load_cols, OSC, OUT, cls, dom, MINO, MING
V="/private/groups/colquittlab/scenicplus/v10nr_clust_public"
TBL=f"{V}/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"
PCT=99.5; PC=0.001; OID=0.70
REG={"KLF12":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA","LIN28B":"Glut-RA",
     "STAT1":"Glut-RA","FOSL2":"Glut-RA","AR":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
CT={os.path.basename(p)[:-3] for p in glob.glob(f"{V}/singletons/*.cb")}
want=collections.defaultdict(set)
for r in csv.DictReader(open(TBL),delimiter="\t"):
    g=r["gene_name"]
    if g not in REG: continue
    if "directly annotated" in r["description"] and "orthologous" not in r["description"]:
        want[g].add(r["#motif_id"])
    else:
        try: oi=float(r["orthologous_identity"])
        except (ValueError,TypeError): oi=0.0
        if oi>=OID: want[g].add(r["#motif_id"])
want={g:(s&CT) for g,s in want.items()}
print("  motifs per TF in the searched collection:",
      {k:len(v) for k,v in sorted(want.items())})
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
for m in {x for s in want.values() for x in s}:
    w=read_cb(f"{V}/singletons/{m}.cb")
    if w is not None and 4<=len(w)<=30: PWM[m]=w
print(f"  matrices usable (width 4-30): {len(PWM):,}")
col=load_cols("data/domain_363.maf"); IDX={"A":0,"C":1,"G":2,"T":3}
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
    seqs.append(("Glut-RA" if cl=="RA_over_Arco" else "GABA-4",o,g))
RC=np.array([3,2,1,0],np.int8)
def scan(seq,lo):
    W=len(lo)
    if len(seq)-W+1<=0: return None
    idx=np.lib.stride_tricks.sliding_window_view(seq,W)
    ok=(idx>=0).all(1); safe=np.where(idx>=0,idx,0)
    f=np.zeros(len(idx)); r=np.zeros(len(idx))
    for k in range(W):
        f+=lo[k][safe[:,k]]; r+=lo[W-1-k][RC[safe[:,k]]]
    out=np.maximum(f,r); out[~ok]=np.nan
    return out
def midpoints(lo):
    rec={"Glut-RA":[],"GABA-4":[]}
    for kls,o,g in seqs:
        a=scan(o,lo); b=scan(g,lo)
        if a is None: continue
        ok=~(np.isnan(a)|np.isnan(b))
        if ok.sum(): rec[kls].append(np.c_[a[ok],b[ok]])
    if not rec["Glut-RA"] or not rec["GABA-4"]: return None
    allv=np.concatenate([x for k in rec for x in rec[k]])
    th=np.percentile(allv.ravel(),PCT); out={}
    for kls in ("Glut-RA","GABA-4"):
        A=np.concatenate(rec[kls]); d=A[:,0]-A[:,1]
        f=d[A[:,0]>=th]; v=d[A[:,1]>=th]
        if len(f)<4 or len(v)<4: return None
        out[kls]=(f.mean()+v.mean())/2
    return out["Glut-RA"],out["GABA-4"]
rng=np.random.default_rng(0); res={}
for tf,ms in sorted(want.items()):
    vals=[midpoints(PWM[m]) for m in sorted(ms) if m in PWM]
    vals=[v for v in vals if v]
    if not vals: print(f"  {tf}: no usable matrices"); continue
    a=np.array(vals); diff=a[:,0]-a[:,1]
    bs=[np.mean(rng.choice(diff,len(diff),replace=True)) for _ in range(4000)] if len(diff)>1 else []
    ci=np.percentile(bs,[2.5,97.5]) if bs else [np.nan,np.nan]
    res[tf]=dict(reg=REG[tf],n=len(vals),n_assigned=len(ms),ra=float(a[:,0].mean()),
                 g4=float(a[:,1].mean()),diff=float(diff.mean()),
                 ci=[float(ci[0]),float(ci[1])])
print(f"\n  {'TF':<9}{'regulon':<9}{'matrices':>10}{'RA mean':>10}{'G4 mean':>10}"
      f"{'RA-G4':>9}{'95% CI (over matrices)':>26}")
for tf,d in sorted(res.items(),key=lambda kv:-kv[1]["diff"]):
    ci=f"[{d['ci'][0]:+.2f}, {d['ci'][1]:+.2f}]" if d["n"]>1 else "single matrix"
    star="*" if d["n"]>1 and (d["ci"][0]>0 or d["ci"][1]<0) else ""
    print(f"  {tf:<9}{d['reg']:<9}{d['n']:>10}{d['ra']:>+10.2f}{d['g4']:>+10.2f}"
          f"{d['diff']:>+9.2f}{ci:>26}{star}")
json.dump(res,open("data/ascertain_tf_ct.json","w"),indent=1)
