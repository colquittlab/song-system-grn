#!/usr/bin/env python
"""One tree per switch column: every site/position where the oscine and
outgroup majority bases differ with >=80% support on each side."""
import json, sys, collections
sys.path.insert(0,"/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/maf_switch")
from treelib import draw_tree
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import per_site  # results/per_site/
D="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/maf_switch"
MINO,MING,THR=int(sys.argv[1]) if len(sys.argv)>1 else 40, 60, 0.80
COMP=str.maketrans("ACGTacgtNn-","TGCAtgcaNn-")
OSC={l.strip() for l in open(f"{D}/data/osc363.txt") if l.strip()}
tree=open(f"{D}/data/tree363.nw").read().strip()
cols={int(k):v for k,v in json.load(open(f"{D}/data/cols363.json")).items()}
R=json.load(open(f"{D}/data/sites363.json"))
idx=[]
for r in R:
    for c in r["cols"]:
        if c["oscN"]<MINO or c["outN"]<MING: continue
        if c["osc"]!=c["out"] and c["oscF"]>=THR and c["outF"]>=THR:
            idx.append((r,c))
print(f"  switch columns to draw: {len(idx)}")
man=[]
for r,c in idx:
    d=dict(cols.get(c["refpos"],{}))
    if r["strand"]=="-": d={s:b.translate(COMP) for s,b in d.items()}
    nm=f"{r['motif']}_{r['site']}_p{c['k']}"
    stem=f"{per_site('switch_trees_363', nm)}/{nm}"
    t=(f"{r['site']}  ({'ARE-full' if r['motif']=='full' else 'JASPAR AR'} position {c['k']})"
       f"  chr20:{c['refpos']+1:,}  {c['out']}→{c['osc']}")
    N,coll,ns=draw_tree(d,tree,stem,OSC,t)
    man.append(dict(file=nm+".pdf",motif=r["motif"],site=r["site"],pos=c["k"],
                    refpos=c["refpos"],change=f"{c['out']}>{c['osc']}",
                    oscF=c["oscF"],outF=c["outF"],oscN=c["oscN"],outN=c["outN"],
                    tips=N,collapsed=coll,substitutions=ns))
    print(f"    {nm:<46} {c['out']}>{c['osc']}  tips {N:3d}  subs {ns}")
json.dump(man,open(f"{D}/trees/manifest.json","w"),indent=1)
print(f"  manifest: {D}/trees/manifest.json")
