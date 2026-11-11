#!/usr/bin/env python
"""Is the AR-site gain rate elevated on the oscine stem?

Presence/absence is read from every node's own sequence, so gains and losses
are observed on branches. Rates are events per unit branch length, summed over
sites.

Ascertainment caveat, and the control for it: sites were called on
lonStrDom2, so only sites present in that genome exist in this set. Every
branch on the root-to-Lonchura path therefore shares the same constraint --
a gain there is over-represented and a loss there is impossible. Comparing the
oscine stem to OTHER BRANCHES ON THAT PATH holds the ascertainment constant;
comparing it to the whole tree does not.
"""
import sys, collections, json, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"..")
from span import elem_coords, in_domain
SPAN=sys.argv[1] if len(sys.argv)>1 else "domain"   # "domain" (default) or "1p5mb"
EL=elem_coords()
from scipy import stats
from treelib import parse_newick, leaves
from gain_trees import F
OSC={l.strip() for l in open("data74/osc74.txt") if l.strip()}
tree=open("data74/tree74.nw").read().strip()
root=parse_newick(tree)
par={}; node={}
def idx(n):
    node[n["name"]]=n
    for c in n["ch"]: par[c["name"]]=n["name"]; idx(c)
idx(root)
# the oscine stem: branch entering the smallest clade holding the songbirds
best_c=[None]
def find(n):
    L={x["name"] for x in leaves(n)}
    if {"Taeniopygia_guttata","Corvus_hawaiiensis","Lonchura_striata"}<=L and "Acanthisitta_chloris" not in L:
        if best_c[0] is None or len(L)<len(best_c[0][1]): best_c[0]=(n,L)
    for c in n["ch"]: find(c)
find(root)
STEM=best_c[0][0]["name"]
print(f"  oscine stem branch -> node {STEM!r}  ({len(best_c[0][1])} tips, length {node[STEM]['bl']:.4f})")
# branches on the root -> Lonchura path
path=[]; cur="Lonchura_striata"
while cur in par: path.append(cur); cur=par[cur]
PATH=set(path)
print(f"  branches on the root→Lonchura path: {len(PATH)}   stem on path: {STEM in PATH}")
# per-site presence: the corrected call -- reference projection, and a node with
# no orthologous sequence across the motif is missing rather than absent
import presence as PR
from rebuild_gain import call
B7=PR.load("data74/ar_windows_gap_74.maf")
sites=[]
for bed,mo in ((f"{F}/tracks/mafb_1p5mb/arefull_sites_MAFB_lsd2.bed","full"),
               (f"{F}/tracks/mafb_1p5mb/jasparAR_sites_MAFB_lsd2.bed","jaspar")):
    for l in open(bed):
        x=l.rstrip("\n").split("\t")
        if SPAN=="domain" and not in_domain(EL[x[10]]): continue
        sites.append(dict(motif=mo,name=x[3],s=int(x[1]),e=int(x[2]),cls=x[9]))
G=collections.Counter(); L_=collections.Counter(); BL=collections.Counter(); NS=0
for st in sites:
    pres=call(B7,st)
    ext=[g for g in pres if not g.startswith("Anc") and pres[g] is not None]
    if len(ext)<50: continue
    NS+=1
    def rec(n):
        for c in n["ch"]:
            a,b2=pres.get(n["name"]),pres.get(c["name"])
            if a is not None and b2 is not None:
                BL[c["name"]]+=c["bl"]
                if a!=b2: (G if b2 else L_)[c["name"]]+=1
            rec(c)
    rec(root)
print(f"  span: {SPAN}   sites with enough coverage: {NS}")
def rate(keys):
    g=sum(G[k] for k in keys); l=sum(L_[k] for k in keys); b=sum(BL[k] for k in keys)
    return g,l,b,(g/b if b else np.nan),(l/b if b else np.nan)
allk=[k for k in BL]
stem=[STEM]; pathk=[k for k in PATH if k in BL]; pathoth=[k for k in pathk if k!=STEM]
oth=[k for k in allk if k not in PATH]
for lab,keys in (("oscine stem",stem),("rest of the root→Lonchura path",pathoth),
                 ("all other branches",oth),("whole tree",allk)):
    g,l,b,gr,lr=rate(keys)
    print(f"  {lab:<34} branches {len(keys):>4}  gains {g:>4}  losses {l:>4}"
          f"  length {b:>8.2f}  gain rate {gr:>7.2f}  loss rate {lr:>7.2f}")
print()
for lab,keys in (("vs rest of the path (ascertainment-matched)",pathoth),
                 ("vs all other branches (ascertainment-inflated)",oth)):
    g0,l0,b0,_,_=rate(stem); g1,l1,b1,_,_=rate(keys)
    exp=b0*(g1/b1)
    p=stats.poisson.sf(g0-1,exp)
    print(f"  {lab}")
    print(f"    oscine stem: {g0} gains over length {b0:.2f}"
          f"   expected {exp:.2f} at the comparison rate   Poisson p = {p:.4f}"
          f"   fold {g0/exp if exp else float('nan'):.2f}x")
json.dump({"stem":STEM,"gains":{k:int(v) for k,v in G.items()},
           "losses":{k:int(v) for k,v in L_.items()},
           "bl":{k:float(v) for k,v in BL.items()},"path":path,"n_sites":NS},
          open(f"data74/rate_test_{SPAN}.json","w"),indent=1)
