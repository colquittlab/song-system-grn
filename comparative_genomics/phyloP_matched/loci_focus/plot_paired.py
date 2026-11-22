#!/usr/bin/env python
"""AR motif minus the remainder of its own enhancer, on both trees.

Pairing within element removes every between-element difference -- GC, repeat
content, alignability, cell-type identity -- so what remains is whether the
binding site is held more tightly than the sequence it sits in.
"""
import numpy as np, collections, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
    "mathtext.it":"Arial:italic","mathtext.bf":"Arial:bold",
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2,
    "legend.handletextpad":.4,"legend.borderpad":.2,"legend.labelspacing":.3})
import matplotlib.pyplot as plt
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
rng=np.random.default_rng(0)
INK="#171B22"; INK2="#454E5C"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
def load(f):
    d={}
    for l in open(f"{F}/{f}"):
        p=l.split("\t"); d[p[0]]=(int(p[1]),int(p[2]),float(p[5]))
    return d
coord={l.split("\t")[3]:1 for l in open(f"{F}/allelem_cls.bed")}
h2e={}; n=0
for s in ("ra","g4"):
    for i,l in enumerate(open(f"{F}/fimo_motif2_{s}.txt")):
        if not i: continue
        e=l.split("\t")[2]
        if e not in coord: continue
        n+=1; h2e[n]=e
def diffs(tree):
    mot=load(f"{tree}_motif2.tab"); rest=load(f"{tree}_gwrest.tab")
    rb=collections.defaultdict(list)
    for k,v in rest.items():
        if v[1]/max(v[0],1)>=.5: rb[k.split("|")[1]].append((v[2],v[0]))
    out=collections.defaultdict(list)
    for k,v in mot.items():
        if v[1]/max(v[0],1)<.5: continue
        p=k.split("|"); e=h2e.get(int(p[3]))
        if not e or e not in rb: continue
        w=np.array(rb[e]); out[(p[1],p[2])].append(v[2]-float(np.average(w[:,0],weights=w[:,1])))
    return out
def ci(d):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(4000)])
    return np.percentile(b,[2.5,97.5])
def pval(d,n=20000):
    b=np.array([np.mean(rng.choice(d,len(d))) for _ in range(n)])
    return max(2*min((b<=0).mean(),(b>=0).mean()),1/n)
fig,axes=plt.subplots(1,2,figsize=(2.30,1.13))
fig.patch.set_facecolor("white")
GR=[("genome-wide",None),("at MAFB","MAFB")]
for ax,(tree,ttl) in zip(axes,(("subtree222","A   oscine subtree"),
                               ("wholetree","B   whole tree (pan-avian)"))):
    D=diffs(tree); NLAB={}
    for i,(glab,loc) in enumerate(GR):
        for j,(cls,col,clab) in enumerate((("ra",RA,"Glut-RA"),("g4",G4,"GABA-4"))):
            d=np.array([x for (c,l),v in D.items() if c==cls and (loc is None or l==loc) for x in v])
            if len(d)<8: continue
            lo,hi=ci(d); x=i+(j-.5)*.22
            ax.plot([x,x],[lo,hi],color=col,lw=1.30,zorder=2,solid_capstyle="round")
            ax.scatter([x],[d.mean()],s=13.8,color=col,zorder=3,linewidths=0)
            pv=pval(d); star="*" if pv<.05 else ""
            dx,ha=(-.06,"right") if j==0 else (.06,"left")
            ax.text(x+dx,d.mean(),f"{d.mean():+.3f}{star}",ha=ha,va="center",
                    fontsize=5,color=col,fontweight="600" if pv<.05 else "400")
            NLAB.setdefault(i,[]).append(len(d))
    for gi,v in NLAB.items():
        ax.text(gi,-.15," / ".join(f"{q:,}" for q in v),ha="center",va="top",
                fontsize=5,color=INK3,transform=ax.get_xaxis_transform())
    ax.axhline(0,color=INK,lw=0.50,zorder=1)
    ax.set_xticks(range(len(GR))); ax.set_xticklabels([g[0] for g in GR],fontsize=5,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600",pad=10)
    ax.set_xlim(-.42,1.42)
    ax.margins(x=.03,y=.05)
    ax.locator_params(nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=INK3,labelsize=5)
axes[0].set_ylabel("motif − enhancer",fontsize=5,color=INK)
axes[0].set_ylim(-.30,.40); axes[1].set_ylim(-1.00,.62)
from matplotlib.patches import Patch
axes[1].legend(handles=[Patch(facecolor=RA,label="Glut-RA DARs"),
                        Patch(facecolor=G4,label="GABA-4 DARs")],
               frameon=False,fontsize=5,loc="lower left",borderaxespad=.1)
fig.subplots_adjust(left=.115,right=.985,top=.90,bottom=.235,wspace=.16)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_paired.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_ar_paired.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_ar_paired.png")
