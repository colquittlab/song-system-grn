#!/usr/bin/env python
"""SCENIC+ MAFB regulators: does each TF's motif show the Glut-RA / GABA-4
asymmetry its regulon assignment predicts?"""
import numpy as np, collections, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "mathtext.fontset":"custom","mathtext.rm":"Arial",
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2,
    "legend.handletextpad":.4,"legend.borderpad":.2,"legend.labelspacing":.3})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
A="/private/warm-archive/colquittlab/assembly/cactus/avian363"
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"; GREY="#9AA3AF"
REG={"KLF12":"Glut-RA","KLF15":"Glut-RA","NR3C2":"Glut-RA","NFIA":"Glut-RA",
     "LIN28B":"Glut-RA","STAT1":"Glut-RA","FOSL2":"Glut-RA",
     "MAF":"GABA-4","SOX6":"GABA-4","TCF12":"GABA-4","ARX":"GABA-4",
     "SATB1":"GABA-4","MEIS2":"GABA-4"}
cls={l.split("\t")[3]:l.rstrip("\n").split("\t")[4] for l in open(f"{F}/allelem_cls.bed")}
def load(tree):
    dar={}
    for l in open(f"{F}/{tree}_am_dars.tab"):
        p=l.split("\t"); dar[p[0]]=(int(p[1]),int(p[2]),float(p[3]))
    agg=collections.defaultdict(lambda:[0,0,0.0])
    for l in open(f"{F}/{tree}_tf13_hits.tab"):
        p=l.split("\t"); m,e,_=p[0].split("@")
        a=agg[(m,e)]; a[0]+=int(p[1]); a[1]+=int(p[2]); a[2]+=float(p[3])
    out=collections.defaultdict(list)
    for (m,e),(msz,mcov,msum) in agg.items():
        if e not in dar or e not in cls: continue
        dsz,dcov,dsum=dar[e]
        rcov=dcov-mcov; rsum=dsum-msum
        if mcov<0.5*msz or rcov<50: continue
        out[(m,cls[e])].append(msum/mcov-rsum/rcov)
    return out
def dodci(a,b,n=3000):
    bs=[np.mean(rng.choice(a,len(a)))-np.mean(rng.choice(b,len(b))) for _ in range(n)]
    return a.mean()-b.mean(), *np.percentile(bs,[2.5,97.5])
AR={"subtree222":0.023,"wholetree":0.086}
fig,axes=plt.subplots(1,2,figsize=(4.9,1.95),sharey=True)
fig.patch.set_facecolor("white")
order=None
for ax,(tree,ttl) in zip(axes,(("subtree222","A   oscine subtree"),
                               ("wholetree","B   whole tree"))):
    D=load(tree); vals={}
    for m in REG:
        a=np.array(D.get((m,"RA_over_Arco"),[])); b=np.array(D.get((m,"GABA4_over_GABA1"),[]))
        if len(a)<40 or len(b)<40: continue
        vals[m]=dodci(a,b)
    if order is None: order=sorted(vals,key=lambda m:-vals[m][0])
    for i,m in enumerate(order):
        d,lo,hi=vals[m]; y=len(order)-1-i
        col=RA if REG[m]=="Glut-RA" else G4
        ax.plot([lo,hi],[y]*2,color=col,lw=.7,alpha=.75,solid_capstyle="round")
        ax.scatter([d],[y],s=9,color=col,zorder=3,linewidths=0)
    ax.axvline(0,color=INK,lw=.5,zorder=1)
    ax.axvline(AR[tree],color=GREY,lw=.7,ls=(0,(2,2)),zorder=1)
    ax.text(AR[tree],len(order)-.35,"AR",ha="center",va="bottom",fontsize=5,color=GREY)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1],fontsize=5,color=INK)
    ax.set_xlabel("Glut-RA − GABA-4",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600")
    ax.tick_params(colors=INK3,labelsize=5)
    ax.margins(y=.04); ax.locator_params(axis="x",nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
axes[1].legend(handles=[Line2D([],[],marker="o",color=RA,ls="",ms=3,label="Glut-RA regulon"),
                        Line2D([],[],marker="o",color=G4,ls="",ms=3,label="GABA-4 regulon")],
               frameon=False,fontsize=5,loc="lower right",borderaxespad=.1)
fig.subplots_adjust(left=.14,right=.985,top=.88,bottom=.20,wspace=.10)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/genome_wide/fig_tf13.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_tf13.png / .pdf")
