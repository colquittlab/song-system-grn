#!/usr/bin/env python
"""SCENIC+ MAFB regulators including AR: cell-type interaction, two views."""
import numpy as np, sys, json, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4,"xtick.major.width":.4,"ytick.major.width":.4,
    "xtick.major.size":1.5,"ytick.major.size":1.5,"xtick.major.pad":1,
    "ytick.major.pad":1,"axes.labelpad":1,"axes.titlepad":2})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
sys.path.insert(0,F)
from tf_core import REG, diffs
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
def ci(a,b,n=3000):
    bs=[np.mean(rng.choice(a,len(a)))-np.mean(rng.choice(b,len(b))) for _ in range(n)]
    return a.mean()-b.mean(), *np.percentile(bs,[2.5,97.5])
Q={f"{t}|{l or 'gw'}":{r["tf"]:r["q_p"] for r in json.load(open(result(f"tf_tab_{t}_{l or 'gw'}.json", l)))}
   for t in ("subtree222","wholetree") for l in (None,"MAFB","domain")}
def panel(ax,tree,loc,mn,ttl,order=None):
    D=diffs(tree,loc); vals={}
    for m in REG:
        a=D.get((m,"RA_over_Arco"),np.array([])); b=D.get((m,"GABA4_over_GABA1"),np.array([]))
        if len(a)>=mn and len(b)>=mn: vals[m]=ci(a,b)
    if order is None: order=sorted(vals,key=lambda m:-vals[m][0])
    q=Q[f"{tree}|{loc or 'gw'}"]
    for i,m in enumerate(order):
        if m not in vals: continue
        d,lo,hi=vals[m]; y=len(order)-1-i
        col=RA if REG[m]=="Glut-RA" else G4
        sig=q.get(m,1)<.05
        ax.plot([lo,hi],[y]*2,color=col,lw=.7,alpha=.8 if sig else .45,solid_capstyle="round")
        ax.scatter([d],[y],s=15 if m=="AR" else 9,color=col,zorder=3,linewidths=0,
                   marker="D" if m=="AR" else "o",alpha=1 if sig else .55)
        if sig: ax.text(hi+.004,y,"*",va="center",ha="left",fontsize=6,color=col,
                        fontweight="600",transform=ax.transData)
    ax.axvline(0,color=INK,lw=.5,zorder=1)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1],fontsize=5,color=INK)
    for t,m in zip(ax.get_yticklabels(),order[::-1]):
        if m=="AR": t.set_fontweight("600")
    ax.set_xlabel("Glut-RA − GABA-4",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600")
    ax.tick_params(colors=INK3,labelsize=5)
    ax.margins(y=.04); ax.locator_params(axis="x",nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    return order
LEG=[Line2D([],[],marker="o",color=RA,ls="",ms=3,label="Glut-RA regulon"),
     Line2D([],[],marker="o",color=G4,ls="",ms=3,label="GABA-4 regulon"),
     Line2D([],[],marker="D",color=RA,ls="",ms=3.5,label="AR")]
# figure 1: oscine vs whole tree, genome-wide
fig,ax=plt.subplots(1,2,figsize=(4.9,1.85),sharey=True); fig.patch.set_facecolor("white")
o=panel(ax[0],"subtree222",None,40,"A   oscine subtree, genome-wide")
panel(ax[1],"wholetree",None,40,"B   whole tree, genome-wide",o)
ax[1].legend(handles=LEG,frameon=False,fontsize=5,loc="lower right",borderaxespad=.1)
fig.subplots_adjust(left=.13,right=.985,top=.88,bottom=.21,wspace=.10)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/genome_wide/fig_tf13.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
# figure 2: genome-wide vs MAFB, oscine
fig,ax=plt.subplots(1,2,figsize=(4.9,1.85),sharey=True); fig.patch.set_facecolor("white")
o=panel(ax[0],"subtree222",None,40,"A   genome-wide")
panel(ax[1],"subtree222","MAFB",8,"B   at MAFB, gene ±1.5 Mb",o)
ax[1].legend(handles=LEG,frameon=False,fontsize=5,loc="lower right",borderaxespad=.1)
fig.subplots_adjust(left=.13,right=.985,top=.88,bottom=.21,wspace=.10)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_tf13_mafb.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
# figure 3: genome-wide vs the DHX35-TOP1 domain, oscine
fig,ax=plt.subplots(1,2,figsize=(4.9,1.85),sharey=True); fig.patch.set_facecolor("white")
o=panel(ax[0],"subtree222",None,40,"A   genome-wide")
panel(ax[1],"subtree222","domain",8,"B   DHX35\u2013TOP1 domain",o)
ax[1].legend(handles=LEG,frameon=False,fontsize=5,loc="lower right",borderaxespad=.1)
fig.subplots_adjust(left=.13,right=.985,top=.88,bottom=.21,wspace=.10)
for e in ("png","pdf"):
    fig.savefig(f"{result('fig_tf13_domain','compare_gw_domain')}.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  rebuilt fig_tf13, fig_tf13_mafb, fig_tf13_domain  (* marks FDR 5%)")
