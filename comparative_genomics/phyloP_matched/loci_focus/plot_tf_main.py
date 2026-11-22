#!/usr/bin/env python
"""Motif minus the rest of its own enhancer at MAFB, ranked on one combined ordering."""
import numpy as np, sys, json, matplotlib
import sys as _s; _s.path.insert(0, "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"); from span import result  # results/<span>/
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
rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
def ci(a,n=3000):
    bs=[np.mean(rng.choice(a,len(a))) for _ in range(n)]
    return a.mean(), *np.percentile(bs,[2.5,97.5])
Q={f'{t}|{l or "gw"}':{r["tf"]:r for r in json.load(open(result(f"tf_tab_{t}_{l or 'gw'}.json", l)))}
   for t in ("subtree222","wholetree") for l in (None,"MAFB","domain")}
def panel(ax,tree,loc,mn,ttl,order=None):
    D=diffs(tree,loc); V={}
    for m in REG:
        a=D.get((m,"RA_over_Arco"),np.array([])); b=D.get((m,"GABA4_over_GABA1"),np.array([]))
        if len(a)>=mn and len(b)>=mn: V[m]=(ci(a),ci(b))
    if order is None: order=sorted(V,key=lambda m:-V[m][0][0])
    q=Q[f'{tree}|{loc or "gw"}']
    for i,m in enumerate(order):
        if m not in V: continue
        (r,rl,rh),(g,gl,gh)=V[m]; y=len(order)-1-i
        qr=q.get(m,{}).get("q_pra",1); qg=q.get(m,{}).get("q_pg4",1)
        ax.plot([min(r,g),max(r,g)],[y]*2,color=RULE,lw=2.0,solid_capstyle="round",zorder=1)
        for v,lo,hi,col,off,qq in ((g,gl,gh,G4,-.19,qg),(r,rl,rh,RA,+.19,qr)):
            ax.plot([lo,hi],[y+off]*2,color=col,lw=.65,alpha=.7,solid_capstyle="round",zorder=2)
            ax.scatter([v],[y+off],s=13 if m=="AR" else 8.5,color=col,zorder=4,linewidths=0,
                       marker="D" if m=="AR" else "o",alpha=1 if qq<.05 else .45)
    if "AR" in order:
        ya=len(order)-1-order.index("AR")
        ax.axhspan(ya-.45,ya+.45,color=RA,alpha=.055,zorder=0,lw=0)
    ax.axvline(0,color=INK,lw=.5,zorder=0)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order[::-1],fontsize=5,color=INK)
    for t,m in zip(ax.get_yticklabels(),order[::-1]):
        if m=="AR": t.set_fontweight("700"); t.set_color(RA)
    ax.set_xlabel("motif − rest of enhancer  (log$_2$ rate)",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600")
    ax.tick_params(colors=INK3,labelsize=5)
    ax.margins(y=.035); ax.locator_params(axis="x",nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    return order
LEG=[Line2D([],[],marker="o",color=RA,ls="",ms=3,label="Glut-RA"),
     Line2D([],[],marker="o",color=G4,ls="",ms=3,label="GABA-4")]
fig,ax=plt.subplots(1,2,figsize=(4.9,1.95),sharey=True); fig.patch.set_facecolor("white")
o=panel(ax[0],"subtree222","MAFB",8,"A   oscine subtree, MAFB ±1.5 Mb")
panel(ax[1],"wholetree","MAFB",8,"B   whole tree, MAFB ±1.5 Mb",o)
lg=ax[0].legend(handles=LEG,fontsize=5,loc="lower left",borderaxespad=.15,
                handletextpad=.3,labelspacing=.22,frameon=True,framealpha=.88,
                facecolor="white",edgecolor="none")
lg.get_frame().set_linewidth(0)
fig.subplots_adjust(left=.13,right=.985,top=.88,bottom=.21,wspace=.08)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/mafb_1p5mb/fig_tf_mafb_main.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("wrote fig_tf_mafb_main.png / .pdf")
fig,ax=plt.subplots(1,2,figsize=(4.9,1.95),sharey=True); fig.patch.set_facecolor("white")
od=panel(ax[0],"subtree222","domain",8,"A   oscine subtree, DHX35\u2013TOP1 domain")
panel(ax[1],"wholetree","domain",8,"B   whole tree, DHX35\u2013TOP1 domain",od)
lg=ax[0].legend(handles=LEG,fontsize=5,loc="lower left",borderaxespad=.15,
                handletextpad=.3,labelspacing=.22,frameon=True,framealpha=.88,
                facecolor="white",edgecolor="none")
lg.get_frame().set_linewidth(0)
fig.subplots_adjust(left=.13,right=.985,top=.88,bottom=.21,wspace=.08)
for e in ("png","pdf"):
    fig.savefig(f"{result('fig_tf_domain_main','domain')}.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("wrote fig_tf_domain_main.png / .pdf")
for m in o:
    d=Q["subtree222|MAFB"][m]
    print("  %-7s Glut-RA %+.3f (q %.3f)   GABA-4 %+.3f (q %.3f)"%(
        m,d["ra"],d["q_pra"],d["g4"],d["q_pg4"]))
