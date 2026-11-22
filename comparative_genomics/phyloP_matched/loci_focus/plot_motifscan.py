#!/usr/bin/env python
"""Where AR sits among 182 motif archetypes, genome-wide."""
import json, numpy as np, matplotlib
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
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#1F3A5F"; G4="#33a02c"; GREY="#B9BFC7"
rows=json.load(open(f"{F}/results/genome_wide/allmotif_subtree222.json"))
x=np.array([r["d_both"] for r in rows]); y=np.array([r["d_diff"] for r in rows])
ra=np.array([r["d_ra"] for r in rows]); g4=np.array([r["d_g4"] for r in rows])
ar=[i for i,r in enumerate(rows) if "MA0007.3" in r["motif"]][0]
both=(ra>0)&(g4>0)
fig,axes=plt.subplots(1,2,figsize=(3.40,1.40),gridspec_kw={"width_ratios":[1.1,1]})
fig.patch.set_facecolor("white")
ax=axes[0]
ax.scatter(x[~both],y[~both],s=3.9,color=GREY,alpha=.65,linewidths=0,label="accelerated in one or both classes")
ax.scatter(x[both],y[both],s=4.4,color=INK3,alpha=.7,linewidths=0,label="constrained in both classes")
ax.scatter([x[ar]],[y[ar]],s=23.1,color=RA,zorder=5,linewidths=0)
ax.annotate("AR  (MA0007.3)",(x[ar],y[ar]),textcoords="offset points",xytext=(12,4),
            fontsize=5,color=RA,fontweight="600")
ax.axhline(0,color=RULE,lw=0.50); ax.axvline(0,color=RULE,lw=0.50)
ax.set_xlabel("motif − enhancer (pooled)",fontsize=6,color=INK)
ax.set_ylabel("Glut-RA − GABA-4",fontsize=6,color=INK)
ax.set_title("A   182 motif archetypes, genome-wide",loc="left",fontsize=7,
             color=INK,fontweight="600",pad=10)
ax.legend(frameon=False,fontsize=5,loc="lower left")
ax=axes[1]
o=np.argsort(-y); rank=np.argsort(o)
ax.scatter(rank[~both]+1,y[~both],s=3.9,color=GREY,alpha=.65,linewidths=0)
ax.scatter(rank[both]+1,y[both],s=4.4,color=INK3,alpha=.7,linewidths=0)
ax.scatter([rank[ar]+1],[y[ar]],s=23.1,color=RA,zorder=5,linewidths=0)
ax.annotate(f"AR — rank {rank[ar]+1}/182",(rank[ar]+1,y[ar]),textcoords="offset points",
            xytext=(16,2),fontsize=5,color=RA,fontweight="600")
lo,hi=np.percentile(y,[2.5,97.5])
ax.axhspan(lo,hi,color=INK3,alpha=.08,linewidth=0,zorder=0)
ax.axhline(0,color=RULE,lw=0.50)
ax.set_xlabel("rank",fontsize=6,color=INK)
ax.set_ylabel("Glut-RA − GABA-4",fontsize=6,color=INK)
ax.set_title("B   the cell-type asymmetry, ranked",loc="left",fontsize=7,
             color=INK,fontweight="600",pad=10)
for a in axes:
    a.margins(x=.03,y=.05)
    a.locator_params(nbins=5)
    for sp in ("top","right"): a.spines[sp].set_visible(False)
    for sp in ("left","bottom"): a.spines[sp].set_color(RULE)
    a.tick_params(colors=INK3,labelsize=5)
fig.subplots_adjust(left=.075,right=.985,top=.90,bottom=.235,wspace=.16)
fig.savefig(f"{F}/results/genome_wide/fig_motifscan.png",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
fig.savefig(f"{F}/results/genome_wide/fig_motifscan.pdf",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
print("  wrote fig_motifscan.png")
