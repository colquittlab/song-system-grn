#!/usr/bin/env python
"""Oscine vs whole-tree difference-of-differences for the 13 SCENIC+ MAFB regulators."""
import numpy as np, json, matplotlib
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
from scipy import stats
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"

def place(x,y,labs,xr,yr,ar):
    """greedy label placement in axis-fraction space; returns (lx,ly,ha) per point"""
    px=(np.array(x)-xr[0])/(xr[1]-xr[0]); py=(np.array(y)-yr[0])/(yr[1]-yr[0])
    W=lambda t:len(t)*0.0255/ar + 0.012      # approx text width, axis fraction
    H=0.052                                   # approx text height
    cand=[(1,0,"left"),(-1,0,"right"),(1,.9,"left"),(-1,.9,"right"),
          (1,-.9,"left"),(-1,-.9,"right"),(0,1.5,"center"),(0,-1.5,"center")]
    boxes=[]; out=[]
    order=sorted(range(len(labs)),key=lambda i:-abs(py[i]-.5))
    res=[None]*len(labs)
    for i in order:
        w=W(labs[i]); best=None
        for dx,dy,ha in cand:
            cx=px[i]+dx*(0.014+ (w/2 if ha=="center" else 0)); cy=py[i]+dy*0.034
            x0=cx if ha=="left" else (cx-w if ha=="right" else cx-w/2)
            box=(x0,cy-H/2,x0+w,cy+H/2)
            pen=0.0
            for b in boxes:
                ox=max(0,min(box[2],b[2])-max(box[0],b[0]))
                oy=max(0,min(box[3],b[3])-max(box[1],b[1]))
                pen+=ox*oy*120
            for j in range(len(labs)):          # keep off the markers
                ox=max(0,min(box[2],px[j]+.012)-max(box[0],px[j]-.012))
                oy=max(0,min(box[3],py[j]+.022)-max(box[1],py[j]-.022))
                pen+=ox*oy*90
            pen+=max(0,-box[0])*8+max(0,box[2]-1)*8
            pen+=max(0,-box[1])*8+max(0,box[3]-1)*8
            pen+=0.004*abs(dy)
            if best is None or pen<best[0]: best=(pen,cx,cy,ha,box)
        boxes.append(best[4]); res[i]=(best[1],best[2],best[3])
    return [(c[0]*(xr[1]-xr[0])+xr[0], c[1]*(yr[1]-yr[0])+yr[0], c[2]) for c in res]

def panel(ax,loc,ttl):
    tag=loc or "gw"
    O={r["tf"]:r for r in json.load(open(result(f"tf_tab_subtree222_{tag}.json", tag)))}
    W={r["tf"]:r for r in json.load(open(result(f"tf_tab_wholetree_{tag}.json", tag)))}
    tfs=[t for t in O if t in W]
    x=np.array([W[t]["dod"] for t in tfs]); y=np.array([O[t]["dod"] for t in tfs])
    r,p=stats.pearsonr(x,y)
    pad=.14
    xr=(x.min()-pad*np.ptp(x), x.max()+pad*np.ptp(x))
    yr=(y.min()-pad*np.ptp(y), y.max()+pad*np.ptp(y))
    ax.axhline(0,color=RULE,lw=.5,zorder=0); ax.axvline(0,color=RULE,lw=.5,zorder=0)
    b,a=np.polyfit(x,y,1)
    xs=np.array(xr); ax.plot(xs,a+b*xs,color=INK3,lw=.5,ls=(0,(3,2)),zorder=1,alpha=.7)
    for t,xi,yi in zip(tfs,x,y):
        col=RA if O[t]["reg"]=="Glut-RA" else G4
        sig=O[t]["q_p"]<.05
        ax.scatter([xi],[yi],s=18 if t=="AR" else 11,color=col,zorder=4,linewidths=0,
                   marker="D" if t=="AR" else "o",alpha=1 if sig else .5)
    ar=(xr[1]-xr[0])/(yr[1]-yr[0])*0.90/1.0
    for t,xi,yi,(tx,ty,ha) in zip(tfs,x,y,place(x,y,tfs,xr,yr,1.0)):
        col=RA if O[t]["reg"]=="Glut-RA" else G4
        ax.plot([xi,tx],[yi,ty],color=col,lw=.35,alpha=.35,zorder=2)
        ax.text(tx,ty,t,fontsize=4.4,color=col,va="center",ha=ha,zorder=5,
                fontweight="700" if t=="AR" else "400")
    ax.set_xlim(*xr); ax.set_ylim(*yr)
    ax.set_xlabel("whole-tree DoD",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600")
    rx,rha=(.97,"right") if loc=="MAFB" else (.03,"left")
    ax.text(rx,.955,"r = %+.2f, p = %.2f"%(r,p),transform=ax.transAxes,
            fontsize=4.8,color=INK3,va="top",ha=rha)
    ax.tick_params(colors=INK3,labelsize=5); ax.locator_params(nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    return r,p,len(tfs)

fig,ax=plt.subplots(1,2,figsize=(4.4,2.0)); fig.patch.set_facecolor("white")
s1=panel(ax[0],None,"A   genome-wide")
s2=panel(ax[1],"MAFB","B   at MAFB, gene ±1.5 Mb")
ax[0].set_ylabel("oscine DoD",fontsize=6,color=INK)
LEG=[Line2D([],[],marker="o",color=RA,ls="",ms=3,label="Glut-RA regulon"),
     Line2D([],[],marker="o",color=G4,ls="",ms=3,label="GABA-4 regulon"),
     Line2D([],[],marker="D",color=RA,ls="",ms=3.3,label="AR")]
lg=ax[1].legend(handles=LEG,fontsize=4.6,loc="upper left",
                bbox_to_anchor=(.015,.88),borderaxespad=0,
                handletextpad=.3,labelspacing=.22,frameon=True,framealpha=.88,
                facecolor="white",edgecolor="none")
lg.get_frame().set_linewidth(0)
fig.subplots_adjust(left=.115,right=.99,top=.88,bottom=.155,wspace=.30)
for e in ("png","pdf"):
    fig.savefig(f"{F}/results/compare_gw_mafb1p5mb/fig_tf_dod_scatter.{e}",dpi=200,facecolor="white",
                bbox_inches="tight",pad_inches=.01)
print("wrote fig_tf_dod_scatter.png / .pdf")
print("  genome-wide  r=%+.3f p=%.3f  n=%d"%s1)
print("  MAFB         r=%+.3f p=%.3f  n=%d"%s2)
fig,ax=plt.subplots(1,2,figsize=(4.4,2.0)); fig.patch.set_facecolor("white")
s1=panel(ax[0],None,"A   genome-wide")
s3=panel(ax[1],"domain","B   DHX35\u2013TOP1 domain")
ax[0].set_ylabel("oscine DoD",fontsize=6,color=INK)
lg=ax[1].legend(handles=LEG,fontsize=4.6,loc="upper left",bbox_to_anchor=(.015,.88),borderaxespad=0,
                handletextpad=.3,labelspacing=.22,frameon=True,framealpha=.88,facecolor="white",edgecolor="none")
lg.get_frame().set_linewidth(0)
fig.subplots_adjust(left=.115,right=.99,top=.88,bottom=.155,wspace=.30)
for e in ("png","pdf"):
    fig.savefig(f"{result('fig_tf_dod_scatter_domain','compare_gw_domain')}.{e}",dpi=200,facecolor="white",
                bbox_inches="tight",pad_inches=.01)
print("  domain       r=%+.3f p=%.3f  n=%d"%s3)
