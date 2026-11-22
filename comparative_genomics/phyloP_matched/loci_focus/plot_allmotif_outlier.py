#!/usr/bin/env python
"""eRegulon TFs against the 182-archetype whole-tree / oscine DoD trend."""
import numpy as np, json, sys, matplotlib
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
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; GREY="#AEB6C2"; RA="#e31a1c"; G4="#2C6E7F"

def place(x,y,labs,xr,yr):
    px=(np.array(x)-xr[0])/(xr[1]-xr[0]); py=(np.array(y)-yr[0])/(yr[1]-yr[0])
    W=lambda t:len(t)*0.0265+0.012; H=0.050
    cand=[(1,0,"left"),(-1,0,"right"),(1,.9,"left"),(-1,.9,"right"),
          (1,-.9,"left"),(-1,-.9,"right"),(0,1.5,"center"),(0,-1.5,"center"),
          (1,1.8,"left"),(-1,1.8,"right"),(1,-1.8,"left"),(-1,-1.8,"right")]
    boxes=[]; res=[None]*len(labs)
    for i in sorted(range(len(labs)),key=lambda i:-abs(py[i]-.5)):
        w=W(labs[i]); best=None
        for dx,dy,ha in cand:
            cx=px[i]+dx*0.016; cy=py[i]+dy*0.036
            x0=cx if ha=="left" else (cx-w if ha=="right" else cx-w/2)
            box=(x0,cy-H/2,x0+w,cy+H/2); pen=.004*abs(dy)
            for b in boxes:
                pen+=max(0,min(box[2],b[2])-max(box[0],b[0]))*\
                     max(0,min(box[3],b[3])-max(box[1],b[1]))*150
            for j in range(len(labs)):
                pen+=max(0,min(box[2],px[j]+.012)-max(box[0],px[j]-.012))*\
                     max(0,min(box[3],py[j]+.022)-max(box[1],py[j]-.022))*110
            pen+=(max(0,-box[0])+max(0,box[2]-1)+max(0,-box[1])+max(0,box[3]-1))*10
            if best is None or pen<best[0]: best=(pen,cx,cy,ha,box)
        boxes.append(best[4]); res[i]=(best[1],best[2],best[3])
    return [(c[0]*(xr[1]-xr[0])+xr[0], c[1]*(yr[1]-yr[0])+yr[0], c[2]) for c in res]

def panel(ax,tagl,ttl,suffix=""):
    D=json.load(open(result(f"allmotif_outlier_{tagl}{suffix}.json", tagl)))
    fit=D["fit"]; A=D["arch"]; E=D["ereg"]
    ax_=np.array([r["x"] for r in A]); ay=np.array([r["y"] for r in A])
    ex=np.array([r["x"] for r in E]); ey=np.array([r["y"] for r in E])
    allx=np.concatenate([ax_,ex]); ally=np.concatenate([ay,ey])
    pad=.07
    xr=(allx.min()-pad*np.ptp(allx), allx.max()+pad*np.ptp(allx))
    yr=(ally.min()-pad*np.ptp(ally), ally.max()+pad*np.ptp(ally))
    ax.axhline(0,color=RULE,lw=.4,zorder=0); ax.axvline(0,color=RULE,lw=.4,zorder=0)
    xs=np.linspace(*xr,50); ln=fit["a"]+fit["b"]*xs
    ax.fill_between(xs,ln-2*fit["sig"],ln+2*fit["sig"],color=GREY,alpha=.18,lw=0,zorder=1)
    ax.plot(xs,ln,color=INK3,lw=.6,zorder=2)
    ax.scatter(ax_,ay,s=3.4,color=GREY,alpha=.6,linewidths=0,zorder=3)
    for r,xi,yi in zip(E,ex,ey):
        col=RA if r["reg"]=="Glut-RA" else G4
        ar=r["tf"].rstrip("*")=="AR"
        ax.scatter([xi],[yi],s=17 if ar else 11,color=col,zorder=5,linewidths=0,
                   marker="D" if ar else "o")
    for r,xi,yi,(tx,ty,ha) in zip(E,ex,ey,place(ex,ey,[r["tf"] for r in E],xr,yr)):
        col=RA if r["reg"]=="Glut-RA" else G4
        ax.plot([xi,tx],[yi,ty],color=col,lw=.35,alpha=.4,zorder=4)
        ax.text(tx,ty,r["tf"],fontsize=4.4,color=col,va="center",ha=ha,zorder=6,
                fontweight="700" if r["tf"].rstrip("*")=="AR" else "400")
    ax.set_xlim(*xr); ax.set_ylim(*yr)
    ax.set_xlabel("whole-tree DoD",fontsize=6,color=INK)
    ax.set_title(ttl,loc="left",fontsize=7,color=INK,fontweight="600")
    ax.text(.025,.965,"r = %+.2f\n$\\tau$ = %.3f\n%d archetypes"%(fit["r"],fit["tau"],len(A)),
            transform=ax.transAxes,fontsize=4.6,color=INK3,va="top",ha="left",linespacing=1.4)
    ax.tick_params(colors=INK3,labelsize=5); ax.locator_params(nbins=4)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)

LEG=[Line2D([],[],marker="o",color=GREY,ls="",ms=2.3,label="motif archetype"),
     Line2D([],[],marker="o",color=RA,ls="",ms=3,label="Glut-RA regulon"),
     Line2D([],[],marker="o",color=G4,ls="",ms=3,label="GABA-4 regulon"),
     Line2D([],[],marker="D",color=INK3,ls="",ms=3,label="AR = SCENIC+ matrix (= NR3C2), AR* = HOMER ARE-full")]


def figure(panels, out, width):
    fig,ax=plt.subplots(1,len(panels),figsize=(width,2.15)); fig.patch.set_facecolor("white")
    for a,args in zip(ax,panels): panel(a,*args)
    ax[0].set_ylabel("oscine DoD",fontsize=6,color=INK)
    lg=ax[0].legend(handles=LEG,fontsize=4.5,loc="lower right",borderaxespad=.12,
                    handletextpad=.3,labelspacing=.2,frameon=True,framealpha=.9,
                    facecolor="white",edgecolor="none")
    lg.get_frame().set_linewidth(0)
    fig.subplots_adjust(left=.105,right=.99,top=.885,bottom=.145,wspace=.26)
    for e in ("png","pdf"):
        fig.savefig(f"{out}.{e}",dpi=200,facecolor="white",bbox_inches="tight",pad_inches=.01)
    print(f"wrote {out}.png / .pdf")


# genome-wide vs the MAFB gene +/-1.5 Mb window
figure([("gw","A   genome-wide"),("MAFB","B   at MAFB, gene \u00b11.5 Mb")],
       f"{F}/results/compare_gw_mafb1p5mb/fig_allmotif_outlier",4.9)
# genome-wide vs the DHX35-TOP1 domain, at the standard and the relaxed element minimum
figure([("gw","A   genome-wide"),
        ("domain","B   DHX35\u2013TOP1 domain, \u22658 per class"),
        ("domain","C   domain, \u22654 per class (sensitivity)","_min4")],
       result("fig_allmotif_outlier_domain","compare_gw_domain"),7.2)
