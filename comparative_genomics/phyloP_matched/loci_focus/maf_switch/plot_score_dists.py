#!/usr/bin/env python
"""Motif-score distributions for oscine vs outgroup consensus, by DAR class.

Top row: the full density, which is the point -- the four series are
superimposed through the bulk. Bottom row: the quantile difference
(oscine minus outgroup score at the same percentile) across the upper tail,
on a -log10(1-p) axis, which is where any separation actually lives.
"""
import numpy as np, matplotlib
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
Z=np.load("data/score_dists.npz"); ZE=np.load("data/score_dists_el.npz")
NB=1000; rng=np.random.default_rng(0)
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"; RA="#e31a1c"; G4="#2C6E7F"
MOT=[("full","ARE-full (12 bp)"),("jaspar","JASPAR AR (16 bp)")]
def kde(v,grid,bw=None):
    v=np.asarray(v); n=len(v)
    bw=bw or 1.06*v.std()*n**-0.2
    out=np.zeros_like(grid)
    for chunk in np.array_split(v,max(1,n//20000)):
        z=(grid[:,None]-chunk[None,:])/bw
        out+=np.exp(-.5*z*z).sum(1)
    return out/(n*bw*np.sqrt(2*np.pi))
fig,axes=plt.subplots(2,2,figsize=(5.2,3.3),
                      gridspec_kw={"hspace":.52,"wspace":.26,"height_ratios":[1,1]})
fig.patch.set_facecolor("white")
for c,(mk,ml) in enumerate(MOT):
    ax=axes[0][c]
    lo=min(Z[f"{mk}|{k}|{d}"].min() for k in ("ra","g4") for d in ("osc","out"))
    hi=max(Z[f"{mk}|{k}|{d}"].max() for k in ("ra","g4") for d in ("osc","out"))
    g=np.linspace(lo,hi,420)
    for k,col in (("ra",RA),("g4",G4)):
        for d,ls in (("osc","-"),("out",(0,(2.4,1.4)))):
            ax.plot(g,kde(Z[f"{mk}|{k}|{d}"],g),color=col,lw=.85,ls=ls,alpha=.9)
    p995=np.percentile(np.concatenate([Z[f"{mk}|{k}|{d}"] for k in ("ra","g4")
                                       for d in ("osc","out")]),99.5)
    ax.axvline(p995,color=INK3,lw=.5,ls=(0,(1.5,1.5)))
    ax.annotate("99.5th pct",xy=(p995,ax.get_ylim()[1]*.86),xytext=(-3,0),
                textcoords="offset points",ha="right",fontsize=4.3,color=INK3)
    ax.set_title(f"{'AB'[c]}   {ml}",loc="left",fontsize=7,color=INK,fontweight="600")
    ax.set_xlabel("PWM log-odds, consensus window",fontsize=6,color=INK)
    if c==0: ax.set_ylabel("density",fontsize=6,color=INK)
    ax.tick_params(colors=INK3,labelsize=5); ax.locator_params(nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    # --- tail: quantile difference
    ax=axes[1][c]
    pct=1-np.logspace(np.log10(.5),np.log10(1e-4),200)
    x=-np.log10(1-pct)
    for k,col,lab in (("ra",RA,"Glut-RA"),("g4",G4,"GABA-4")):
        o=ZE[f"{mk}|{k}|osc"]; g=ZE[f"{mk}|{k}|out"]; el=ZE[f"{mk}|{k}|elem"]
        ax.plot(x,np.quantile(o,pct)-np.quantile(g,pct),color=col,lw=1.0,label=lab,zorder=3)
        # bootstrap over elements: windows inside one element overlap heavily
        uniq=np.unique(el); idx={e:np.flatnonzero(el==e) for e in uniq}
        bs=np.empty((NB,len(pct)))
        for i in range(NB):
            sel=np.concatenate([idx[e] for e in rng.choice(uniq,len(uniq),replace=True)])
            bs[i]=np.quantile(o[sel],pct)-np.quantile(g[sel],pct)
        lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
        ax.fill_between(x,lo,hi,color=col,alpha=.15,lw=0,zorder=1)
    ax.axhline(0,color=INK,lw=.5)
    for p,lab in ((99.0,"99"),(99.5,"99.5"),(99.9,"99.9")):
        ax.axvline(-np.log10(1-p/100),color=RULE,lw=.5,ls=(0,(1.5,1.5)))
        ax.annotate(lab,xy=(-np.log10(1-p/100),ax.get_ylim()[0]),xytext=(1.5,1.5),
                    textcoords="offset points",fontsize=4.2,color=INK3)
    ax.set_xlabel("percentile  ($-\\log_{10}(1-p)$)",fontsize=6,color=INK)
    if c==0: ax.set_ylabel("oscine − outgroup\nscore at that quantile",fontsize=6,color=INK)
    ax.tick_params(colors=INK3,labelsize=5); ax.locator_params(axis="y",nbins=5)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(RULE)
    if c==1:
        ax.legend(frameon=False,fontsize=4.8,loc="upper left",borderaxespad=.1,
                  handletextpad=.4,handlelength=1.3)
axes[0][1].legend(handles=[Line2D([],[],color=INK3,lw=.9,label="oscine consensus"),
                           Line2D([],[],color=INK3,lw=.9,ls=(0,(2.4,1.4)),label="outgroup consensus"),
                           Line2D([],[],color=RA,lw=.9,label="Glut-RA DARs"),
                           Line2D([],[],color=G4,lw=.9,label="GABA-4 DARs")],
                  frameon=False,fontsize=4.6,loc="upper left",borderaxespad=.1,
                  handletextpad=.4,handlelength=1.3,labelspacing=.3)
fig.subplots_adjust(left=.105,right=.985,top=.925,bottom=.115)
for e in ("png","pdf"):
    fig.savefig(f"../results/mafb_domain/fig_score_dists.{e}",dpi=220,facecolor="white",
                bbox_inches="tight",pad_inches=.015)
print("  wrote fig_score_dists (with element-bootstrap bands)")
for mk,ml in MOT:
    print(f"\n  {ml}: oscine minus outgroup quantile")
    for p in (50,90,99,99.5,99.9,99.99):
        r=np.quantile(Z[f"{mk}|ra|osc"],p/100)-np.quantile(Z[f"{mk}|ra|out"],p/100)
        g=np.quantile(Z[f"{mk}|g4|osc"],p/100)-np.quantile(Z[f"{mk}|g4|out"],p/100)
        print(f"    P{p:<6} Glut-RA {r:+6.2f}   GABA-4 {g:+6.2f}   diff {r-g:+6.2f}")
