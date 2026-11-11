#!/usr/bin/env python
"""A targeted alignment column mapped onto the species tree, leaves and
reconstructed internal nodes alike.

States come straight out of the HAL: cactus stores a reconstructed sequence for
every ancestral node, so the base at an internal node is read, not inferred
here. Clades that are monomorphic and contain no oscine are collapsed to a
wedge so the substitution is the thing you see.

usage: plot_maf_tree.py <col.json> <tree.nw> <out_stem> "<title>"
"""
import json, re, sys, collections, numpy as np, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
    "pdf.fonttype":42,"ps.fonttype":42,
    "font.size":5,"axes.labelsize":6,"axes.titlesize":7,
    "xtick.labelsize":5,"ytick.labelsize":5,"legend.fontsize":5,
    "axes.linewidth":.4})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
INK="#171B22"; INK3="#6E7887"; RULE="#DDE1E7"
BASE={"A":"#2E8B57","C":"#e31a1c","G":"#B8860B","T":"#2C6E7F","-":"#C9CED6","N":"#C9CED6"}
col=json.load(open(sys.argv[1] if len(sys.argv)>1 else f"{F}/results/mafb_1p5mb/maf_target_col.json"))
tree=open(sys.argv[2] if len(sys.argv)>2 else "/tmp/tree.nw").read().strip().rstrip(";")
STEM=sys.argv[3] if len(sys.argv)>3 else f"{F}/results/mafb_1p5mb/fig_maf_tree"
TITLE=sys.argv[4] if len(sys.argv)>4 else ""


def parse(s):
    pos=[0]
    def node():
        ch=[]
        if s[pos[0]]=="(":
            pos[0]+=1
            while True:
                ch.append(node())
                if pos[0]<len(s) and s[pos[0]]==",": pos[0]+=1
                else: break
            pos[0]+=1
        m=re.match(r"[^,()\[\]:;]*",s[pos[0]:]); name=m.group(0); pos[0]+=len(name)
        bl=0.0
        if pos[0]<len(s) and s[pos[0]]==":":
            m=re.match(r":([-\d.eE+]*)",s[pos[0]:]); bl=float(m.group(1) or 0); pos[0]+=len(m.group(0))
        return dict(name=name,bl=bl,ch=ch)
    return node()
root=parse(tree)
def leaves(n): return [n] if not n["ch"] else [x for c in n["ch"] for x in leaves(c)]
# oscines = smallest clade holding the songbird anchors but not Acanthisitta,
# taken from THIS tree so the highlight can never disagree with what is drawn
ANCHOR={"Taeniopygia_guttata","Corvus_hawaiiensis","Lonchura_striata"}
OSC=set(); _best=[None]
def _osc(n):
    L={x["name"] for x in leaves(n)}
    if ANCHOR<=L and "Acanthisitta_chloris" not in L:
        if _best[0] is None or len(L)<len(_best[0]): _best[0]=L
    for c in n["ch"]: _osc(c)
_osc(root); OSC=_best[0] or set()
def states(n): return {col.get(x["name"]) for x in leaves(n)}|{col.get(n["name"])}
def has_osc(n): return any(x["name"] in OSC for x in leaves(n))
# collapse maximal monomorphic, oscine-free clades
def mark(n):
    st=states(n)-{None}
    n["collapse"]= (not has_osc(n)) and len(st)==1 and len(leaves(n))>1
    if not n["collapse"]:
        for c in n["ch"]: mark(c)
mark(root)
def rotate(n):
    """put the oscine-bearing child last at every node, so Passeri sits at the
    bottom of the drawing (layout assigns y in visit order, y inverted)"""
    if not n["ch"]: return
    for c in n["ch"]: rotate(c)
    bear=[c for c in n["ch"] if any(x["name"] in OSC for x in leaves(c))]
    if len(bear)==1:
        n["ch"].remove(bear[0]); n["ch"].append(bear[0])
rotate(root)
order=[]
def setx(n,x):                      # every node gets x, collapsed ones included
    n["x"]=x+n["bl"]
    for c in n["ch"]: setx(c,n["x"])
setx(root,0.0)
def layout(n):                      # only drawn tips get a y slot
    if n["collapse"] or not n["ch"]:
        n["y"]=len(order); order.append(n); return
    for c in n["ch"]: layout(c)
    n["y"]=float(np.mean([c["y"] for c in n["ch"]]))
layout(root)
N=len(order)
fig,ax=plt.subplots(figsize=(4.15,max(2.0,N*0.105+0.55)))
fig.patch.set_facecolor("white")
def draw(n):
    for c in n["ch"] if not n["collapse"] else []:
        sp,sc=col.get(n["name"]),col.get(c["name"])
        chg = sp and sc and sp!=sc
        ax.plot([n["x"],n["x"]],[n["y"],c["y"]],color=RULE,lw=.55,zorder=1,solid_capstyle="round")
        ax.plot([n["x"],c["x"]],[c["y"],c["y"]],
                color="#C7521B" if chg else RULE,lw=1.5 if chg else .55,zorder=2,
                solid_capstyle="round")
        if chg:
            ax.annotate(f"{sp}→{sc}",xy=((n["x"]+c["x"])/2,c["y"]),
                        xytext=(0,3.1),textcoords="offset points",ha="center",
                        va="bottom",fontsize=4.4,color="#8A3A12",fontweight="700")
        draw(c)
    b=col.get(n["name"])
    if n["collapse"]:
        L=leaves(n); y0,y1=n["y"]-.36,n["y"]+.36
        xm=max(x["x"] for x in L)
        ax.fill([n["x"],xm,xm],[n["y"],y0,y1],color=BASE.get(b,"#C9CED6"),alpha=.30,
                lw=0,zorder=1)
        ax.text(xm+.012,n["y"],f"  {len(L)} taxa",va="center",ha="left",
                fontsize=4.2,color=INK3,style="italic")
    if b:
        ax.scatter([n["x"]],[n["y"]],s=11 if not n["ch"] else 7.5,
                   color=BASE.get(b,"#C9CED6"),zorder=4,linewidths=0,
                   marker="o" if n["ch"] else "s")
    if not n["ch"]:
        nm=n["name"].replace("_"," ")
        ax.text(n["x"]+.012,n["y"],nm,va="center",ha="left",fontsize=4.2,
                color="#8A1416" if n["name"] in OSC else INK,
                fontweight="700" if n["name"] in OSC else "400",style="italic")
draw(root)
ax.set_ylim(-.9,N-.1); ax.invert_yaxis()
xs=[n["x"] for n in order]; ax.set_xlim(-.02,max(xs)*1.37)
ax.axis("off")
if TITLE: ax.set_title(TITLE,loc="left",fontsize=6.5,color=INK,fontweight="600",pad=3)
seen=[b for b in "ACGT" if b in set(col.values())]
ax.legend(handles=[Line2D([],[],marker="s",color=BASE[b],ls="",ms=3,label=b) for b in seen]
          +[Line2D([],[],color="#C7521B",lw=1.4,label="substitution"),
            Line2D([],[],color="none",label="bold = oscine (Passeri)")],
          frameon=False,fontsize=4.6,loc="lower left",ncol=len(seen)+2,
          bbox_to_anchor=(0,-.045),handletextpad=.3,columnspacing=.8)
print(f"  oscines highlighted: {len(OSC)}")
fig.subplots_adjust(left=.01,right=.99,top=.955,bottom=.045)
for e in ("png","pdf"):
    fig.savefig(f"{STEM}.{e}",dpi=240,facecolor="white",bbox_inches="tight",pad_inches=.015)
print(f"  wrote {STEM}  ({N} tips drawn, {sum(1 for n in order if n['collapse'])} collapsed clades)")
