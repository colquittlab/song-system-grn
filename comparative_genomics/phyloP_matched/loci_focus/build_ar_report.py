#!/usr/bin/env python
"""Generate the AR-motif report from the consolidated JSON, so every number is traceable."""
import json, base64, os, collections
import numpy as np
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus"
M=f"{F}/maf_switch"; OUT=f"{F}/../report/ar_motif_report.html"
A=json.load(open(f"{F}/results/compare_mafb_spans/ar_report_part1.json")); B=json.load(open(f"{F}/results/compare_mafb_spans/ar_report_part23.json"))
def num(x,d=3,sign=True):
    s=f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("-","−")
def pv(p):
    c="sig" if p<.05 else ("near" if p<.1 else "ns")
    t="&lt;0.001" if p<1e-3 else (f"{p:.3f}" if p<.01 else f"{p:.2f}")
    return f'<td class="n {c}">{t}</td>'
def img(path,alt,cap):
    if not os.path.exists(path): return ""
    b=base64.b64encode(open(path,"rb").read()).decode()
    return (f'<figure><img alt="{alt}" src="data:image/png;base64,{b}">'
            f'<figcaption>{cap}</figcaption></figure>')
MN={"full":"ARE-full","half":"ARE-half","jaspar":"JASPAR AR"}
MW={"full":"12 bp, HOMER de novo","half":"8 bp, HOMER de novo","jaspar":"16 bp, MA0007.3"}

# ---------- part 1 tables ----------
def t_gw(tree):
    r=[]
    for cols,cl in (("all","all columns"),("hi","high-IC only"),("lo","low-IC only")):
        for mk in ("full","half","jaspar"):
            k=f"{tree}|{cols}|{mk}"
            if k not in A["gw"]: continue
            d=A["gw"][k]
            r.append(f'<tr><td>{cl}</td><td>{MN[mk]}</td><td class="n ra">{num(d["ra"])}</td>'
                     f'{pv(d["p_ra"])}<td class="n g4">{num(d["g4"])}</td>{pv(d["p_g4"])}'
                     f'<td class="n hl">{num(d["dod"])}</td>{pv(d["p_dod"])}'
                     f'<td class="n sm">{d["n_ra"]:,} / {d["n_g4"]:,}</td></tr>')
    return "\n".join(r)
def t_dist():
    lab={"domain":"DHX35–TOP1 domain (996 kb)","500kb":"gene ±500 kb",
         "1Mb":"gene ±1 Mb","1.5Mb":"gene ±1.5 Mb (published)"}
    r=[]
    for s in ("domain","500kb","1Mb","1.5Mb"):
        for i,mk in enumerate(("full","half","jaspar")):
            k=f"{s}|{mk}"
            if k not in A["dist"]: continue
            d=A["dist"][k]
            first=f'<td rowspan="3" class="grp">{lab[s]}<span class="sm2">{A["sets"][s]["n"]} DARs &middot; {A["sets"][s]["ra"]} RA / {A["sets"][s]["g4"]} G4</span></td>' if i==0 else ""
            r.append(f'<tr>{first}<td>{MN[mk]}</td><td class="n ra">{num(d["ra"])}</td>{pv(d["p_ra"])}'
                     f'<td class="n g4">{num(d["g4"])}</td><td class="n hl">{num(d["dod"])}</td>'
                     f'{pv(d["p_dod"])}<td class="n sm">{d["n_ra"]} / {d["n_g4"]}</td></tr>')
    return "\n".join(r)
def t_domic():
    r=[]
    for cols,cl in (("all","all columns"),("hi","high-IC only"),("lo","low-IC only")):
        for mk in ("full","jaspar"):
            k=f"{cols}|{mk}"
            if k not in A["dom_ic"]: continue
            d=A["dom_ic"][k]
            r.append(f'<tr><td>{cl}</td><td>{MN[mk]}</td><td class="n ra">{num(d["ra"])}</td>'
                     f'<td class="n g4">{num(d["g4"])}</td><td class="n hl">{num(d["dod"])}</td>'
                     f'{pv(d["p_dod"])}<td class="n sm">{d["n_ra"]} / {d["n_g4"]}</td></tr>')
    return "\n".join(r)
def t_dens():
    r=[]
    for mk in ("full","half","jaspar"):
        d=A["density"][mk]
        r.append(f'<tr><td>{MN[mk]}</td><td class="n">{d["ra"]:.1f}</td><td class="n">{d["g4"]:.1f}</td>'
                 f'<td class="n hl">{d["ra"]/d["g4"]:.2f}&times;</td></tr>')
    return "\n".join(r)
def t_ic():
    r=[]
    for mk in ("full","half","jaspar"):
        v=A["ic"][mk]
        cells="".join(f'<span class="icb {"hi" if x>=1 else ""}" style="--h:{min(x,2)/2*100:.0f}%" title="col {i+1}: {x:.2f} bits"></span>' for i,x in enumerate(v))
        hi=sum(1 for x in v if x>=1)
        r.append(f'<tr><td>{MN[mk]}<span class="sm2">{MW[mk]}</span></td><td class="icrow">{cells}</td>'
                 f'<td class="n sm">{hi} / {len(v)}</td></tr>')
    return "\n".join(r)

# ---------- part 2 tables ----------
def t_switch():
    lab={"domain":"DHX35–TOP1 domain","1.5Mb":"gene ±1.5 Mb"}
    al={"avian363":"avian363","74way":"74-way"}
    r=[]
    for a_ in ("avian363","74way"):
        for s in ("domain","1.5Mb"):
            for mk in ("full","jaspar"):
                k=f"{a_}|{s}|{mk}"
                if k not in B["switch"]: continue
                d=B["switch"][k]
                r.append(f'<tr><td>{al[a_]}<span class="sm2">{B["split"][a_]}</span></td><td>{lab[s]}</td>'
                         f'<td>{MN[mk]}</td>'
                         f'<td class="n ra">{d["sw_ra"]:.0f} / {d["te_ra"]:.0f}</td>'
                         f'<td class="n g4">{d["sw_g4"]:.0f} / {d["te_g4"]:.0f}</td>'
                         f'<td class="n hl">{num(d["diff"],4)}</td>{pv(d["p"])}</tr>')
    return "\n".join(r)
def t_swic():
    r=[]
    for mk in ("full","jaspar"):
        d=B["switch_ic"][mk]
        r.append(f'<tr><td>{MN[mk]}</td><td class="n">{d["sw"]:.2f} bits</td>'
                 f'<td class="n sm">{d["n_sw"]}</td><td class="n">{100*d["hi_sw"]:.0f}%</td>'
                 f'<td class="n">{d["ns"]:.2f} bits</td><td class="n sm">{d["n_ns"]}</td>'
                 f'<td class="n">{100*d["hi_ns"]:.0f}%</td></tr>')
    return "\n".join(r)
def t_swph():
    tr={"subtree222":"oscine","wholetree":"whole tree"}
    r=[]
    for mk in ("full","jaspar"):
        for t in ("subtree222","wholetree"):
            k=f"{mk}|{t}"
            if k not in B["switch_phylop"]: continue
            d=B["switch_phylop"][k]
            r.append(f'<tr><td>{MN[mk]}</td><td>{tr[t]}</td><td class="n">{num(d["sw"])}</td>'
                     f'<td class="n sm">{d["n_sw"]}</td><td class="n">{num(d["ns"])}</td>'
                     f'<td class="n sm">{d["n_ns"]}</td><td class="n hl">{num(d["diff"])}</td>'
                     f'{pv(d["p"])}</tr>')
    return "\n".join(r)
def t_pwm():
    lab={"domain":"DHX35–TOP1 domain","1.5Mb":"gene ±1.5 Mb"}
    r=[]
    for s in ("domain","1.5Mb"):
        for mk in ("full","jaspar"):
            k=f"{s}|{mk}"
            if k not in B["pwm"]: continue
            d=B["pwm"][k]
            for c,cc in (("Glut-RA","ra"),("GABA-4","g4")):
                v=d[c]
                r.append(f'<tr><td>{lab[s]}</td><td>{MN[mk]}</td><td class="{cc}">{c}</td>'
                         f'<td class="n">{num(v["site"],2)}</td><td class="n ctrl">{num(v["ctrl"],2)}</td>'
                         f'<td class="n hl">{num(v["paired"],2)}</td>{pv(v["p"])}'
                         f'<td class="n sm">{v["n"]}</td></tr>')
            r.append(f'<tr class="dodrow"><td></td><td></td><td>Glut-RA − GABA-4</td><td></td><td></td>'
                     f'<td class="n hl">{num(d["dod"]["d"],2)}</td>{pv(d["dod"]["p"])}<td></td></tr>')
    return "\n".join(r)

# ---------- part 3 tables ----------
def t_tf(key,cols):
    rows=sorted(B["tf"][key],key=lambda r:-r["dod"])
    r=[]
    for d in rows:
        hl=' class="arrow"' if d["tf"]=="AR" else ""
        c=[f'<tr{hl}><td>{d["tf"]}</td><td class="{"ra" if d["reg"]=="Glut-RA" else "g4"}">{d["reg"]}</td>']
        if "ra" in cols:
            c.append(f'<td class="n">{num(d["ra"])}</td>{pv(d["q_pra"])}')
            c.append(f'<td class="n">{num(d["g4"])}</td>{pv(d["q_pg4"])}')
        c.append(f'<td class="n hl">{num(d["dod"])}</td>{pv(d["q_p"])}'
                 f'<td class="n sm">{d["na"]:,} / {d["nb"]:,}</td></tr>')
        r.append("".join(c))
    return "\n".join(r)
def t_out():
    r=[]
    for d in B["outlier"]["gw"]["ereg"]:
        hl=' class="arrow"' if d["tf"].startswith("AR") else ""
        r.append(f'<tr{hl}><td>{d["tf"]}</td><td class="{"ra" if d["reg"]=="Glut-RA" else "g4"}">{d["reg"]}</td>'
                 f'<td class="n">{num(d["x"])}</td><td class="n">{num(d["y"])}</td>'
                 f'<td class="n hl">{num(d["res"])}</td><td class="n">{num(d["z_meta"],2)}</td>'
                 f'{pv(d["p_meta"])}{pv(d["q_meta"])}</tr>')
    return "\n".join(r)
fit=B["outlier"]["gw"]["fit"]
import sys as _sys
_sys.path.insert(0,F)
from span import elem_coords as _ec, in_domain as _ind, DS as _DS, DE as _DE
_EL=_ec()
GAIN=json.load(open(f"{M}/data/site_gain.json"))
SWEEP=json.load(open(f"{M}/data/gain_sweep.json"))
TREES=json.load(open(f"{M}/trees_gain/manifest.json"))
RATE=json.load(open(f"{M}/data74/rate_test_domain.json"))
RATEW=json.load(open(f"{M}/data74/rate_test_1p5mb.json"))
from scipy import stats as _st
for g in GAIN:
    g["gain"]=g["f_osc"]-g["f_out"]
    g["dom"]=_ind(_EL[g["elem"]])
    g["p"]=_st.fisher_exact([[g["p_osc"],g["n_osc"]-g["p_osc"]],
                             [g["p_out"],g["n_out"]-g["p_out"]]])[1]
GAIN_D=[g for g in GAIN if g["dom"]]
TREES_D=[t for t in TREES if _ind(_EL[t["elem"]])]
def t_gain():
    r=[]
    for g in sorted(GAIN_D,key=lambda x:-x["gain"]):
        r.append(f'<tr><td>{g["name"]}</td><td>{MN[g["motif"]]}</td>'
                 f'<td class="{"ra" if g["cls"]=="Glut-RA" else "g4"}">{g["cls"]}</td>'
                 f'<td class="n">{g["p_osc"]} / {g["n_osc"]}</td>'
                 f'<td class="n">{g["p_out"]} / {g["n_out"]}</td>'
                 f'<td class="n hl">{num(g["gain"],2)}</td>{pv(g["p"])}</tr>')
    return "\n".join(r)
def t_sweep2():
    r=[]
    for i,d in enumerate(SWEEP):
        first=f'<td rowspan="2" class="grp">{d["q"]}th</td>' if i%2==0 else ""
        r.append(f'<tr>{first}<td>{MN[d["motif"]]}</td><td class="n">{num(d["osc"],3)}</td>'
                 f'<td class="n">{num(d["out"],3)}</td><td class="n hl">{num(d["diff"],3)}</td>'
                 f'<td class="n sm">~{d["win"]:,}</td></tr>')
    return "\n".join(r)
def t_trees():
    r=[]
    for t in sorted(TREES_D,key=lambda x:(x["gains"],-(x["f_osc"]-x["f_out"]))):
        c=CNT.get((t["motif"],t["name"]),{})
        cls=' class="arrow"' if c.get("tier")=="oscine-stem" else ""
        loc=", ".join(c.get("where",[])) or "&mdash;"
        r.append(f'<tr{cls}><td>{t["name"]}</td><td>{MN[t["motif"]]}</td>'
                 f'<td class="{"ra" if t["cls"]=="Glut-RA" else "g4"}">{t["cls"]}</td>'
                 f'<td class="n">{t["f_osc"]:.2f}</td><td class="n">{t["f_out"]:.2f}</td>'
                 f'<td class="n hl">{t["gains"]}</td><td class="n">{t["losses"]}</td>'
                 f'<td class="sm">{loc}</td></tr>')
    return "\n".join(r)
def _rate_rows(R,span,first):
    G=R["gains"]; L=R["losses"]; BL=R["bl"]; stem=R["stem"]; path=set(R["path"])
    rows=[("oscine stem",[stem],1),
          ("rest of the root&rarr;<i>Lonchura</i> path",
           [k for k in BL if k in path and k!=stem],1),
          ("all branches off that path",[k for k in BL if k not in path],0)]
    out=[]
    for i,(lab,keys,hi) in enumerate(rows):
        g=sum(G.get(k,0) for k in keys); l=sum(L.get(k,0) for k in keys)
        b=sum(BL.get(k,0) for k in keys)
        sp=(f'<td rowspan="3" class="grp">{span}<span class="sm2">{R["n_sites"]} sites</span></td>'
            if i==0 else "")
        tr='<tr class="arrow">' if (hi and first) else '<tr>'
        out.append(f'{tr}{sp}<td>{lab}</td>'
                   f'<td class="n">{len(keys):,}</td><td class="n">{g}</td>'
                   f'<td class="n">{l}</td><td class="n">{b:,.1f}</td>'
                   f'<td class="n hl">{g/b:.2f}</td><td class="n">{l/b:.2f}</td></tr>')
    return out
def t_rate():
    return "\n".join(_rate_rows(RATE,"DHX35&ndash;TOP1 domain",True)
                     +_rate_rows(RATEW,"gene &plusmn;1.5 Mb",False))


# ---------- part 3: generated from the data ----------
CHOICE=json.load(open(f"{F}/results/regulator_motifs/tf_scenic_choice.json"))
def t_motifs():
    r=[]
    for tf,c in CHOICE.items():
        hl=' class="arrow"' if tf in ("AR","NR3C2") else ""
        r.append(f'<tr{hl}><td>{tf}</td><td class="{"ra" if c["regulon"]=="Glut-RA" else "g4"}">{c["regulon"]}</td>'
                 f'<td class="mono sm">{c["cluster"]}</td><td class="n">{c["nes"]:.1f}</td>'
                 f'<td class="sm">{c["how"]}</td><td class="mono">{c["consensus"]}</td></tr>')
    return "\n".join(r)
def _rank(rows,tf):
    o=sorted(rows,key=lambda r:-r["dod"]); n=[r["tf"] for r in o]
    return (n.index(tf)+1 if tf in n else None), len(o)
def part3():
    T=B["tf"]; G=T["subtree222_gw"]; W=T["wholetree_gw"]; Mx=T["subtree222_MAFB"]; Dm=T["subtree222_domain"]
    g={r["tf"]:r for r in G}; mx={r["tf"]:r for r in Mx}; dm={r["tf"]:r for r in Dm}
    sigG=[r["tf"] for r in G if r["q_p"]<.05]; sigW=[r["tf"] for r in W if r["q_p"]<.05]
    rG,nG=_rank(G,"AR"); rW,nW=_rank(W,"AR"); rM,nM=_rank(Mx,"AR"); rD,nD=_rank(Dm,"AR")
    spanG=max(r["dod"] for r in G)-min(r["dod"] for r in G)
    O=B["outlier"]["gw"]; ar={r["tf"]:r for r in O["ereg"]}
    msd=np.sqrt(np.median([a["vy"]+O["fit"]["b"]**2*a["vx"]-2*O["fit"]["b"]*a["cxy"] for a in O["arch"]]))
    od={k:{r["tf"]:r for r in B["outlier"][k]["ereg"]} for k in ("MAFB","domain")}
    def lst(x): return ", ".join(x) if x else "none"
    def arz(d,k): return (f'{k}: residual <span class="mono">{num(d[k]["res"])}</span>, z = {d[k]["z_meta"]:+.2f}, '
                          f'p = {d[k]["p_meta"]:.2f}') if k in d else f"{k}: not tested"
    return f"""<div class="defn">
<p>The same paired statistic applied to the motif of every TF that SCENIC+ nominates as an input
to MAFB &mdash; seven assigned to the Glut-RA regulon (AR, KLF, NR3C2, NFIA, LIN28B, STAT1,
FOSL2) and six to GABA-4 (MAF, SOX6, TCF12, ARX, SATB1, MEIS2) &mdash; asking whether AR's
cell-type asymmetry is distinctive among the factors proposed to act at this gene.</p>
<p><b>Matrices.</b> Each regulator, AR included, is represented by one matrix chosen from the motif
evidence SCENIC+ itself used for its cistrome, in the run that nominated it (glut-ra config4 for
the Glut-RA regulators, the GABA run's config4 for GABA-4; the two share every parameter): the
highest-NES motif annotated to the TF at the level of its MAFB eRegulon &mdash; direct for all but
SATB1, whose link is extended-only &mdash; and, where that motif is a metacluster, its medoid
member. One matrix per TF rather than a cluster's full membership keeps every TF at the same
per-base false-positive rate. All are scanned and measured by the same pipeline, at the element
level. <b>NR3C2 and AR are backed by the same cluster, metacluster_185.1, and so get the same
matrix:</b> no motif-based analysis can separate them, and NR3C2 is reported with AR. KLF12 is
reported as KLF. q is Benjamini-Hochberg within each table.</p>
<p><b>What changed.</b> An earlier version used a hand-assembled matrix set in which several
matrices belonged to another TF family (TCF12 a forkhead motif, STAT1 a TCF/LEF motif, SOX6 a
T-box motif), every matrix had lost its first column, and AR alone was measured with a different
motif and unit (HOMER ARE-full, per instance). That set is kept in
<span class="mono">superseded/tf13_motifs_2026-11-09/</span>. The ARE-full measurement appears
below as <b>AR*</b>, a labelled reference outside the comparison set.</p>
</div>
<div class="tw"><table>
<caption>the matrix representing each regulator</caption>
<thead><tr><th>TF</th><th>regulon</th><th>SCENIC+ motif</th><th class="n">NES</th><th>matrix</th>
<th>consensus</th></tr></thead>
<tbody>__MOTIFS__</tbody></table></div>
<p class="tnote">Two follow SCENIC+'s evidence to another family's motif: MAF's strongest cluster is
AP-1-like rather than a MARE (its MARE cluster ranks third), and STAT1's only direct evidence is an
AP-1 composite. LIN28B is an RNA-binding protein with no well-defined DNA motif. Several medoids are
short core motifs (7&ndash;8 bp for KLF, MEIS2, ARX, FOSL2, MAF, SOX6), which yields fewer sites at
FIMO p &lt; 1e-4.</p>

<h3>Genome-wide, oscine subtree</h3>
<div class="tw"><table>
<caption>motif minus rest of enhancer, by SCENIC+ regulon assignment</caption>
<thead><tr><th>TF</th><th>regulon</th><th class="n">Glut-RA</th><th class="n">q</th>
<th class="n">GABA-4</th><th class="n">q</th><th class="n">DoD</th><th class="n">q</th>
<th class="n">n elements</th></tr></thead>
<tbody>__TFO__</tbody></table></div>
<p class="tnote">AR ranks {rG} of {nG} on the Glut-RA&minus;GABA-4 difference, at
<span class="mono">{num(g["AR"]["dod"])}</span> (q = {g["AR"]["q_p"]:.2f}); the {nG} values span
{spanG:.3f}. Clearing FDR 5%: {lst(sigG)}. Under the superseded set AR ranked first at q = 0.032;
that value came from the ARE-full motif measured per instance, not from a matrix measured like
the others.</p>
<h3>Whole tree</h3>
<div class="tw"><table>
<caption>the same comparison computed pan-avian</caption>
<thead><tr><th>TF</th><th>regulon</th><th class="n">Glut-RA</th><th class="n">q</th>
<th class="n">GABA-4</th><th class="n">q</th><th class="n">DoD</th><th class="n">q</th>
<th class="n">n elements</th></tr></thead>
<tbody>__TFW__</tbody></table></div>
<p class="tnote">Pan-avian, AR ranks {rW} of {nW}; clearing FDR 5%: {lst(sigW)}.</p>
__FIG_TF__
<h3>In the DHX35&ndash;TOP1 domain</h3>
<div class="tw"><table>
<caption>default span, oscine subtree, regulators with &ge;8 elements per class</caption>
<thead><tr><th>TF</th><th>regulon</th><th class="n">DoD</th><th class="n">q</th>
<th class="n">n elements</th></tr></thead>
<tbody>__TFD__</tbody></table></div>
<p class="tnote">Only {nD} regulators have enough domain elements to test. AR ranks {rD} of {nD}, at
<span class="mono">{num(dm["AR"]["dod"])}</span> (q = {dm["AR"]["q_p"]:.2f}, {dm["AR"]["na"]}/{dm["AR"]["nb"]} elements).
Over the gene &plusmn;1.5 Mb window {nM} regulators qualify and AR ranks {rM}
(<span class="mono">{num(mx["AR"]["dod"])}</span>, q = {mx["AR"]["q_p"]:.2f}). Nothing clears correction
in either; the smallest q is {min(r["q_p"] for r in Dm):.2f} in the domain and
{min(r["q_p"] for r in Mx):.2f} over &plusmn;1.5 Mb.</p>
<h3>Against the 182-archetype background</h3>
<p>Widening the comparison set to 182 motif archetypes fixes a trend line between pan-avian and
oscine interaction (<span class="mono">r = __FITR__, p = __FITP__</span>), which the regulators
alone could not. Each TF is then tested as an outlier from that line, leave-one-out, weighted by its
own sampling error plus the between-motif dispersion
(<span class="mono">&tau; = __FITT__</span> against a median sampling SE of
<span class="mono">{msd:.4f}</span>).</p>
<div class="tw"><table>
<caption>genome-wide, by meta-regression p; AR = SCENIC+ matrix (= NR3C2), AR* = HOMER ARE-full per instance</caption>
<thead><tr><th>TF</th><th>regulon</th><th class="n">whole-tree</th><th class="n">oscine</th>
<th class="n">residual</th><th class="n">z</th><th class="n">p</th><th class="n">q</th></tr></thead>
<tbody>__OUT__</tbody></table></div>
<p class="tnote">No regulator is an outlier: the smallest q is {min(r["q_meta"] for r in O["ereg"]):.2f}.
Genome-wide, {arz(ar,"AR")} &mdash; on the opposite side of the line from the ARE-full
reference, {arz(ar,"AR*")}. The same holds over &plusmn;1.5 Mb ({arz(od["MAFB"],"AR")}) and in the
domain ({arz(od["domain"],"AR")}), though at the locus there is no usable trend: it collapses to
r = {B["outlier"]["MAFB"]["fit"]["r"]:+.2f} over &plusmn;1.5 Mb and
r = {B["outlier"]["domain"]["fit"]["r"]:+.2f} in the domain, with &tau; = 0, so sampling error
swamps any between-motif signal.</p>
__FIG_OUT__"""

# ---------- part 4: generated, so the prose cannot drift from the data ----------
CNT={(c["motif"],c["name"]):c for c in json.load(open(f"{M}/trees363/counts_ar.json"))}
R363=json.load(open(f"{M}/data363/rate_test363.json"))
def t_place():
    """AR sites whose principal oscine gain sits at or near the oscine stem (avian363)."""
    rows=sorted([c for c in CNT.values() if c["tier"] in ("stem","near-stem")],
                key=lambda c:(not c["dom"],c["depth"],c["nonosc_gains"]))
    r=[]
    for c in rows:
        hl=' class="arrow"' if c["nonosc_gains"]==0 else ""
        where_="stem" if c["depth"]==0 else f'near-stem, {c["depth"]} node{"s" if c["depth"]>1 else ""} in'
        r.append(f'<tr{hl}><td>{c["name"].split("_",1)[1]}</td><td>{MN[c["motif"]]}</td>'
                 f'<td class="{"ra" if c["cls"]=="Glut-RA" else "g4"}">{c["cls"]}</td>'
                 f'<td>{"domain" if c["dom"] else "&plusmn;1.5 Mb"}</td><td>{where_}</td>'
                 f'<td class="n">{100*c["osc_cov"]:.0f}%</td>'
                 f'<td class="n">{c["f_osc"]:.2f}</td><td class="n">{c["f_out"]:.2f}</td>'
                 f'<td class="n">{c["other_osc_gains"]}</td><td class="n hl">{c["nonosc_gains"]}</td></tr>')
    return "\n".join(r)
def t_rate363():
    r=[]
    for fam,lab in (("AR","AR"),("TF","other regulators")):
        if fam not in R363["rate"]: continue
        for sp,spl in (("domain","DHX35&ndash;TOP1 domain"),("1p5mb","gene &plusmn;1.5 Mb")):
            d=R363["rate"][fam][sp]
            for k,kl in (("stem","oscine stem"),("stem+backbone","stem + near-stem backbone")):
                x=d[k]; hl=' class="arrow"' if fam=="AR" else ""
                r.append(f'<tr{hl}><td>{lab}</td><td>{spl}<span class="sm2">{d["n_sites"]} sites</span></td><td>{kl}</td>'
                         f'<td class="n">{x["gains"]}</td><td class="n">{x["length"]:.2f}</td>'
                         f'<td class="n">{x["exp_matched"]:.2f}</td><td class="n hl">{x["fold_matched"]:.2f}&times;</td>'
                         f'{pv(x["p_matched"])}</tr>')
    return "\n".join(r)
def _bh(pv):
    pv=np.asarray(pv,float); o=np.argsort(pv); q=np.empty(len(pv))
    q[o]=np.minimum(np.minimum.accumulate((pv[o]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1],1)
    return q
for g,q in zip(GAIN_D,_bh([g["p"] for g in GAIN_D])): g["q_dom"]=q
for g,q in zip(GAIN,_bh([g["p"] for g in GAIN])): g["q_win"]=q
def _tier(g):
    c=CNT.get((g["motif"],g["name"])); return c["tier"] if c else "no coverage"
def _rate(R):
    from scipy import stats as st_
    G=R["gains"]; L=R["losses"]; BL=R["bl"]; stem=R["stem"]; path=set(R["path"])
    rest=[k for k in BL if k in path and k!=stem]; off=[k for k in BL if k not in path]
    g0=G.get(stem,0); b0=BL[stem]
    gr=sum(G.get(k,0) for k in rest)/sum(BL[k] for k in rest)
    go=sum(G.get(k,0) for k in off)/sum(BL[k] for k in off)
    e=b0*gr
    return dict(n=R["n_sites"],g0=g0,exp=e,fold=g0/e,p=st_.poisson.sf(g0-1,e),
                naive=g0/(b0*go),
                lpath=sum(L.get(k,0) for k in rest)+L.get(stem,0),
                loff=sum(L.get(k,0) for k in off))
def _pfmt(p): return f"{p:.2g}" if p>=0.001 else "&lt; 0.001"
def part4():
    cand_d=[g for g in GAIN_D if g["gain"]>=.5 and g["q_dom"]<.05]
    cand_w=[g for g in GAIN if g["gain"]>=.5 and g["q_win"]<.05]
    tiers=collections.Counter(_tier(g) for g in cand_d)
    stem_d=[g for g in cand_d if _tier(g) in ("stem","near-stem")]
    zero=[g for g in cand_d if g["p_out"]==0]
    extra=len(cand_w)-len(cand_d)
    x598=[g for g in GAIN if g["name"]=="RA_Glut-RA.395598_1" and g["motif"]=="full"]
    # how much the corrected call moved
    old={(g["motif"],g["name"]):g for fn in ("site_gain_concat","tf_site_gain_concat")
         for g in json.load(open(f"{M}/data/{fn}.json"))}
    cov=json.load(open(f"{M}/data/site_cov.json"))
    both=[c for c in cov.values() if (c["motif"],c["name"]) in old]
    moved=sum(1 for c in both if abs(c["p_osc"]-old[(c["motif"],c["name"])]["p_osc"])
              +abs(c["p_out"]-old[(c["motif"],c["name"])]["p_out"])>=5)
    nocov=[g["n_out_nocov"]/(g["n_out"]+g["n_out_nocov"]) for g in GAIN if g["n_out"]+g["n_out_nocov"]]
    RA_=R363["rate"]["AR"]; RT_=R363["rate"].get("TF",{})
    def rr(fam,sp,k): return R363["rate"][fam][sp][k] if fam in R363["rate"] else None
    clean=[c for c in CNT.values() if c["tier"] in ("stem","near-stem") and c["nonosc_gains"]==0]
    sw={(r["q"],r["motif"]):r for r in SWEEP}
    jd=[r["diff"] for r in SWEEP if r["motif"]=="jaspar"]
    def nm(g): return g["name"].split("_",1)[1].rsplit("_",1)[0]
    def stemlist():
        if not stem_d: return "none"
        return ", ".join(f'<span class="mono">{nm(g)}</span> ({MN[g["motif"]]}, '
                         f'{"stem" if CNT[(g["motif"],g["name"])]["depth"]==0 else str(CNT[(g["motif"],g["name"])]["depth"])+" nodes in"})'
                         for g in stem_d)
    tsum=", ".join(f"{v} {k}" for k,v in sorted(tiers.items(),key=lambda kv:-kv[1]))
    out=f"""<div class="defn">
<p><b>The question changes.</b> The analyses above ask whether bases are held still. This one
asks whether a scoring AR site <i>exists</i> in each lineage. That makes indels a mechanism of
gain rather than a confound: a deletion that brings two half-sites into register creates a site
exactly as a substitution does.</p>
<p><b>Method.</b> For each AR site called on lonStrDom2, each species' sequence across a
&plusmn;30 bp window is projected onto the reference &mdash; every reference position taken once,
from the block where that species aligns best, with insertions kept between consecutive
positions &mdash; and scanned with the PWM on both strands. A species has the site if its best hit
clears the 99.99th percentile of the matrix's uniform-sequence null, the stringency FIMO's
p &lt; 1e-4 implies. A species with bases at fewer than 80% of the motif's positions has no
orthologous sequence there and is <b>missing, not absent</b>: it is dropped from the
frequencies, and on the tree no gain or loss is called on any branch that touches it.</p>
<p><b>Correction.</b> An earlier version joined every aligned row of every overlapping block.
That spliced duplicated spans and paralogous rows into sequence no genome has &mdash; it can make
motif hits and break real ones &mdash; and it scored genomes with no sequence at the site as
lacking it. Across all {len(both)} motif sites at the locus, AR and the other predicted inputs together,
presence moved by five or more genomes at {moved}. Among the AR sites scored here, non-oscines
with no sequence at the motif were a median {100*np.median(nocov):.0f}% of the outgroup per site,
up to {100*max(nocov):.0f}%. Every number in this section is on the corrected
call.</p>
</div>

<h3>Candidate gains</h3>
<div class="tw"><table>
<caption>every site in the default span, avian363, ranked by gain</caption>
<thead><tr><th>site</th><th>motif</th><th>class</th><th class="n">oscine</th>
<th class="n">outgroup</th><th class="n">gain</th><th class="n">p</th></tr></thead>
<tbody>__GAIN__</tbody></table></div>
<p class="tnote">{len(GAIN_D)} sites in the default span, ranked; {len(cand_d)} clear gain
&ge; 0.50 and within-domain q &lt; 0.05, {len(zero)} of them absent from every non-oscine with
sequence. The reference is itself an oscine and sites were called on it, so oscine presence is
inflated before any biology; the ranking is informative but the fractions are not an unbiased
gain rate. Fisher p treats genomes as independent and only ranks. Widening to &plusmn;1.5 Mb adds
{extra} further candidates""" + (f""", among them the strongest in the whole set
(<span class="mono">Glut-RA.395598</span>, ARE-full, {x598[0]["p_osc"]}/{x598[0]["n_osc"]} against
{x598[0]["p_out"]}/{x598[0]["n_out"]}, {_tier(x598[0])}), but those sit beyond TOP1 and are not
MAFB-domain results.</p>""" if x598 else ".</p>") + f"""

<h3>Reverse-ascertainment calibration</h3>
<p>Letting each of {sw[("99.99","full")]["n_osc"]+sw[("99.99","full")]["n_out"]} genomes
({sw[("99.99","full")]["n_osc"]} oscine, {sw[("99.99","full")]["n_out"]} non-oscine) in turn play
the role lonStrDom2 played gives the expected
asymmetry under no real gain. Near and far sets are matched on <i>patristic distance</i> rather
than on species count, because the familiar outgroups
(<i>Gallus</i>, <i>Struthio</i>, <i>Falco</i>) are too isolated in this tree to supply a matched
comparison; the outgroup arm is therefore seabirds and waterbirds, the only densely sampled
non-oscine neighbourhoods available. Coverage applies per window, in the ascertainment genome
and in both arms.</p>
<div class="tw"><table>
<caption>oscine-ascertained minus outgroup-ascertained gain, by site-calling threshold</caption>
<thead><tr><th>threshold</th><th>motif</th><th class="n">oscine arm</th>
<th class="n">outgroup arm</th><th class="n">difference</th><th class="n">windows/genome</th></tr></thead>
<tbody>__SWEEP2__</tbody></table></div>
<p class="tnote">ARE-full's apparent excess is <span class="mono">{sw[("99.99","full")]["diff"]:+.3f}</span>
at the strictest threshold, on ~{sw[("99.99","full")]["win"]} windows per genome, and falls to
<span class="mono">{sw[("99.9","full")]["diff"]:+.3f}</span> once the sample grows to
~{sw[("99.9","full")]["win"]} &mdash; the behaviour of an estimate riding on a handful of extreme
windows, and ARE-full is the matrix discovered on Glut-RA sequence. JASPAR, which carries no
discovery bias, stays between <span class="mono">{min(jd):+.3f}</span> and
<span class="mono">{max(jd):+.3f}</span> at every threshold. At the clade level, songbirds are not
unusually prone to gaining AR sites. Whether the effect is specific to the strongest-scoring
tier is not resolvable at this sample size.</p>

<h3>Where the gains sit on the full avian tree</h3>
<div class="defn">
<p><b>Placement is read on the avian363 tree.</b> An earlier version used the 74-way
alignment, whose seventeen oscines are all Passerides or <i>Corvus</i>: it has no basal oscines,
so its "oscine stem" was really the core-oscine stem. avian363 samples the basal lineages
&mdash; lyrebird and scrubbird, treecreeper and bowerbird, fairy-wren, pardalotes and honeyeaters,
logrunner and Australasian babbler &mdash; so a gain can be placed on the true stem or a few
nodes inside it.</p>
<p><b>Near the stem counts.</b> Song and the song system are unlikely to have arisen in one
genetic event, so a gain a few nodes inside the stem is as relevant as one exactly on it. The
deep oscine backbone is defined by coverage: the four nodes between the oscine ancestor and the
core oscines each still hold &ge;90% of oscines (99, 97, 94 and 92%), and coverage falls to 71% at
the next node, where the Corvides split. Each site is placed by its <i>principal</i> oscine gain,
the one with most oscines below it; independent gains outside oscines are counted as homoplasy,
not used to discard the site.</p>
</div>
__FIG_SITEALN__
<div class="tw"><table>
<caption>AR sites whose principal oscine gain is at or near the oscine stem (avian363); shaded = no gain outside oscines</caption>
<thead><tr><th>site</th><th>motif</th><th>class</th><th>span</th><th>gained at</th>
<th class="n">oscines below</th><th class="n">oscine</th><th class="n">outgroup</th>
<th class="n">re-gains in oscines</th><th class="n">gains outside oscines</th></tr></thead>
<tbody>__PLACE__</tbody></table></div>
<p class="tnote">Of the {len(cand_d)} domain candidates the placements are: {tsum}. At or near the
oscine stem: {stemlist()}. Across all AR sites, {len(clean)} gained at or near the stem with no
independent gain outside oscines: {", ".join(f'<span class="mono">{c["name"].split("_",1)[1]}</span>' for c in clean) or "none"}.
Several stem-placed sites are oscine-enriched but not oscine-specific &mdash; present in nearly half
of non-oscines too &mdash; which the outgroup column shows.</p>

<h3>Is the gain rate elevated at the base of oscines?</h3>
<div class="tw"><table>
<caption>gains per unit branch length, summed over sites, against the rest of the root-to-<i>Lonchura</i> path (ascertainment-matched)</caption>
<thead><tr><th>set</th><th>span</th><th>branches</th><th class="n">gains</th><th class="n">length</th>
<th class="n">expected</th><th class="n">fold</th><th class="n">p</th></tr></thead>
<tbody>__RATE363__</tbody></table></div>
<p class="tnote">Only the ascertainment-matched comparison means anything: every branch on the
path to the reference must end in a site-bearing genome, so gains there are over-represented and
losses nearly impossible, and comparing the stem with the rest of the tree is
{RA_["domain"]["stem+backbone"]["fold_off"]:.0f}&times; inflated in the domain for that reason alone.
Against the matched branches, AR gains on the stem and backbone together are
{RA_["domain"]["stem+backbone"]["fold_matched"]:.2f}&times; expectation in the domain
(p = {_pfmt(RA_["domain"]["stem+backbone"]["p_matched"])}) and
{RA_["1p5mb"]["stem+backbone"]["fold_matched"]:.2f}&times; over &plusmn;1.5 Mb
(p = {_pfmt(RA_["1p5mb"]["stem+backbone"]["p_matched"])}); on the stem branch alone,
{RA_["domain"]["stem"]["gains"]} gains in the domain and {RA_["1p5mb"]["stem"]["gains"]} over &plusmn;1.5 Mb.""" + (f"""
The other regulators, in the same test, give
{RT_["domain"]["stem+backbone"]["fold_matched"]:.2f}&times; (domain) and
{RT_["1p5mb"]["stem+backbone"]["fold_matched"]:.2f}&times; (&plusmn;1.5 Mb).""" if RT_ else "") + """</p>

<div class="verdict">
<h3><span class="chip open">no excess at the base of oscines</span>The rate result</h3>
<p>The 74-way test found AR gains concentrated on its stem branch (2.0&times;, p = 0.011 over
&plusmn;1.5 Mb). With the basal oscines sampled, that branch turns out to be the core-oscine stem, and
the gains are distributed down the backbone at the rate the matched branches predict. The base of
the oscine radiation does not gain AR sites faster than the rest of the reference lineage. What
remains are individual sites gained there &mdash; candidates for function, not evidence of an
elevated rate.</p>
</div>
"""
    return out

SY=json.load(open(f"{M}/data/ascertain_sym2.json"))
SWP=json.load(open(f"{M}/data/ascertain_sweep.json"))
QU=json.load(open(f"{M}/data/score_quantiles.json"))
def t_asc():
    r=[]
    for mk in ("full","jaspar"):
        for c,cc in (("Glut-RA","ra"),("GABA-4","g4")):
            k=f"{mk}|{c}"
            if k not in SY: continue
            d=SY[k]
            lo,hi=d["ci_mid"]
            r.append(f'<tr><td>{MN[mk]}</td><td class="{cc}">{c}</td>'
                     f'<td class="n">{num(d["fwd"],2)}</td><td class="n">{num(d["rev"],2)}</td>'
                     f'<td class="n hl">{num(d["mid"],2)}</td>'
                     f'<td class="n sm">[{num(lo,2)}, {num(hi,2)}]</td>'
                     f'<td class="n sm">{d["n_fwd"]} / {d["n_rev"]}</td></tr>')
    return "\n".join(r)
def t_sweep():
    r=[]
    for pct in ("98.0","99.0","99.5","99.8","99.9"):
        for i,mk in enumerate(("full","jaspar")):
            k=f"{pct}|{mk}"
            if k not in SWP: continue
            d=SWP[k]; lo,hi=d["ci"]
            sig=' class="arrow"' if lo>0 else ""
            first=f'<td rowspan="2" class="grp">{pct}th</td>' if i==0 else ""
            r.append(f'<tr{sig}>{first}<td>{MN[mk]}</td><td class="n sm">{num(d["cut"],2)}</td>'
                     f'<td class="n ra">{num(d["ra"],2)}</td><td class="n g4">{num(d["g4"],2)}</td>'
                     f'<td class="n hl">{num(d["diff"],2)}</td>'
                     f'<td class="n sm">[{num(lo,2)}, {num(hi,2)}]</td>'
                     f'<td class="n sm">{d["n"][0]}/{d["n"][1]} &middot; {d["n"][2]}/{d["n"][3]}</td></tr>')
    return "\n".join(r)
def t_quant():
    r=[]
    for mk in ("full","jaspar"):
        for i,p in enumerate(("50","90","99","99.5","99.9","99.99")):
            d=QU[f"{mk}|{p}"]
            first=f'<td rowspan="6" class="grp">{MN[mk]}<span class="sm2">{QU[f"{mk}|n"]:,} windows</span></td>' if i==0 else ""
            r.append(f'<tr>{first}<td class="n">P{p}</td><td class="n ra">{num(d["ra"],2)}</td>'
                     f'<td class="n g4">{num(d["g4"],2)}</td>'
                     f'<td class="n hl">{num(d["diff"],2)}</td></tr>')
    return "\n".join(r)

HTML = """<meta name="claude:viewport-fit-sides" content="cover">
<title>AR Motif Evolution at MAFB</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Serif:ital,wght@0,400;0,500;1,400&display=swap">
<style>
/* Single column, prose at 68ch, tables allowed to run wider in their own scrollers. */
:root{
  --paper:#FBFAF8; --ink:#1A1D23; --ink2:#4A515E; --ink3:#767E8C;
  --rule:#E2E0DB; --panel:#F3F1EC; --ra:#C8161A; --g4:#2C6E7F; --warn:#9A6B1E;
  --sans:"IBM Plex Sans",system-ui,sans-serif;
  --serif:"IBM Plex Serif",Georgia,serif;
  --mono:"IBM Plex Mono",ui-monospace,monospace;
  color-scheme:light;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --paper:#15171B; --ink:#E8E6E1; --ink2:#AEB3BC; --ink3:#848B96;
  --rule:#2B2F36; --panel:#1D2026; --ra:#F2575A; --g4:#5FA8BA; --warn:#D7A64B;
  color-scheme:dark;}}
:root[data-theme="dark"]{
  --paper:#15171B; --ink:#E8E6E1; --ink2:#AEB3BC; --ink3:#848B96;
  --rule:#2B2F36; --panel:#1D2026; --ra:#F2575A; --g4:#5FA8BA; --warn:#D7A64B;
  color-scheme:dark;}
body{background:var(--paper);color:var(--ink);font-family:var(--serif);
  font-size:16px;line-height:1.62;margin:0}
.wrap{max-width:1040px;margin:0 auto;
  padding-left:max(16px,env(safe-area-inset-left,0px));
  padding-right:max(16px,env(safe-area-inset-right,0px));
  padding-block:40px 72px}
h1,h2,h3,.eyebrow,th,.chip{font-family:var(--sans)}
h1{font-size:2.05rem;line-height:1.14;font-weight:700;letter-spacing:-.02em;
  margin:.1em 0 .35em;text-wrap:balance}
h2{font-size:1.3rem;font-weight:650;letter-spacing:-.01em;margin:0 0 .15em;text-wrap:balance}
h3{font-size:1.02rem;font-weight:600;margin:2.2em 0 .5em;text-wrap:balance}
p{margin:0 0 1em;max-width:68ch}
.eyebrow{font-size:.705rem;text-transform:uppercase;letter-spacing:.11em;
  color:var(--ink3);font-weight:600}
.lede{font-size:1.06rem;color:var(--ink2);max-width:66ch}
.sect{border-top:2px solid var(--ink);padding-top:14px;margin-top:56px}
.sect:first-of-type{margin-top:40px}
.secthead{display:flex;gap:18px;align-items:baseline;flex-wrap:wrap;margin-bottom:.9em}
.secthead .no{font-family:var(--mono);font-size:.78rem;color:var(--ink3);
  border:1px solid var(--rule);border-radius:3px;padding:2px 7px;flex:none}
.defn{background:var(--panel);border-left:3px solid var(--ink3);
  padding:14px 18px;margin:0 0 1.4em;border-radius:0 3px 3px 0}
.defn p{margin:0;font-size:.93rem;color:var(--ink2);max-width:72ch}
.defn p+p{margin-top:.6em}
.tw{overflow-x:auto;margin:0 0 .6em;min-width:0}
table{border-collapse:collapse;width:100%;font-family:var(--sans);font-size:.83rem}
caption{font-family:var(--sans);font-size:.78rem;color:var(--ink3);text-align:left;
  padding:0 0 7px;font-weight:500}
th{font-weight:600;font-size:.715rem;text-transform:uppercase;letter-spacing:.055em;
  color:var(--ink3);text-align:left;padding:0 10px 6px 0;border-bottom:1px solid var(--rule);
  white-space:nowrap}
td{padding:5px 10px 5px 0;border-bottom:1px solid var(--rule);vertical-align:baseline}
th.n,td.n{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono);
  font-size:.79rem;white-space:nowrap}
td.hl{font-weight:600}
td.sm{color:var(--ink3);font-size:.72rem}
.sm2{display:block;color:var(--ink3);font-size:.7rem;font-weight:400;
  font-family:var(--sans);letter-spacing:0;text-transform:none}
td.grp{font-weight:600;border-right:1px solid var(--rule);padding-right:14px;vertical-align:top}
td.ra,span.ra{color:var(--ra)} td.g4,span.g4{color:var(--g4)}
td.sig{color:var(--ink);font-weight:600} td.near{color:var(--warn)} td.ns{color:var(--ink3)}
td.ctrl{color:var(--ink3)}
tr.arrow td{background:color-mix(in srgb,var(--ra) 7%,transparent)}
tr.dodrow td{border-bottom:2px solid var(--rule);font-style:italic;color:var(--ink2)}
.tnote{font-family:var(--sans);font-size:.76rem;color:var(--ink3);max-width:76ch;
  margin:.2em 0 2.2em;line-height:1.55}
figure{margin:1.6em 0 2em}
figure img{display:block;width:100%;height:auto;border:1px solid var(--rule);border-radius:2px;
  background:#fff}
figcaption{font-family:var(--sans);font-size:.78rem;color:var(--ink3);margin-top:8px;
  max-width:76ch;line-height:1.55}
.icrow{display:flex;gap:2px;align-items:flex-end;height:26px}
.icb{width:11px;background:var(--rule);height:var(--h);min-height:2px;border-radius:1px}
.icb.hi{background:var(--ink2)}
.key{display:flex;gap:22px;flex-wrap:wrap;font-family:var(--sans);font-size:.76rem;
  color:var(--ink3);margin:0 0 2.4em;padding-top:4px}
.key b{color:var(--ink);font-weight:600}
.verdict{border:1px solid var(--rule);border-radius:3px;padding:16px 20px;margin:1.4em 0 2em;
  background:var(--panel)}
.verdict h3{margin:0 0 .4em;font-size:.95rem}
.verdict p{font-size:.93rem;color:var(--ink2);margin:0 0 .7em;max-width:72ch}
.verdict p:last-child{margin:0}
.chip{display:inline-block;font-size:.67rem;text-transform:uppercase;letter-spacing:.07em;
  font-weight:600;padding:2px 7px;border-radius:2px;margin-right:8px;vertical-align:1px}
.chip.hold{background:var(--g4);color:var(--paper)}
.chip.fall{background:var(--ink3);color:var(--paper)}
.chip.cand{background:var(--ra);color:var(--paper)}
code,.mono{font-family:var(--mono);font-size:.86em}
.foot{border-top:1px solid var(--rule);margin-top:56px;padding-top:16px;
  font-family:var(--sans);font-size:.75rem;color:var(--ink3);line-height:1.7}
.foot code{font-size:.72rem}
@media (max-width:560px){
  body{font-size:15px} h1{font-size:1.6rem} .wrap{padding-block:28px 48px}
  .secthead{gap:10px}
}
</style>

<div class="wrap">
<p class="eyebrow">Songbird enhancer evolution &middot; supporting analysis</p>
<h1>AR Motif Evolution at MAFB</h1>
<p class="lede">Three motif-level analyses of androgen-receptor binding sites in Glut-RA and
GABA-4 differentially accessible regions, run genome-wide and across progressively tighter
definitions of the MAFB regulatory locus.</p>

<div class="key">
  <span><b>Glut-RA</b> RA over Arcopallium</span>
  <span><b>GABA-4</b> GABA-4 over GABA-1</span>
  <span><b>oscine</b> birdAnc222, 146 species</span>
  <span><b>outgroup</b> 217 species</span>
  <span><b>MAFB</b> chr20:13,365,864&ndash;13,369,534 (&minus;)</span>
</div>

<section class="sect">
<div class="secthead"><span class="no">1</span><div>
  <h2>Motif versus the rest of its enhancer</h2>
  <p class="eyebrow">phyloP contrast, paired within element</p></div></div>

<div class="defn">
<p><b>Metric.</b> For each motif site, the mean phyloP across the motif match minus the mean
phyloP across the rest of the same DAR, where &ldquo;rest&rdquo; is the element minus every base any
hit of that motif touches. Pairing within the element removes between-element differences in
GC, repeat content, alignability and cell-type identity. Positive means the motif is held more
tightly than the sequence around it.</p>
<p><b>Unit.</b> One value per merged site, so the two strand reports of a palindromic ARE
collapse into one. <b>DoD</b> is the Glut-RA value minus the GABA-4 value. Significance from
20,000 sign-flip (per class) or label permutations (DoD).</p>
</div>

<h3>Site density: the half site is the control</h3>
<div class="tw"><table>
<caption>merged motif sites per Mb of DAR sequence</caption>
<thead><tr><th>motif</th><th class="n">Glut-RA</th><th class="n">GABA-4</th><th class="n">ratio</th></tr></thead>
<tbody>__DENS__</tbody></table></div>
<p class="tnote">AR binds as a dimer on the inverted repeat, so only full sites are AR-specific;
half sites are shared across the nuclear-receptor family. The full-site motifs are enriched in
Glut-RA DARs, the half site is not.</p>

<h3>Information content per column</h3>
<div class="tw"><table>
<caption>bars are per-column IC, 0&ndash;2 bits; filled bars are &ge;1 bit</caption>
<thead><tr><th>motif</th><th>columns</th><th class="n">high-IC</th></tr></thead>
<tbody>__IC__</tbody></table></div>
<p class="tnote">ARE-half has no column reaching 1 bit, so the IC-restricted analyses below
cannot be run on it.</p>

<h3>Genome-wide, oscine subtree</h3>
<div class="tw"><table>
<caption>motif minus rest of enhancer, by which motif columns enter the numerator</caption>
<thead><tr><th>columns</th><th>motif</th><th class="n">Glut-RA</th><th class="n">p</th>
<th class="n">GABA-4</th><th class="n">p</th><th class="n">DoD</th><th class="n">p</th>
<th class="n">n sites</th></tr></thead>
<tbody>__GWO__</tbody></table></div>
<p class="tnote">The contrast lives entirely in the informative columns: restricting to
high-IC bases raises the DoD for both full-site motifs, and the degenerate columns go flat.
That is the behaviour a genuine binding-site constraint signal should show.</p>

<h3>Genome-wide, whole tree</h3>
<div class="tw"><table>
<caption>the same statistic computed pan-avian</caption>
<thead><tr><th>columns</th><th>motif</th><th class="n">Glut-RA</th><th class="n">p</th>
<th class="n">GABA-4</th><th class="n">p</th><th class="n">DoD</th><th class="n">p</th>
<th class="n">n sites</th></tr></thead>
<tbody>__GWW__</tbody></table></div>
<p class="tnote">ARE-full is negative pan-avian in both classes while JASPAR is positive; that
disagreement between two full-site matrices is unexplained and is the reason the
tree-to-tree sign reversal should not be presented as a property of AR sites.</p>

__FIG_FOLD__

<h3>Narrowing the locus definition</h3>
<div class="tw"><table>
<caption>oscine subtree, all motif columns, by element set</caption>
<thead><tr><th>element set</th><th>motif</th><th class="n">Glut-RA</th><th class="n">p</th>
<th class="n">GABA-4</th><th class="n">DoD</th><th class="n">p</th><th class="n">n sites</th></tr></thead>
<tbody>__DIST__</tbody></table></div>
<p class="tnote">The published set is the gene &plusmn;1.5 Mb, which is a &plusmn;500 kb window
expanded by 1 Mb for statistical power. The DHX35&ndash;TOP1 domain is the gene-bounded
alternative: 996 kb of intergenic sequence containing MAFB and two putative lncRNAs, and
balanced in class composition. The ARE-full Glut-RA value declines monotonically as the window
is pulled in toward the gene, reaching significance only at the widest definition.</p>

<h3>Column restriction inside the domain</h3>
<div class="tw"><table>
<caption>DHX35&ndash;TOP1 domain, oscine subtree</caption>
<thead><tr><th>columns</th><th>motif</th><th class="n">Glut-RA</th><th class="n">GABA-4</th>
<th class="n">DoD</th><th class="n">p</th><th class="n">n sites</th></tr></thead>
<tbody>__DOMIC__</tbody></table></div>
<p class="tnote">Inside the domain the pattern inverts: the high-IC columns give the weaker
contrast and the degenerate columns the stronger one, the opposite of the genome-wide result.
On 4&ndash;11 sites per class nothing is significant, but the direction argues against the
domain signal being a weak version of the genome-wide one.</p>

<div class="verdict">
<h3><span class="chip hold">holds</span>Genome-wide, class-level</h3>
<p>AR full sites are enriched in Glut-RA DARs and more constrained relative to their own
enhancer in Glut-RA than in GABA-4. The effect concentrates on the high-information columns,
and the motif-matched but non-AR-specific half site shows neither the enrichment nor the
asymmetry.</p>
</div>
<div class="verdict">
<h3><span class="chip fall">does not hold</span>Locus-level</h3>
<p>The contrast weakens as the locus is defined more precisely and is weakest inside MAFB's
own regulatory domain, where the column restriction also runs the wrong way. The published
&plusmn;1.5 Mb figures should not be described as a MAFB-locus result.</p>
</div>
</section>

<section class="sect">
<div class="secthead"><span class="no">2</span><div>
  <h2>Nucleotide switches and motif scoring</h2>
  <p class="eyebrow">lineage-specific change, not constraint</p></div></div>

<div class="defn">
<p><b>Switch.</b> An alignment column where the oscine majority base differs from the outgroup
majority base, with &ge;80% support on each side and adequate coverage (&ge;40 oscine / &ge;60
outgroup for avian363; &ge;6 / &ge;15 for the 74-way). Gaps and uncalled bases are excluded from
both the numerator and the denominator. Rates are switches per <i>testable column</i>, so sites
with poorer alignment coverage do not distort them.</p>
<p><b>Motif scoring.</b> The graded alternative: build the oscine and outgroup consensus
sequences, score each as PWM log-odds against a uniform background on the better strand, and
take the difference. A change at a high-information column moves the score a lot, one at a
degenerate column barely at all, and no arbitrary threshold is needed.</p>
</div>

<h3>Switch rate, Glut-RA versus GABA-4</h3>
<div class="tw"><table>
<caption>switches / testable columns; p from 20,000 permutations of DAR class across sites</caption>
<thead><tr><th>alignment</th><th>element set</th><th>motif</th><th class="n">Glut-RA</th>
<th class="n">GABA-4</th><th class="n">rate diff</th><th class="n">p</th></tr></thead>
<tbody>__SW__</tbody></table></div>
<p class="tnote">The Glut-RA excess is significant across the &plusmn;1.5 Mb set in three of
four alignment &times; motif combinations and in none of the domain-restricted tests. The domain
sets are roughly half the size, so loss of significance is partly a power effect; what
distinguishes ARE-full is that the point estimate also falls by half or more, which power loss
alone does not explain.</p>

__FIG_SW__

<h3>Switches happen at unconserved bases</h3>
<div class="tw"><table>
<caption>per-base phyloP at switching versus non-switching motif columns, avian363 calls</caption>
<thead><tr><th>motif</th><th>tree</th><th class="n">switch cols</th><th class="n">n</th>
<th class="n">other cols</th><th class="n">n</th><th class="n">difference</th><th class="n">p</th></tr></thead>
<tbody>__SWPH__</tbody></table></div>
<p class="tnote">Switch columns sit near zero on both trees and below non-switch columns in all
four comparisons. This is close to definitional &mdash; a switch requires the base to have
changed, and phyloP scores change against a neutral model &mdash; and it means the constraint
and switch analyses measure opposite things at the same position. In the worked example below,
the switching base is the one unconserved column in a site whose pan-avian mean is
<span class="mono">+1.95</span>.</p>

<h3>Switching columns are not degenerate columns</h3>
<div class="tw"><table>
<caption>information content of switching versus non-switching columns</caption>
<thead><tr><th>motif</th><th class="n">switch IC</th><th class="n">n</th><th class="n">%&ge;1 bit</th>
<th class="n">other IC</th><th class="n">n</th><th class="n">%&ge;1 bit</th></tr></thead>
<tbody>__SWIC__</tbody></table></div>
<p class="tnote">Switches fall on bases AR reads about as often as anywhere else in the motif,
so these are substitutions at functional positions that phyloP necessarily scores as
unconserved.</p>

<h3>Change in PWM match, with a within-element control</h3>
<div class="tw"><table>
<caption>oscine consensus minus outgroup consensus, log-odds; avian363</caption>
<thead><tr><th>element set</th><th>motif</th><th>class</th><th class="n">AR site</th>
<th class="n">control</th><th class="n">paired</th><th class="n">p</th><th class="n">n</th></tr></thead>
<tbody>__PWM__</tbody></table></div>
<p class="tnote">Control windows tile the rest of each element at the same width, excluding
anything overlapping an AR hit (median 33 windows per site). They sit at zero in every cell,
which confirms that the large positive site values are ascertainment: these sites were called
by FIMO on an oscine reference, so the oscine sequence sits at a selected optimum and the
outgroup sequence does not. The control establishes the bias rather than removing it. The symmetric-ascertainment
analysis below measures and corrects it; within this table the interpretable quantity remains
the Glut-RA minus GABA-4 difference, positive in all four cells and significant in none.</p>

<h3>Symmetric ascertainment: measuring the winner's curse</h3>
<p>The within-element control shows that the oscine advantage is confined to the called sites
and is absent at unselected positions in the same sequence. It does not correct the bias,
because the control windows were never ascertained &mdash; their value near zero is the
unselected baseline, so subtracting it leaves the site value essentially unchanged.</p>
<p>The correction requires ascertaining symmetrically. FIMO selected positions where the
<i>oscine</i> sequence matches the PWM, and any such position sits at a local maximum for that
clade, so the other clade can only regress downward whether or not anything evolved. Scanning
the same elements for positions where the <i>outgroup</i> consensus matches reverses the
selection. If the effect is purely ascertainment the two are equal and opposite; the midpoint
is the corrected estimate and the asymmetry is what remains. Thresholds are set at the same
percentile of each motif's own score distribution, because a fixed log-odds cut is far more
stringent for the 16 bp JASPAR matrix than for the 12 bp ARE-full one.</p>
<div class="tw"><table>
<caption>DHX35&ndash;TOP1 domain, matched cuts at the 99.5th percentile; CIs bootstrap over elements</caption>
<thead><tr><th>motif</th><th>class</th><th class="n">ascert. oscine</th>
<th class="n">ascert. outgroup</th><th class="n">midpoint</th><th class="n">95% CI</th>
<th class="n">n windows</th></tr></thead>
<tbody>__ASC__</tbody></table></div>
<p class="tnote">Forward and reverse are close to equal and opposite, so roughly 70&ndash;85% of
the raw advantage is winner's curse. The GABA-4 ARE-full midpoint lands at
<span class="mono">&minus;0.03</span> &mdash; the class where no directional change is predicted
shows none, which is the control behaving correctly. Glut-RA ARE-full is
<span class="mono">+1.11</span>.</p>

<h3>Threshold sensitivity</h3>
<div class="tw"><table>
<caption>corrected midpoints across matched percentile cuts; shaded rows have a CI excluding zero</caption>
<thead><tr><th>percentile</th><th>motif</th><th class="n">cut</th><th class="n">Glut-RA</th>
<th class="n">GABA-4</th><th class="n">difference</th><th class="n">95% CI</th>
<th class="n">n RA &middot; G4</th></tr></thead>
<tbody>__SWEEP__</tbody></table></div>
<p class="tnote">The Glut-RA midpoint rises monotonically with threshold, consistent with the
signal being concentrated in genuine matches. But the class difference clears zero at one of
five cuts, at a threshold chosen after seeing the data, so it should not be reported as a
result. What the sweep supports is that the corrected Glut-RA value is positive at every cut
and exceeds GABA-4 at four of five, with GABA-4 at or near zero wherever the cut is strict.</p>

<h3>Where in the distribution the signal lives</h3>
__FIG_DIST__
<div class="tw"><table>
<caption>oscine minus outgroup consensus score at matched quantiles</caption>
<thead><tr><th>motif</th><th class="n">quantile</th><th class="n">Glut-RA</th>
<th class="n">GABA-4</th><th class="n">difference</th></tr></thead>
<tbody>__QUANT__</tbody></table></div>
<p class="tnote">Through the bulk of the distribution the two clades match the AR matrices
equally well; medians differ by 0.02&ndash;0.07 log-odds. Separation appears only past the 99th
percentile, and the element-bootstrap bands overlap across the whole range, so the Glut-RA
excess visible in the ARE-full panel is not resolved. This is why the correction can only be estimated in the top fraction of a percent
of windows, where the counts fall to tens &mdash; the restriction is forced by the data rather
than chosen. The P99.99 values rest on two or three windows each and are not usable.</p>

<div class="verdict">
<h3><span class="chip fall">underpowered</span>Both switch metrics</h3>
<p>The switch-rate asymmetry and the PWM-score difference point the same way as the constraint
analysis, and neither resolves at these sample sizes &mdash; 9 to 25 sites per class inside the
MAFB sets, 3 to 11 inside the domain. These analyses are better read as candidate generation
than as tests.</p>
</div>

__FIG_TREE__
</section>

<section class="sect">
<div class="secthead"><span class="no">3</span><div>
  <h2>AR against the other predicted MAFB inputs</h2>
  <p class="eyebrow">cross-motif comparison</p></div></div>

__PART3__

</section>

<section class="sect">
<div class="secthead"><span class="no">4</span><div>
  <h2>AR site gain</h2>
  <p class="eyebrow">presence per lineage, not conservation</p></div></div>

<div class="defn">
<p><b>Default span.</b> Everything below is the gene-bounded
<b>DHX35&ndash;TOP1 domain</b>, <span class="mono">chr20:12,551,036&ndash;13,547,533</span>,
996 kb of intergenic sequence holding MAFB and two putative lncRNAs; a DAR belongs if it
overlaps the interval. The gene &plusmn;1.5 Mb window appears beside it where the comparison is
useful, labelled as such. It was chosen for statistical power rather than biology, and since
the AR effect grows with window width, results stated on it are not MAFB-locus claims.</p>
</div>

__PART4__

<div class="verdict">
<h3><span class="chip cand">candidate</span>Glut-RA.395131</h3>
<p>One element carries converging evidence rather than a test. SCENIC+ links
<span class="mono">chr20:13,241,507&ndash;13,242,007</span> to MAFB through NR3C2, which is
motif-indistinguishable from AR; ArchR links the same element independently; it is a Glut-RA
DAR 124 kb from the gene; it carries a JASPAR AR site with pan-avian phyloP
<span class="mono">+1.95</span>; and position 5 of that site is an oscine-specific
T&rarr;C, carried by 94% of 146 oscines against 92% T in 217 outgroup species.</p>
<p>The tree below shows that column mapped onto leaves and reconstructed internal nodes. The
oscine-stem gain is real but is one of 17 substitutions at this column across the 363-species
tree, with two further independent gains in <i>Cariama</i> and <i>Acanthisitta</i> and a
reversal in <i>Erithacus</i>. This is a site worth testing functionally, not evidence of
lineage-specific selection.</p>
</div>
</section>

<div class="foot">
<p><b>Alignments.</b> avian363 (<code>363-avian-2020_lsd.hal</code>, 364 extant + 361 ancestral
genomes), which is the alignment the phyloP tracks were computed on. Internal-node states are
read from cactus reconstructions, not inferred here. The 74-way
(<code>lonStrDom2/birds.hal</code>, 73 extant) is shown only where noted.</p>
<p><b>Motifs.</b> ARE-full and ARE-half are HOMER <i>de novo</i> motifs 2 and 4 from the
Glut-RA versus background scan; JASPAR AR is MA0007.3. Sites are FIMO hits at p &lt; 1e-4
within DARs, merged across strands.</p>
<p><b>Code.</b> <code>ar_core.py</code> &middot; <code>ar_highic.py</code> &middot;
<code>ar_distance.py</code> &middot; <code>maf_switch/sites363.py</code> &middot;
<code>maf_switch/perm_asym.py</code> &middot; <code>maf_switch/pwm_delta.py</code> &middot;
<code>maf_switch/pwm_control.py</code> &middot; <code>tf_core.py</code> &middot;
<code>allmotif_outlier.py</code> &middot; <code>maf_switch/ascertain_sym2.py</code> &middot;
<code>maf_switch/site_gain.py</code> &middot; <code>maf_switch/gain_calib_radius.py</code> &middot;
<code>maf_switch/presence.py</code> (the presence call) &middot; <code>maf_switch/rebuild_gain.py</code> &middot; <code>maf_switch/gain_trees.py</code> &middot;
<code>maf_switch/rate_test.py</code> &middot; <code>maf_switch/calib_proj.py</code> &middot; <code>maf_switch/plot_site_alignment.py</code>. Every figure and table here regenerates from
<code>consolidate_ar.py</code> and <code>consolidate_ar2.py</code>.</p>
</div>
</div>
"""
rep={"__PART3__":part3(),"__MOTIFS__":t_motifs(),"__TFD__":t_tf("subtree222_domain",[]),"__PART4__":part4(),"__DENS__":t_dens(),"__IC__":t_ic(),"__GWO__":t_gw("subtree222"),
     "__GWW__":t_gw("wholetree"),"__DIST__":t_dist(),"__DOMIC__":t_domic(),
     "__SW__":t_switch(),"__SWPH__":t_swph(),"__SWIC__":t_swic(),"__PWM__":t_pwm(),
     "__TFO__":t_tf("subtree222_gw",["ra"]),"__TFW__":t_tf("wholetree_gw",["ra"]),
     "__TFM__":t_tf("subtree222_MAFB",[]),"__OUT__":t_out(),
     "__FITR__":f'{fit["r"]:+.2f}'.replace("-","−"),
     "__FITP__":f'{fit["p"]:.0e}',"__FITT__":f'{fit["tau"]:.4f}',
     "__FIG_FOLD__":img(f"{F}/results/compare_gw_mafb1p5mb/fig_ar3_fold.png","Motif versus enhancer levels for three AR motifs",
        "<b>Motif and enhancer phyloP, three AR motifs.</b> Open circles are the rest of the "
        "enhancer, filled circles the motif, the arrow the gap between them, solid where "
        "p &lt; 0.05. Merged-site unit. Genome-wide rows use all DARs; MAFB rows use the "
        "published ±1.5 Mb set."),
     "__FIG_SW__":img(f"{M}/../results/mafb_domain/fig_switches363_domain.png","Switch heatmap restricted to the domain",
        "<b>Oscine/outgroup switches inside the DHX35–TOP1 domain.</b> One row per AR site, "
        "one cell per motif column, shaded by how much the oscine major allele gains over the "
        "outgroup. Labelled cells meet the switch criteria. Row labels are red for Glut-RA DARs "
        "and teal for GABA-4; pale grey cells have too few aligned genomes."),
     "__ASC__":t_asc(),"__SWEEP__":t_sweep(),"__QUANT__":t_quant(),
     "__FIG_DIST__":img(f"{M}/../results/mafb_domain/fig_score_dists.png","Motif score distributions",
        "<b>Consensus motif-score distributions.</b> Top: full density for oscine (solid) and "
        "outgroup (dashed) consensus windows in Glut-RA (red) and GABA-4 (teal) DARs \u2014 the "
        "four series are superimposed through the bulk. Bottom: the difference between the "
        "oscine and outgroup score at the same quantile, across the upper tail, with 95% "
        "bands bootstrapped over elements (1,000 resamples). The curves hug zero until about "
        "the 99th percentile and separate only beyond it, where the bands widen enough to "
        "overlap throughout."),
     "__GAIN__":t_gain(),"__SWEEP2__":t_sweep2(),"__PLACE__":t_place(),"__RATE363__":t_rate363(),
     "__FIG_SITEALN__":img(f"{F}/results/per_site/site_alignments_363/mafb_domain/aln363_jaspar_RA_Glut-HVC-2.394967_1.png",
        "Aligned sequence at the AR site gain",
        "<b>The oscine-stem gain at <span class='mono'>Glut-HVC-2.394967</span>, base by base.</b> "
        "Nineteen species from eight clades; dots match lonStrDom2, letters differ, and the logo "
        "is the JASPAR AR matrix in reference orientation. Rows a and b are cactus's "
        "reconstructed Eupasseres and oscine ancestors: two substitutions on the branch between "
        "them, both at high-information positions, create the site. Grey, dotted branches have no "
        "orthologous sequence and carry no presence state."),
     "__FIG_GTREE__":img(f"{F}/results/per_site/gain_trees_ar/mafb_domain/jaspar_RA_Glut-HVC-2.394967_1.png",
        "AR site gain mapped onto the tree",
        "<b>A single AR-site gain on the songbird stem.</b> Presence is read from each node's "
        "own sequence, cactus ancestral reconstructions included, so gains and losses are "
        "observed on branches rather than inferred from tips. Squares are leaves, circles "
        "internal nodes, red present and grey absent. This site is absent from all 56 outgroup "
        "genomes including <i>Acanthisitta</i> and the suboscine <i>Chiroxiphia</i>, gained once "
        "on the oscine stem, and subsequently lost twice within songbirds."),
     "__FIG_TREE__":img(f"{F}/results/per_site/switch_trees_363/mafb_domain/jaspar_RA_Glut-RA.395131_1_p5.png",
        "Target column mapped onto the species tree",
        "<b>The Glut-RA.395131 switch column on the 363-species tree.</b> Squares are leaves, "
        "circles reconstructed internal nodes, orange branches substitutions. Monomorphic "
        "oscine-free clades are collapsed to wedges. Bold names are oscines."),
     "__FIG_TF__":img(f"{F}/results/genome_wide/fig_tf13.png","Cell-type interaction for the thirteen MAFB regulators",
        "<b>Cell-type interaction across the SCENIC+ MAFB regulators.</b> Glut-RA minus GABA-4 "
        "for each TF's paired motif-versus-enhancer value, with 95% bootstrap CIs. Asterisks "
        "mark FDR 5% within each panel. Rows are ordered by the oscine value. Every regulator, "
        "AR included, uses the SCENIC+-derived matrix in the table above."),
     "__FIG_OUT__":img(f"{F}/results/compare_gw_mafb1p5mb/fig_allmotif_outlier.png","Regulators against the archetype trend",
        "<b>The regulators against 182 motif archetypes.</b> Grey points are archetypes, the "
        "line their OLS fit, the band ±2 robust SD of their residuals. AR (diamond) is the "
        "SCENIC+ matrix, shared with NR3C2; AR* is the HOMER ARE-full motif, a reference."),
     }
html=HTML
for k,v in rep.items(): html=html.replace(k,v)
open(OUT,"w").write(html)
print(f"  wrote {OUT}  ({len(html)/1e6:.2f} MB)")
