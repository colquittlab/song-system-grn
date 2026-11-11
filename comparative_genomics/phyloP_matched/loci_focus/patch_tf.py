import re,base64,sys
R="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/report/song_enhancer_report.html"
h=open(R).read(); n0=len(h)
def sub(old,new,label):
    global h
    if h.count(old)!=1:
        sys.exit("ANCHOR %s: %d matches"%(label,h.count(old)))
    h=h.replace(old,new)

def tbody(caption,path,label):
    global h
    i=h.index(caption); j=h.index("<tbody>",i)+len("<tbody>"); k=h.index("</tbody>",j)
    h=h[:j]+"\n"+open(path).read()+"\n    "+h[k:]

# ---- 1. earlier callout -------------------------------------------------
sub("""genome-wide, <span class="mono">+0.201, p = 0.0033</span> at MAFB)""",
    """genome-wide, <span class="mono">+0.201, p = 0.0067</span> at MAFB)""","callout-p")
sub("""The reversal has since been tested against the thirteen SCENIC+ MAFB
    regulators and is <em>not</em> unique to AR &mdash; KLF12 and KLF15 show it more strongly.""",
    """The reversal has since been tested against the other twelve SCENIC+ MAFB
    regulators and is <em>not</em> unique to AR &mdash; the KLF family shows it more strongly.""","callout-klf")

# ---- 2. section intro ---------------------------------------------------
sub("""<p>A later SCENIC+ run nominates seven further inputs to MAFB in Glut-RA (KLF12, KLF15, NR3C2,
  NFIA, LIN28B, STAT1, FOSL2) and six in GABA-4 (MAF, SOX6, TCF12, ARX, SATB1, MEIS2). Running
  each through the paired statistic used for AR asks whether the regulon split leaves a mark on
  constraint, and gives AR a purpose-built comparison set rather than the 182 arbitrary
  archetypes.</p>""",
    """<p>A later SCENIC+ run nominates six further inputs to MAFB in Glut-RA (KLF, NR3C2, NFIA,
  LIN28B, STAT1, FOSL2) and six in GABA-4 (MAF, SOX6, TCF12, ARX, SATB1, MEIS2). AR belongs to
  the Glut-RA regulon and is treated here as one of its thirteen members rather than as an
  outside reference, so it sits inside the multiple-testing correction with the rest. Running
  each through the paired statistic asks whether the regulon split leaves a mark on constraint,
  and gives AR a purpose-built comparison set rather than the 182 arbitrary archetypes.</p>""","intro")

# ---- 3. finding: oscine ------------------------------------------------
sub("""  <div class="finding falls">
    <div class="vhead"><span class="chip">no oscine signal</span>
      <p class="claim">None of the thirteen regulators shows a Glut-RA / GABA-4 constraint
      asymmetry on the oscine lineage once multiple testing is accounted for.</p></div>
    <p class="evidence">Genome-wide, across 1,207&ndash;6,863 elements per motif per class,
    <span class="mono">0 of 13</span> survive FDR 5% on the oscine subtree. Three are nominally
    significant &mdash; TCF12 <span class="mono">+0.022</span>, ARX
    <span class="mono">&minus;0.019</span>, MAF <span class="mono">&minus;0.017</span> &mdash;
    and TCF12 runs <em>opposite</em> to its GABA-4 assignment. The sign matches the SCENIC+
    regulon in <span class="mono">7 of 13</span> cases, which is chance. AR's
    <span class="mono">+0.023</span> would rank <span class="mono">1 of 14</span>, but against a
    distribution this flat (median <span class="mono">+0.004</span>, sd
    <span class="mono">0.013</span>) that is a weak distinction.</p>
  </div>""",
    """  <div class="finding holds">
    <div class="vhead"><span class="chip">AR alone on the oscine lineage</span>
      <p class="claim">AR is the only one of the thirteen regulators whose Glut-RA / GABA-4
      constraint asymmetry on the oscine lineage survives correction.</p></div>
    <p class="evidence">Genome-wide, across 1,207&ndash;8,927 elements per motif per class,
    <span class="mono">1 of 13</span> clears FDR 5% on the oscine subtree: AR, at
    <span class="mono">+0.023</span> (<span class="mono">p = 0.0025, q = 0.033</span>), the
    largest interaction in the set. Three others are nominally significant and none of them
    survives &mdash; TCF12 <span class="mono">+0.022</span>, ARX
    <span class="mono">&minus;0.019</span>, MAF <span class="mono">&minus;0.017</span> &mdash;
    and TCF12 runs <em>opposite</em> to its GABA-4 assignment. The sign matches the SCENIC+
    regulon in <span class="mono">7 of 13</span> cases, which is chance.</p>
    <p class="evidence">The distinction is real but narrow. The thirteen values span a band only
    0.042 wide (median <span class="mono">+0.004</span>, sd <span class="mono">0.015</span>),
    so AR leads a flat field: it is 1.5 standard deviations above the median, not an outlier of
    a different kind. Its advantage comes as much from sample size &mdash; 8,927 Glut-RA
    elements, the largest here &mdash; as from effect size, and the test is better powered for
    AR than for any other motif in the set. Ranking first among thirteen pre-specified
    regulators is nonetheless a stronger statement than ranking seventh among 182 arbitrary
    archetypes.</p>
  </div>""","finding-oscine")

# ---- 4. finding: whole tree --------------------------------------------
sub("""    direction tracks the regulon: the Glut-RA set averages <span class="mono">+0.025</span>
    against <span class="mono">&minus;0.069</span> for the GABA-4 set. KLF12
    (<span class="mono">+0.119</span>), KLF15 (<span class="mono">+0.092</span>) and NFIA
    (<span class="mono">+0.069</span>) are Glut-RA-biased as predicted; TCF12""",
    """    direction tracks the regulon: the Glut-RA set averages <span class="mono">+0.024</span>
    against <span class="mono">&minus;0.069</span> for the GABA-4 set. KLF
    (<span class="mono">+0.119</span>), AR (<span class="mono">+0.086</span>) and NFIA
    (<span class="mono">+0.069</span>) are Glut-RA-biased as predicted; TCF12""","finding-tree-a")
sub("""clear exception, GABA-4-assigned but <span class="mono">+0.139</span>. AR ranks
    <span class="mono">4 of 14</span> here.</p>""",
    """clear exception, GABA-4-assigned but <span class="mono">+0.139</span> &mdash; the largest
    value in the set and in the wrong direction. AR ranks <span class="mono">3 of 13</span>
    here, and the sign matches the regulon in <span class="mono">8 of 13</span> cases.</p>""","finding-tree-b")

# ---- 5/6. table 1 and its note -----------------------------------------
tbody("<caption>motif minus its own enhancer, genome-wide, by SCENIC+ regulon</caption>",
      "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/t_dod.html","t1")
sub("""  <p class="tnote">Positive = the motif is held more tightly, relative to its own enhancer, in
  Glut-RA than in GABA-4 DARs. p from 20,000 permutations; the significance marks are
  BH-corrected within each tree. Nine of the thirteen PWMs are family-level metaclusters, so
  "FOSL2" is AP-1, "MAF" the MAF family, "KLF12" the KLF family and so on &mdash; KLF12 and KLF15
  in particular are near-duplicate measurements. LIN28B is an RNA-binding protein with no
  validated DNA motif; its CIS-BP matrix is inferred and the row should be discounted. At MAFB
  itself no TF reaches significance (6&ndash;48 elements per class, all p &gt; 0.08).</p>""",
    """  <p class="tnote">Values are BH-corrected q, from 20,000 label permutations within each tree
  across all thirteen regulators including AR. Positive = the motif is held more tightly,
  relative to its own enhancer, in Glut-RA than in GABA-4 DARs. Eight of the thirteen PWMs are
  family-level metaclusters, so "FOSL2" is AP-1, "MAF" the MAF family, "KLF" the KLF family and
  so on; KLF12 and KLF15 were nominated separately by SCENIC+ but map to one metacluster and are
  reported once. LIN28B is an RNA-binding protein with no validated DNA motif; its CIS-BP matrix
  is inferred and the row should be discounted. At MAFB itself no TF reaches significance
  (8&ndash;48 elements per class, all p &gt; 0.09).</p>""","t1note")

# ---- 7. AR mid-pack finding --------------------------------------------
sub("""    <p class="evidence">In Glut-RA DARs on the oscine subtree, <span class="mono">11 of 13</span>
    regulators are positive and almost all clear FDR 5%. AR's <span class="mono">+0.044</span>
    ranks <span class="mono">6 of 14</span>, well behind FOSL2
    (<span class="mono">+0.131</span>), MAF (<span class="mono">+0.086</span>), NFIA
    (<span class="mono">+0.057</span>), MEIS2 (<span class="mono">+0.055</span>) and SOX6
    (<span class="mono">+0.051</span>). Pan-avian it ranks <span class="mono">11 of 14</span>.""",
    """    <p class="evidence">In Glut-RA DARs on the oscine subtree, <span class="mono">11 of 13</span>
    regulators are positive and almost all clear FDR 5%. AR's <span class="mono">+0.044</span>
    ranks <span class="mono">6 of 13</span>, well behind FOSL2
    (<span class="mono">+0.131</span>), MAF (<span class="mono">+0.086</span>), NFIA
    (<span class="mono">+0.057</span>), MEIS2 (<span class="mono">+0.055</span>) and SOX6
    (<span class="mono">+0.051</span>). Pan-avian it ranks <span class="mono">11 of 13</span>.""","midpack")

# ---- 8/9. table 2 and note ---------------------------------------------
tbody("<caption>motif minus the rest of its own enhancer, Glut-RA DARs</caption>",
      "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/t_ra.html","t2")
sub("""  <p class="tnote">Positive = the motif is more constrained than the remainder of the same
  element. BH-corrected within each tree across the thirteen regulators; AR is shown for
  reference and was not part of that correction.</p>""",
    """  <p class="tnote">Positive = the motif is more constrained than the remainder of the same
  element. BH-corrected q within each tree across all thirteen regulators, AR included. This is
  the main-effect column; the interaction in the table above is the difference between it and
  the matching GABA-4 column.</p>""","t2note")

# ---- 10. KLF reversal finding ------------------------------------------
sub("""      <p class="claim">The sign reversal between trees is not unique to AR &mdash; the KLF
      family does it more strongly.</p></div>
    <p class="evidence">Earlier this section offered the tree-to-tree reversal as one of the
    properties setting AR apart. Of the thirteen regulators, two others reverse the same way in
    Glut-RA DARs, and by a wider margin: KLF12 swings from
    <span class="mono">&minus;0.186</span> pan-avian to <span class="mono">+0.038</span> in
    oscines and KLF15 from <span class="mono">&minus;0.167</span> to
    <span class="mono">+0.038</span>, against AR's <span class="mono">&minus;0.033</span> to
    <span class="mono">+0.044</span>. KLF12 and KLF15 are near-duplicate measurements, so they
    count as roughly one observation &mdash; but the pattern is clearly not AR's alone. What AR
    retains is that it is the only motif here to combine the reversal with a positive cell-type
    interaction; KLF12 and KLF15 end up at <span class="mono">+0.014</span> and
    <span class="mono">+0.004</span> on that measure.</p>""",
    """      <p class="claim">The sign reversal between trees is not unique to AR &mdash; the KLF
      family does it more strongly.</p></div>
    <p class="evidence">Earlier this section offered the tree-to-tree reversal as one of the
    properties setting AR apart. Of the thirteen regulators, two reverse the same way in Glut-RA
    DARs and the KLF family does it by a wider margin, swinging from
    <span class="mono">&minus;0.186</span> pan-avian to <span class="mono">+0.038</span> in
    oscines against AR's <span class="mono">&minus;0.033</span> to
    <span class="mono">+0.044</span>. LIN28B reverses too but is negative in both and should be
    discounted. So the reversal is not AR's alone. What AR retains is that it is the only motif
    to combine the reversal with a cell-type interaction that survives correction; KLF reaches
    <span class="mono">+0.014</span> on that measure (<span class="mono">q = 0.44</span>).</p>""","klfrev")

# ---- 11. fig_tf13 caption ----------------------------------------------
sub("""    Each TF's motif scored by the paired statistic used for AR, with the Glut-RA minus GABA-4
    interaction plotted; bars are 95% bootstrap CIs, colour is the SCENIC+ regulon assignment,
    and the dashed line marks AR. Rows are ordered by the oscine value and held in that order in
    both panels. On the oscine subtree (A) the regulators sit in a narrow band around zero and
    nothing survives correction. On the whole tree (B) they spread out and sort largely by
    assignment &mdash; so what the regulon split predicts is ancient constraint, not
    oscine-specific constraint.</figcaption>""",
    """    Each TF's motif scored by the paired statistic, with the Glut-RA minus GABA-4 interaction
    plotted; bars are 95% bootstrap CIs, colour is the SCENIC+ regulon assignment, AR is the
    diamond, and an asterisk marks FDR 5% within the panel. Rows are ordered by the oscine value
    and held in that order in both panels. On the oscine subtree (A) the regulators sit in a
    narrow band around zero and only AR clears correction. On the whole tree (B) they spread out
    an order of magnitude wider and sort largely by assignment &mdash; so most of what the
    regulon split predicts is ancient constraint rather than oscine-specific
    constraint.</figcaption>""","cap1")

# ---- 12. MAFB finding ---------------------------------------------------
sub("""    <p class="evidence">AR scores <span class="mono">+0.201 (p = 0.003)</span> in Glut-RA DARs
    against <span class="mono">&minus;0.025 (p = 0.84)</span> in GABA-4, an interaction of
    <span class="mono">+0.226</span> that ranks <span class="mono">1 of 12</span> &mdash; the
    other eleven span only <span class="mono">&minus;0.144</span> to
    <span class="mono">+0.060</span>, so AR is nearly four times the next largest positive
    value. But it reaches only <span class="mono">p = 0.094</span> on 26 and 29 elements, and
    <span class="mono">0 of 12</span> motifs clear significance on the interaction at this
    locus.</p>""",
    """    <p class="evidence">AR scores <span class="mono">+0.201 (p = 0.0065)</span> in Glut-RA DARs
    against <span class="mono">&minus;0.025 (p = 0.83)</span> in GABA-4, an interaction of
    <span class="mono">+0.226</span> that ranks <span class="mono">1 of 11</span> &mdash; the
    other ten span only <span class="mono">&minus;0.144</span> to
    <span class="mono">+0.060</span>, so AR is nearly four times the next largest positive
    value. But it reaches only <span class="mono">p = 0.098</span> on 26 and 29 elements, and
    <span class="mono">0 of 11</span> motifs clear correction on the interaction at this locus
    &mdash; the smallest q here is <span class="mono">0.97</span>.</p>""","mafbfind")

# ---- 13/14. MAFB table and note ----------------------------------------
tbody("<caption>at MAFB, oscine subtree: motif minus its own enhancer</caption>",
      "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/t_mafb.html","t3")
sub("""  <p class="tnote">LIN28B (7 / 8 elements) and SATB1 (6 / 8) fall below the minimum and are not
  shown. On the whole tree at MAFB nothing is significant either: ARX has the largest interaction
  (<span class="mono">+1.338, p = 0.083</span>) and AR ranks 5 of 12
  (<span class="mono">+0.325, p = 0.43</span>).</p>""",
    """  <p class="tnote">Uncorrected p from 20,000 permutations; with eleven motifs tested the
  smallest interaction q is <span class="mono">0.97</span>, so no column survives correction at
  this locus. LIN28B (7 / 8 elements) and SATB1 (6 / 8) fall below the eight-element minimum and
  are not shown. On the whole tree at MAFB nothing is significant either: ARX has the largest
  interaction (<span class="mono">+1.338, p = 0.086</span>) and AR ranks 5 of 11
  (<span class="mono">+0.325, p = 0.43</span>).</p>""","t3note")

# ---- 15. MAFB fig caption ----------------------------------------------
sub("""    Glut-RA minus GABA-4 for each regulator's paired motif-vs-enhancer value, oscine subtree,
    with AR as the grey diamond; bars are 95% bootstrap CIs and rows are ordered by the
    genome-wide value. Genome-wide (A) the regulators sit in a band a few hundredths wide and AR
    is barely separated from them. At MAFB (B) the band widens with the smaller samples, but AR
    is the only motif displaced clearly to the right; LIN28B and SATB1 fall below the minimum
    element count there and are blank.</figcaption>""",
    """    Glut-RA minus GABA-4 for each regulator's paired motif-vs-enhancer value, oscine subtree,
    AR as the diamond; bars are 95% bootstrap CIs, an asterisk marks FDR 5% and rows are ordered
    by the genome-wide value. Genome-wide (A) the regulators sit in a band a few hundredths wide
    and AR sits at its right edge, the only one to clear correction. At MAFB (B) the band widens
    ten-fold with the smaller samples and nothing is significant, but AR is the only motif
    displaced clearly to the right; LIN28B and SATB1 fall below the minimum element count there
    and are blank.</figcaption>""","cap2")

# ---- 16. re-embed both figures -----------------------------------------
F="/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/"
for alt,png in (("SCENIC+ MAFB regulators tested with the paired statistic","results/genome_wide/fig_tf13.png"),
                ("SCENIC+ regulators genome-wide versus at MAFB","results/compare_gw_mafb1p5mb/fig_tf13_mafb.png")):
    b=base64.b64encode(open(F+png,"rb").read()).decode()
    pat=re.compile(r'(<img alt="%s" src="data:image/png;base64,)[A-Za-z0-9+/=]+(")'%re.escape(alt))
    h,k=pat.subn(lambda m:m.group(1)+b+m.group(2),h)
    if k!=1: sys.exit("img %s: %d"%(png,k))

open(R,"w").write(h)
print("patched  %d -> %d bytes"%(n0,len(h)))
