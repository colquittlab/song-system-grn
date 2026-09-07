"""Class-level accuracy benchmark for cross-species matches.

The earlier gold set relied on marker-named reference labels (Sst, Pvalb, VLMC ...) and
so only covered 15 mostly non-neuronal clusters. This one instead asks a coarser but
far more general question that covers EVERY finch cluster: does the matched reference
label belong to the right broad cell class?

Finch cluster name -> expected reference class, from the finch naming scheme itself:
  PC-*            progenitor / radial glia
  *-NB, *-IP*     neuroblast / intermediate progenitor
  GABA-*, Glut-*  neuron
  Astro           astrocyte lineage
  OPC, Oligo      oligodendrocyte lineage
  Micro           immune
  Epen            ependymal
  ChP             choroid plexus / ependymal
  Endo-*          vascular / mesenchymal

This is a weak test -- a method can get the class right and the subtype badly wrong --
so it measures gross mis-assignment, not fine accuracy. It is nonetheless the only
label-independent ground truth available across all clusters, and gross errors are
exactly what a composite is supposed to suppress.
"""
import re
import numpy as np
import pandas as pd

FINCH_EXPECT = [
    (re.compile(r"^PC-"), "progenitor"),
    (re.compile(r"(-NB(-\d+)?$)|IP"), "neuroblast"),
    (re.compile(r"^(GABA|Glut)-"), "neuron"),
    (re.compile(r"^Astro(-\d+)?$"), "astro"),
    (re.compile(r"^(OPC|Oligo)(-\d+)?$"), "oligo"),
    (re.compile(r"^Micro$"), "immune"),
    (re.compile(r"^Epen$|^Epen-"), "ependymal"),
    (re.compile(r"^ChP$"), "ependymal"),
    (re.compile(r"^Endo(-\d+)?$"), "vascular"),
]

# reference-label -> class, matched against the reference annotation text
REF_PAT = {
    "progenitor": r"Radial glia|\bRgl|Neural progenitor|\bNPC",
    "neuroblast": r"Neuroblast|\bNbl|IMN\b",
    "neuron":     r"Neuron|Gaba|Glut|Chol|Dopa|\bNeur\d|IT |ET |CTX|Sst|Pvalb|Vip|Lamp5|MSN|D1|D2",
    "astro":      r"Astro|Glioblast|\bGbl",
    "oligo":      r"Oligo|\bOPC|COP |MOL|NFOL|Committed oligodendrocyte",
    "immune":     r"Immune|Microglia|\bMgl|\bPvm|BAM |DC NN|Macrophage",
    # NB: must not match Yao's combined class string "30 Astro-Epen", hence no bare
    # "Epen" -- and 'astro' is tested BEFORE 'ependymal' below for the same reason.
    "ependymal":  r"Ependymal|Epen\d|^Epen|Tanycyte|CHOR|Chpl|Choroid|Hypendymal",
    "vascular":   r"Vascular|Endo|VLMC VLMC|VLMC|Peri|Peric|Pia|Meninges|Arachnoid|Dura|Fibro|Vendo|Angiob|ABC NN|Mesenchyme|SMC|Vsm",
}


def expected_class(finch_cluster: str):
    for pat, cls in FINCH_EXPECT:
        if pat.search(finch_cluster):
            return cls
    return None


def ref_class(label: str, annot: pd.DataFrame | None = None):
    """Classify a reference label, preferring the label's OWN name over its annotation.

    Two-stage because the annotation text is sometimes ambiguous by construction: Yao's
    class field for both ependymal and astrocyte labels is the single string
    "30 Astro-Epen", so scoring "1175 Ependymal NN_1" off the annotation gives astro and
    scoring "1166 Astro-OLF NN_1" off it gives ependymal, depending only on test order.
    The label name itself ("Ependymal NN_1", "Astro-OLF NN_1", "CHOR NN_1") is specific,
    so it is tried alone first; the annotation is consulted only when the name is
    uninformative, which is the case for La Manno's numeric Neur###/Nbl###/Rgl### names.
    """
    order = ["immune", "ependymal", "astro", "oligo", "vascular", "progenitor",
             "neuroblast", "neuron"]
    for text in ([str(label)] if annot is None else
                 [str(label),
                  str(label) + " " + " ".join(
                      str(annot.loc[label, c]) for c in annot.columns)
                  if label in annot.index else str(label)]):
        for cls in order:
            if re.search(REF_PAT[cls], text, re.I):
                return cls
    return "unknown"


def score(pred: pd.Series, annot: pd.DataFrame | None = None, verbose=True):
    rows = []
    for c, p in pred.items():
        exp = expected_class(c)
        if exp is None:
            continue
        got = ref_class(p, annot)
        rows.append(dict(cluster=c, expected=exp, got=got, ok=(got == exp), pred=p))
    D = pd.DataFrame(rows)
    if verbose and len(D):
        bad = D[~D.ok]
        if len(bad):
            print(f"    misclassified {len(bad)}/{len(D)}:")
            for _, r in bad.head(12).iterrows():
                print(f"      {r.cluster:20s} expected {r.expected:11s} got {r.got:11s}  {str(r.pred)[:44]}")
    return (D.ok.mean() if len(D) else np.nan), D
