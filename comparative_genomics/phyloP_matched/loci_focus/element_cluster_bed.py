#!/usr/bin/env python
"""Browser BED9 tracks for the MAFB-domain accessible elements.

Reads results/mafb_domain/element_clusters_k4.tsv (plot_element_heatmaps.R; lobes with Ward
cluster, rule class and scores) and element_p2g_mafb.tsv (plot_element_scatter.R). One line per lobe.
Rule-class tracks (element_rule_classes.py) draw each lobe's full scored extent; cluster and class
tracks draw the SnapATAC2 peak with the summit +/-100 bp thick (clipped to the peak). Class is RA,
G4 or both (an RA and a G4 peak overlapping each other, or a peak overlapping DARs of both contrasts).

    python element_cluster_bed.py
Writes results/mafb_domain/
  element_rules.bed and element_rule_{deep_oscine,deep,oscine,weak}.bed   rule classes
  element_classes.bed                                                    all lobes, coloured by class
  element_clusters_{C1_C2,C1,C2}.bed                                     Ward clusters C1 and C2
"""
import pandas as pd

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
LAB = {1: "C1_deep", 2: "C2_oscine"}
COL = {1: "200,22,26", 2: "27,175,122"}
CLS_COL = {"RA": "235,104,52", "G4": "42,120,214", "both": "120,60,160"}
RULE_COL = {"deep+oscine": "120,60,160", "deep": "200,22,26", "oscine": "27,175,122", "weak": "150,150,150"}
CORE = 100


def main():
    d = pd.read_csv(f"{R}/element_clusters_k4.tsv", sep="\t").merge(
        pd.read_csv(f"{R}/element_p2g_mafb.tsv", sep="\t")[["element", "p2g_r"]], on="element")
    d = d[d.cluster.isin(LAB)].sort_values(["chrom", "peak_start", "lobe"])
    rows = []
    for r in d.itertuples(index=False):
        p2g = "NA" if pd.isna(r.p2g_r) else f"{r.p2g_r:.2f}"
        nm = (f"{r.cls}|{LAB[r.cluster]}|whole%>1={r.whole_pct1:.0f}|osc%>1={r.osc_pct1:.0f}|oscblock={r.osc_block:.0f}|p2g={p2g}"
              + (f"|lobe{r.lobe}/{r.n_lobes}" if r.n_lobes > 1 else ""))
        rows.append((r.chrom, r.peak_start, r.peak_end, nm, 0, ".", max(r.peak_start, r.summit - CORE), min(r.peak_end, r.summit + CORE), COL[r.cluster], r.cluster))
    B = pd.DataFrame(rows)
    # all elements, coloured by class
    a = pd.read_csv(f"{R}/element_clusters_k4.tsv", sep="\t").merge(
        pd.read_csv(f"{R}/element_p2g_mafb.tsv", sep="\t")[["element", "p2g_r"]], on="element").sort_values(["chrom", "peak_start", "lobe"])
    with open(f"{R}/element_classes.bed", "w") as fh:
        fh.write(f'track name="MAFB_elem_class" description="MAFB domain elements by class: RA, G4, both (thick = summit +/-{CORE} bp, clipped to peak)" itemRgb="On"\n')
        for r in a.itertuples(index=False):
            p2g = "NA" if pd.isna(r.p2g_r) else f"{r.p2g_r:.2f}"
            nm = (f"{r.cls}|C{r.cluster}|{r.rule}|whole%>1={r.whole_pct1:.0f}|osc%>1={r.osc_pct1:.0f}|p2g={p2g}"
                  + (f"|lobe{r.lobe}/{r.n_lobes}" if r.n_lobes > 1 else ""))
            fh.write("\t".join(map(str, (r.chrom, r.peak_start, r.peak_end, nm, 0, ".", max(r.peak_start, r.summit - CORE),
                                          min(r.peak_end, r.summit + CORE), CLS_COL[r.cls]))) + "\n")
    # rule classes: full extent, coloured by rule
    with open(f"{R}/element_rules.bed", "w") as fh0:
        fh0.write(f'track name="MAFB_elem_rule" description="MAFB domain elements by constraint rule class (full scored extent)" itemRgb="On"\n')
        for rc, col in RULE_COL.items():
            x = a[a.rule == rc]; tag = rc.replace("+", "_")
            with open(f"{R}/element_rule_{tag}.bed", "w") as fh:
                fh.write(f'track name="MAFB_{tag}" description="MAFB domain elements, rule class {rc} (full scored extent)" itemRgb="On"\n')
                for r in x.itertuples(index=False):
                    p2g = "NA" if pd.isna(r.p2g_r) else f"{r.p2g_r:.2f}"
                    nm = (f"{r.cls}|{rc}|whole_block={r.whole_ext_block:.0f}|osc%>1={r.osc_ext_pct1:.0f}|p2g={p2g}"
                          + (f"|lobe{r.lobe}/{r.n_lobes}" if r.n_lobes > 1 else ""))
                    line = "\t".join(map(str, (r.chrom, r.ext_start, r.ext_end, nm, 0, ".", r.ext_start, r.ext_end, col))) + "\n"
                    fh.write(line); fh0.write(line)
            print(f"  element_rule_{tag}.bed: {len(x)} lobes (" + ", ".join(f"{k} {(x.cls == k).sum()}" for k in CLS_COL) + ")")
    print(f"  element_classes.bed: {len(a)} lobes (" + ", ".join(f"{k} {(a.cls == k).sum()}" for k in CLS_COL) + ")")
    for keep, tag in (((1, 2), "C1_C2"), ((1,), "C1"), ((2,), "C2")):
        x = B[B[9].isin(keep)].iloc[:, :9]
        with open(f"{R}/element_clusters_{tag}.bed", "w") as fh:
            fh.write(f'track name="MAFB_{tag}" description="MAFB domain RA/G4 elements, constraint cluster '
                     f'{tag.replace("_", " + ")} (thick = summit +/-{CORE} bp, clipped to peak)" itemRgb="On"\n')
            x.to_csv(fh, sep="\t", header=False, index=False)
        sub = d[d.cluster.isin(keep)]
        print(f"  element_clusters_{tag}.bed: {len(x)} lobes (" + ", ".join(f"{k} {(sub.cls == k).sum()}" for k in ("RA", "G4", "both")) + ") "
              f"from {sub.drop_duplicates(['cls', 'peak_start']).shape[0]} peaks")


if __name__ == "__main__":
    main()
