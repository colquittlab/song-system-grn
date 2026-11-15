#!/usr/bin/env python
"""Browser BED9 tracks for the constraint clusters C1 (deep-conserved) and C2 (oscine-lineage constrained).

Reads results/mafb_domain/element_clusters_k4.tsv (plot_element_heatmaps.R) and element_p2g_mafb.tsv
(plot_element_scatter.R). One line per lobe: chromStart-chromEnd = the SnapATAC2 peak, thickStart-
thickEnd = the scored summit +/-100 bp, so multi-part peaks show one line per lobe with the lobe's core
drawn thick. Name = class|cluster|whole % > 1|oscine % > 1|oscine block %|MAFB p2g r[|lobe i/n].

    python element_cluster_bed.py
Writes results/mafb_domain/element_clusters_{C1_C2,C1,C2}.bed.
"""
import pandas as pd

R = "/private/warm-archive/colquittlab/assembly/cactus/avian363/phyloP_matched/loci_focus/results/mafb_domain"
LAB = {1: "C1_deep", 2: "C2_oscine"}
COL = {1: "200,22,26", 2: "27,175,122"}
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
        rows.append((r.chrom, r.peak_start, r.peak_end, nm, 0, ".", r.summit - CORE, r.summit + CORE, COL[r.cluster], r.cluster))
    B = pd.DataFrame(rows)
    for keep, tag in (((1, 2), "C1_C2"), ((1,), "C1"), ((2,), "C2")):
        x = B[B[9].isin(keep)].iloc[:, :9]
        with open(f"{R}/element_clusters_{tag}.bed", "w") as fh:
            fh.write(f'track name="MAFB_{tag}" description="MAFB domain RA/G4 elements, constraint cluster '
                     f'{tag.replace("_", " + ")} (thick = summit +/-{CORE} bp)" itemRgb="On"\n')
            x.to_csv(fh, sep="\t", header=False, index=False)
        sub = d[d.cluster.isin(keep)]
        print(f"  element_clusters_{tag}.bed: {len(x)} lobes (RA {(sub.cls == 'RA').sum()}, G4 {(sub.cls == 'G4').sum()}) "
              f"from {sub.drop_duplicates(['cls', 'peak_start']).shape[0]} peaks")


if __name__ == "__main__":
    main()
