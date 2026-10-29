#!/usr/bin/env python3
"""Effect of the region-to-gene search space on the SCENIC+ results of glut-ra_glut-c1h-1_allruns_hybrid: config2 (150 kb, default), config3 (250 kb),
config4 (500 kb); otherwise identical (see ../scenicplus/make_configs.py).

Per config:
  search space        region-gene pairs considered, by distance band (<= 150 kb, 150-250 kb, 250-500 kb)
  region-to-gene      pairs fitted
  eRegulons           +/+ eRegulons (direct and extended), TFs, unique TF-gene and TF-region links, median targets per eRegulon; and the distance band of
                      the region-gene links they use (joined to the config's region-to-gene table)
  per TF (+/+)        targets and regions, Jaccard overlap of target genes with config2, Pearson r of cell-level gene-based AUC with config2, and the
                      RA vs C1H-1 difference in gene-based AUC (standardized mean difference over cells; cells are not independent samples, so this is
                      an effect size for comparing configs, not a test)
  run time            per Snakemake rule (configs 3 and 4; config2 reused config1's search space and region-to-gene fit)
config2 used config1's search space and region-to-gene table (shared/outs).

Run in the scenicplus environment:
    python search_space_comparison.py
Outputs (next to this script): search_space_summary.csv, search_space_link_distance.csv, search_space_per_tf.csv, search_space_runtime.csv
"""
import re
from datetime import datetime
from pathlib import Path

import anndata
import mudata
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/glut-ra_glut-c1h-1_allruns_hybrid")
RES = R / "results"
CONFIGS = {"config2": 150, "config3": 250, "config4": 500}             # search-space limit, kb
UPSTREAM = {"config2": RES / "shared" / "outs", "config3": RES / "config3" / "outs", "config4": RES / "config4" / "outs"}   # search space + region-to-gene
BANDS = [(0, 150_000, "<=150 kb"), (150_000, 250_000, "150-250 kb"), (250_000, 500_001, "250-500 kb")]
RA, C1H1 = "Glut-CACNA1H-RA", "Glut-CACNA1H-1"


def band(d):
    out = pd.Series(pd.NA, index=d.index, dtype="object")
    for lo, hi, name in BANDS:
        out[(d >= lo) & (d < hi)] = name
    return out


def distance(col):
    return col.astype(str).str.strip("[]").str.split(",").str[0].astype(float).abs()


def read_eregulons(c):
    e = pd.concat([pd.read_csv(RES / c / "outs" / f, sep="\t") for f in ("eRegulon_direct.tsv", "eRegulons_extended.tsv")])
    return e


def runtime(c):
    """Duration of each Snakemake rule from the job's .err log (start of a rule to the start of the next block)."""
    logs = sorted((RES / c / "logs").glob("*.err"))
    if not logs:
        return []
    lines = logs[0].read_text().splitlines()
    stamps, rows = [], []
    for i, l in enumerate(lines):
        m = re.match(r"\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]", l)
        if m:
            t = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
            nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
            rule = re.match(r"(?:local)?rule (\w+):", nxt)
            stamps.append((t, rule.group(1) if rule else None))
    for (t, rule), (t2, _) in zip(stamps, stamps[1:]):
        if rule and rule != "all":
            rows.append({"config": c, "rule": rule, "minutes": round((t2 - t).total_seconds() / 60, 1)})
    if stamps:
        rows.append({"config": c, "rule": "TOTAL", "minutes": round((stamps[-1][0] - stamps[0][0]).total_seconds() / 60, 1)})
    return rows


def gene_auc(c, names):
    """Cell x TF gene-based AUC of the TF's +/+ eRegulon (direct if it exists, else extended)."""
    out = {}
    for kind in ("extended", "direct"):                 # direct overwrites extended
        m = mudata.read(str(RES / c / "outs" / f"AUCell_{kind}.h5mu"))
        g = m.mod["Gene_based"]
        for v in g.var_names:
            mm = re.match(r"(.+)_(direct|extended)_\+/\+_\(\d+g\)$", v)
            if mm and mm.group(1) in names:
                out[mm.group(1)] = pd.Series(np.asarray(g[:, v].X).ravel(), index=g.obs_names)
    return pd.DataFrame(out)


def main():
    labels = anndata.read_h5ad(R / "anndata_rna" / "adata.h5ad", backed="r").obs["cluster"].astype(str)
    summary, dist_rows, per_tf, auc = [], [], {}, {}
    for c, kb in CONFIGS.items():
        ss = pd.read_csv(UPSTREAM[c] / "search_space.tsv", sep="\t", usecols=["Name", "Gene", "Distance"])
        ss["band"] = band(distance(ss["Distance"]))
        r2g = pd.read_csv(UPSTREAM[c] / "region_to_gene_adj.tsv", sep="\t", usecols=["target", "region", "Distance"])
        r2g["band"] = band(distance(r2g["Distance"]))
        e = read_eregulons(c)
        pp = e[e["eRegulon_name"].str.endswith("+/+")]
        links = pp[["Gene", "Region"]].drop_duplicates().merge(r2g.rename(columns={"target": "Gene", "region": "Region"}), on=["Gene", "Region"], how="left")
        row = {"config": c, "search_space_kb": kb,
               "search_space_pairs": len(ss), **{f"search_space_{b}": int((ss["band"] == b).sum()) for _, _, b in BANDS},
               "region_to_gene_pairs": len(r2g),
               "eRegulons_all": e["eRegulon_name"].nunique(), "eRegulons_pp": pp["eRegulon_name"].nunique(),
               "eRegulons_pp_direct": pp.loc[~pp["is_extended"].astype(bool), "eRegulon_name"].nunique(), "TFs_pp": pp["TF"].nunique(),
               "TF_gene_links_pp": len(pp[["TF", "Gene"]].drop_duplicates()), "TF_region_links_pp": len(pp[["TF", "Region"]].drop_duplicates()),
               "region_gene_links_pp": len(links),
               "median_genes_per_eRegulon_pp": float(pp.groupby("eRegulon_name")["Gene"].nunique().median())}
        summary.append(row)
        for _, _, b in BANDS:
            dist_rows.append({"config": c, "band": b, "region_gene_links_pp": int((links["band"] == b).sum()),
                              "fraction": round(float((links["band"] == b).mean()), 4)})
        g = pp.groupby("TF").agg(genes=("Gene", lambda s: set(s)), regions=("Region", lambda s: set(s)))
        per_tf[c] = g
        auc[c] = gene_auc(c, set(g.index))
        print(f"{c}: {row['eRegulons_pp']} +/+ eRegulons, {row['TFs_pp']} TFs, {row['TF_gene_links_pp']} TF-gene links", flush=True)

    S = pd.DataFrame(summary)
    S.to_csv(HERE / "search_space_summary.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(HERE / "search_space_link_distance.csv", index=False)

    # per TF
    tfs = sorted(set().union(*[set(v.index) for v in per_tf.values()]))
    rows = []
    for tf in tfs:
        r = {"TF": tf}
        g2 = per_tf["config2"]["genes"].get(tf, set())
        for c in CONFIGS:
            g = per_tf[c]["genes"].get(tf, set()); rg = per_tf[c]["regions"].get(tf, set())
            r[f"genes_{c}"] = len(g); r[f"regions_{c}"] = len(rg)
            if c != "config2":
                r[f"jaccard_genes_vs_config2_{c}"] = round(len(g & g2) / len(g | g2), 3) if (g or g2) else np.nan
                if tf in auc[c] and tf in auc["config2"]:
                    a, b = auc[c][tf].align(auc["config2"][tf], join="inner")
                    r[f"auc_r_vs_config2_{c}"] = round(float(np.corrcoef(a, b)[0, 1]), 3)
            if tf in auc[c]:
                x = auc[c][tf]; lab = labels.reindex(x.index)
                xa, xb = x[lab == RA], x[lab == C1H1]
                sd = np.sqrt((xa.var() + xb.var()) / 2)
                r[f"ra_vs_c1h1_smd_{c}"] = round(float((xa.mean() - xb.mean()) / sd), 3) if sd > 0 else np.nan
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(HERE / "search_space_per_tf.csv", index=False)
    rt = pd.DataFrame([x for c in ("config3", "config4") for x in runtime(c)])
    rt.to_csv(HERE / "search_space_runtime.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n", S.set_index("config").T.to_string())
    print("\nregion-gene links used by +/+ eRegulons, by distance:\n", pd.DataFrame(dist_rows).pivot(index="band", columns="config", values="region_gene_links_pp").to_string())
    print("\nrun time (min):\n", rt.pivot(index="rule", columns="config", values="minutes").to_string() if len(rt) else "no logs")
    for c in ("config3", "config4"):
        both = T[(T["genes_config2"] > 0) & (T[f"genes_{c}"] > 0)]
        print(f"\n{c} vs config2: TFs in both {len(both)}, only config2 {int(((T['genes_config2'] > 0) & (T[f'genes_{c}'] == 0)).sum())}, "
              f"only {c} {int(((T['genes_config2'] == 0) & (T[f'genes_{c}'] > 0)).sum())}; median Jaccard of target genes {both[f'jaccard_genes_vs_config2_{c}'].median():.2f}; "
              f"median AUC r {both[f'auc_r_vs_config2_{c}'].median():.2f}; "
              f"RA vs C1H-1 effect size: r across TFs with config2 {both[[f'ra_vs_c1h1_smd_{c}', 'ra_vs_c1h1_smd_config2']].corr().iloc[0, 1]:.2f}")
    key = T[T["TF"].isin(["MAFB", "AR", "NR3C2", "MEF2C", "SOX5", "NFIA"])]
    print("\nselected TFs:\n", key.set_index("TF")[[c for c in T.columns if c.startswith(("genes_", "jaccard", "auc_r", "ra_vs"))]].to_string())


if __name__ == "__main__":
    main()
