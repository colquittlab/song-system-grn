#!/usr/bin/env python3
"""Summarize the SCENIC+ parameter sweep: one row per config, from the tables all.qmd exports.

Reads <out_dir>/configN/scenicplus_eRegulons.txt for every config found, joins the swept parameters from
scenicplus/config_parameters.tsv, and writes scenicplus_sweep_summary.csv next to this script (small, so
tracked) and into out_dir. Overlap columns compare each config with REF (config37).

    python scenicplus_sweep_summary.py [out_dir]
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else
           "~/ssd/rstudio/multiome/motor-pathway/scenicplus/motor-pathway_scenicplus_v2_hybrid_all").expanduser()
PARAMS = HERE.parent / "ra-arco-hvc-nc" / "scenicplus" / "config_parameters.tsv"
REF = "config37"   # main comparison config: 40 topics, new cisTarget database on the current consensus, ArchR song-pair DARs


def load(cfg):
    d = pd.read_csv(OUT / cfg / "scenicplus_eRegulons.txt", sep="\t", usecols=["TF", "Gene", "eRegulon_name"])
    d["sign"] = d.eRegulon_name.str.extract(r"_(\+/\+|-/-|\+/-|-/\+)")[0]
    d["type"] = d.eRegulon_name.str.extract(r"_(direct|extended)_")[0]
    return d


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else float("nan")


def tf_status(d, tf):
    """e.g. 'direct +/+ (38), extended -/- (10)' or '-' -- targets per eRegulon of this TF."""
    s = d[d.TF == tf]
    if s.empty:
        return "-"
    return ", ".join(f"{t} {sg} ({g.Gene.nunique()})" for (t, sg), g in s.groupby(["type", "sign"]))


cfgs = sorted((p.name for p in OUT.glob("config[0-9]*") if (p / "scenicplus_eRegulons.txt").exists()),
              key=lambda c: int(c[6:]))
data = {c: load(c) for c in cfgs}
pos = {c: d[d.sign == "+/+"] for c, d in data.items()}
tfs = {c: set(p.TF) for c, p in pos.items()}
edges = {c: set(map(tuple, p[["TF", "Gene"]].drop_duplicates().values)) for c, p in pos.items()}

rows = []
for c in cfgs:
    names = data[c].drop_duplicates("eRegulon_name")
    tg = pos[c].groupby("eRegulon_name").Gene.nunique()
    rows.append({
        "config": c,
        "eRegulons": len(names),
        "pos_eRegulons": int((names.sign == "+/+").sum()),
        "neg_eRegulons": int((names.sign == "-/-").sum()),
        "pos_TFs": len(tfs[c]),
        "median_targets_per_pos_eRegulon": float(tg.median()),
        "pos_TF_gene_edges": len(edges[c]),
        f"TF_jaccard_vs_{REF}": round(jaccard(tfs[c], tfs[REF]), 3) if REF in tfs else float("nan"),
        f"edge_jaccard_vs_{REF}": round(jaccard(edges[c], edges[REF]), 3) if REF in edges else float("nan"),
        "MAFB": tf_status(data[c], "MAFB"),
        "EMX2": tf_status(data[c], "EMX2"),
    })
out = pd.DataFrame(rows)
if PARAMS.exists():
    p = pd.read_csv(PARAMS, sep="\t")[["config", "changed_from_config1"]]
    out = out.merge(p, on="config", how="left")

pd.set_option("display.width", 250, "display.max_columns", 30)
print(out.to_string(index=False))
out.to_csv(HERE / "scenicplus_sweep_summary.csv", index=False)
out.to_csv(OUT / "scenicplus_sweep_summary.csv", index=False)
