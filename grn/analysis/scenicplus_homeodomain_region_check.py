#!/usr/bin/env python3
"""Do homeodomain eRegulons' regions pull toward interneuron (LGE) chromatin, and does removing GABA region sets change it?

For each config (default: every config with results), over ALL sign classes and both cistrome types:
  * share of homeodomain-TF eRegulons whose REGION activity peaks in GABA-LGE-1/-2, and in any GABA cluster
  * mean z-score (across clusters, within each regulon) of region activity in GABA-LGE-1 for homeodomain vs other TFs
  * eRegulons whose regions peak in LGE but whose genes peak elsewhere, and the homeodomain share among them
  * the focal TFs (ALX4, EMX2, AR, LHX2): every eRegulon with its sign, size, gene-activity peak, region-activity peak
Homeodomain = DBD class "Homeodomain" in the Lambert et al. TF table (data/lambert_tfs/tfs.csv).

Writes scenicplus_homeodomain_check.csv (one row per config) and scenicplus_homeodomain_focal_tfs.csv.

    python scenicplus_homeodomain_region_check.py [config ...]
"""
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RES = Path("/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/"
           "motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/ra-arco-hvc-nc_hybrid/results")
FOCAL = ("ALX4", "EMX2", "AR", "LHX2")
LGE = ("GABA-LGE-1", "GABA-LGE-2")
s = lambda x: x.decode() if isinstance(x, bytes) else x

dbd = dict(zip(*pd.read_csv(REPO / "data/lambert_tfs/tfs.csv")[["Name", "DBD"]].values.T))
is_hd = lambda tf: "homeodomain" in str(dbd.get(tf, "")).lower()


def profiles(cfg):
    """{'gene'|'region': DataFrame regulon x cluster of mean AUC} for all regulons, direct + extended."""
    outs = RES / cfg / "outs"
    full = outs / "scplusmdata.h5mu"
    f = h5py.File(full if full.exists() else outs / "scplusmdata_slim.h5mu", "r")
    o = f["mod"]["scRNA_counts"]["obs"]
    c2c = dict(zip((s(x) for x in o["_index"][:]), np.array([s(c) for c in o["cluster"]["categories"][:]])[o["cluster"]["codes"][:]]))
    out = {}
    for basis in ("gene", "region"):
        frames = []
        for kind in ("direct", "extended"):
            m = f["mod"][f"{kind}_{basis}_based_AUC"]
            names = [s(x) for x in m["var"]["_index"][:]]
            cl = np.array([c2c[s(x)] for x in m["obs"]["Cell"][:]])
            X = np.asarray(m["X"][:], dtype=np.float32)
            clusters = sorted(set(cl))
            frames.append(pd.DataFrame(np.vstack([X[cl == c].mean(0) for c in clusters]).T, index=names, columns=clusters))
        out[basis] = pd.concat(frames)
    f.close()
    return out


def rho(a, b):
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


cfgs = sys.argv[1:] or sorted((p.name for p in RES.glob("config*") if (p / "outs" / "scplusmdata.h5mu").exists()
                               or (p / "outs" / "scplusmdata_slim.h5mu").exists()), key=lambda c: int(c[6:]))
rows, focal_rows = [], []
for cfg in cfgs:
    P = profiles(cfg)
    G, R = P["gene"], P["region"]
    sizes = {n.rsplit("_(", 1)[0]: n for n in G.index}          # keep the gene-signature name for the target-gene count
    G.index = [n.rsplit("_(", 1)[0] for n in G.index]            # gene and region names differ only in their _(Ng)/_(Nr) suffix
    R.index = [n.rsplit("_(", 1)[0] for n in R.index]
    assert G.index.is_unique and set(G.index) == set(R.index), "gene/region regulon sets differ"
    R = R.reindex(G.index)
    tf = np.array([n.split("_")[0] for n in G.index])
    hd = np.array([is_hd(t) for t in tf])
    gpk, rpk = G.idxmax(1).values, R.idxmax(1).values
    Z = R.sub(R.mean(1), axis=0).div(R.std(1) + 1e-9, axis=0)
    in_lge = np.isin(rpk, LGE)
    in_gaba = np.array([str(c).startswith("GABA") for c in rpk])
    disc = in_lge & ~np.array([str(c).startswith("GABA") for c in gpk])
    rows.append({
        "config": cfg, "n_eregulons": len(G), "n_homeodomain": int(hd.sum()),
        "hd_region_peak_in_LGE": round(float(in_lge[hd].mean()), 3), "other_region_peak_in_LGE": round(float(in_lge[~hd].mean()), 3),
        "hd_region_peak_in_any_GABA": round(float(in_gaba[hd].mean()), 3),
        "hd_mean_z_LGE1": round(float(Z.loc[hd, "GABA-LGE-1"].mean()), 2), "other_mean_z_LGE1": round(float(Z.loc[~hd, "GABA-LGE-1"].mean()), 2),
        "n_region_LGE_gene_elsewhere": int(disc.sum()), "homeodomain_share_of_those": round(float(hd[disc].mean()), 3) if disc.any() else np.nan})
    for name in G.index:
        t = name.split("_")[0]
        if t in FOCAL:
            sign = re.search(r"_(\+/\+|-/-|\+/-|-/\+)", name).group(1)
            focal_rows.append({"config": cfg, "TF": t, "eRegulon": name, "sign": sign,
                               "n_genes": int(re.search(r"[(](\d+)g", sizes[name]).group(1)),
                               "gene_peak": G.loc[name].idxmax(), "region_peak": R.loc[name].idxmax(),
                               "rho_gene_vs_region": round(rho(G.loc[name].values, R.loc[name].values), 3)})
    print(cfg, "done", flush=True)

def merge_write(new, path):
    """Replace the rows of the configs just computed, keep the others, so runs accumulate."""
    if path.exists():
        old = pd.read_csv(path)
        new = pd.concat([old[~old.config.isin(new.config.unique())], new])
    new = new.sort_values("config", key=lambda c: c.str[6:].astype(int), kind="stable")
    new.to_csv(path, index=False)
    return new


summ = merge_write(pd.DataFrame(rows), HERE / "scenicplus_homeodomain_check.csv")
merge_write(pd.DataFrame(focal_rows), HERE / "scenicplus_homeodomain_focal_tfs.csv")
pd.set_option("display.width", 250, "display.max_columns", 20)
print("\n" + summ.to_string(index=False))
