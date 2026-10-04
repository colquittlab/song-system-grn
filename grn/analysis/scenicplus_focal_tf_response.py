#!/usr/bin/env python3
"""How does the presence of a TF's eRegulon respond to the motif-stage parameters across the correct-input configs?

Uses scenicplus_coherence_eregulons.csv (every eRegulon, all sign classes) and config_parameters.tsv. Configs: the
original sweep (1-11) and the looser-regime sweep (16-28 that have results). Configs 9 and 10 change non-motif
parameters (extended annotation, search space), and 12/13/14/15 and 29-33 are controls, so they are left out of the
parameter association (kept in the presence table). Per focal TF and sign class:
  * in how many configs the eRegulon exists
  * for each of the six swept parameters, AUC = P(a config WITH the eRegulon has a higher value than one WITHOUT)
    (0.5 = no relation; direction matters: lower NES / higher rank, p-value, fdr and lower log2FC are "looser")
  * a combined looseness score (mean z of -NES, rank, p-value, -log2FC, fdr; the AUC threshold is not directional and
    is left out of it) and its AUC
With ~20 configs these are descriptive: the design is a Latin hypercube, so parameters are roughly independent, but one
run per config and few configs mean only large effects mean anything.

    python scenicplus_focal_tf_response.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
E = pd.read_csv(HERE / "scenicplus_coherence_eregulons.csv")
E["sign"] = E["sign"].astype(str)
P = pd.read_csv(HERE.parent / "ra-arco-hvc-nc/scenicplus/config_parameters.tsv", sep="\t").set_index("config")
FOCAL = ["AR", "ALX4", "EMX2", "LHX2", "MAFB"]
PARAMS = ["ctx_nes_threshold", "ctx_auc_threshold", "ctx_rank_threshold", "dem_adj_pval_thr", "dem_log2fc_thr", "motif_similarity_fdr"]

have = set(E.config.unique())
corr_input = [c for c in P.index if int(c[6:]) in list(range(1, 12)) + list(range(16, 29)) and c in have]
assoc = [c for c in corr_input if c not in ("config9", "config10")]
print(f"configs with results in the correct-input design: {len(corr_input)}; used for the parameter association: {len(assoc)}")

X = P.loc[assoc, PARAMS].astype(float)
Z = (X - X.mean()) / X.std()
loose = (-Z["ctx_nes_threshold"] + Z["ctx_rank_threshold"] + Z["dem_adj_pval_thr"] - Z["dem_log2fc_thr"] + Z["motif_similarity_fdr"]) / 5


def auc(values, present):
    a, b = values[present], values[~present]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    return float(np.mean([(x > y) + 0.5 * (x == y) for x in a for y in b]))


rows, pres_rows = [], []
for tf in FOCAL:
    for sign in ("any", "+/+", "+/-"):
        t = E[(E.TF == tf) & E.config.isin(corr_input)]
        if sign != "any":
            t = t[t.sign == sign]
        present_all = set(t.config)
        pres_rows.append({"TF": tf, "sign": sign, "n_present": len(present_all), "n_configs": len(corr_input),
                          "configs": " ".join(sorted(present_all, key=lambda c: int(c[6:])))})
        present = pd.Series([c in present_all for c in assoc], index=assoc)
        if present.sum() in (0, len(assoc)):
            continue
        r = {"TF": tf, "sign": sign, "n_present": int(present.sum()), "n_configs": len(assoc), "AUC_looseness": round(auc(loose.values, present.values), 2)}
        for p in PARAMS:
            r[f"AUC_{p}"] = round(auc(X[p].values, present.values), 2)
        rows.append(r)

pres = pd.DataFrame(pres_rows)
resp = pd.DataFrame(rows)
pres.to_csv(HERE / "scenicplus_focal_tf_presence.csv", index=False)
resp.to_csv(HERE / "scenicplus_focal_tf_response.csv", index=False)
pd.set_option("display.width", 260, "display.max_columns", 20, "display.max_colwidth", 90)
print("\nPresence of an eRegulon by sign class (configs 1-11 + 16-28 that have results):")
print(pres[["TF", "sign", "n_present", "n_configs"]].to_string(index=False))
print("\nAUC of each parameter for 'eRegulon present' (>0.5: present configs have the higher value; <0.5: the lower):")
print(resp.drop(columns=["n_configs"]).to_string(index=False))
