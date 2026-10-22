"""figures/figure2/common.py -- inputs and shared helpers for the Figure 2 Python panel scripts (E, F, G).

Every Python panel script imports this module and nothing else in this folder, so each panel can be regenerated on its
own. Inputs are the label-transfer and joint-embedding tables written by the Seurat RPCA runs (README.md), tracked in
the repository under snrna/integration/gg_glutonly_hybrid/rpca/results/cca/; set FIG2_RPCA_RESULTS to read them from
elsewhere. Panels are written to FIG2_OUT (default: figures/figure2/output/, not tracked).

The pieces the original scripts imported from outside this repository are copied here, so nothing else is needed:
  - symmetrize() and best_per_cluster() from finch-integration-toolkit/reciprocal_score.py (a private lab repository);
  - the position palette and fixed chicken-column order from snrna/integration/position_palette.py.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RPCA_RESULTS = Path(os.environ.get("FIG2_RPCA_RESULTS", REPO_ROOT / "snrna/integration/gg_glutonly_hybrid/rpca/results/cca"))
FIG2_OUT = Path(os.environ.get("FIG2_OUT", REPO_ROOT / "figures/figure2/output"))
FIG2_OUT.mkdir(parents=True, exist_ok=True)

# Finch glutamatergic subset with Glut-SATB2-1 kept, chicken excitatory reference without the precursor types
# (Ex_Pre_KCNH7, Ex_Pre_SATB2); RPCA on SCTransform data, 40 dimensions, k.filter 200 (nominal; forced to NA under SCT).
RUN_TAG = "gg_glutonly_hybrid_plusSATB2_noPre"
KF, DIMS = 200, 40
KA_VALUES = [2, 5, 10, 20, 30, 50, 75, 100]


def matrix_paths(ka, tag=RUN_TAG):
    """Forward (finch cells labelled from the chicken reference) and reverse label-transfer matrices of one k.anchor:
    rows = query cluster, columns = reference cluster, values = fraction of the query cluster's cells assigned to the
    reference cluster x their mean prediction score (Seurat TransferData). 'mouse' in the file names is the chicken."""
    f = RPCA_RESULTS / f"{tag}_rpca_SCT_finch_from_mouse_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    r = RPCA_RESULTS / f"{tag}_rpca_SCT_mouse_from_finch_ka{ka}_kf{KF}_d{DIMS}_matrix.csv"
    return f, r


# ---- reciprocal transfer score (finch-integration-toolkit/reciprocal_score.py) -----------------------------------------
def symmetrize(fwd_csv, rev_csv):
    """fwd: index=finch cluster, columns=chicken cluster (frac*mean_score, finch<-chicken).
    rev: index=chicken cluster, columns=finch cluster (frac*mean_score, chicken<-finch).
    Returns S: index=finch cluster (fwd's own set only), columns=chicken cluster (union),
    S[i,j] = sqrt(fwd[i,j] * rev[j,i]), 0 where either side lacks that pair."""
    F = pd.read_csv(fwd_csv, index_col=0)
    R = pd.read_csv(rev_csv, index_col=0).T  # -> index=finch cluster, columns=chicken cluster
    all_chicken = F.columns.union(R.columns)
    F2 = F.reindex(columns=all_chicken, fill_value=0.0)
    R2 = R.reindex(index=F.index, columns=all_chicken, fill_value=0.0)
    S = np.sqrt(F2.clip(lower=0).values * R2.clip(lower=0).values)
    return pd.DataFrame(S, index=F.index, columns=all_chicken)


def best_per_cluster(fwd_csv, rev_csv):
    S = symmetrize(fwd_csv, rev_csv)
    return S.max(axis=1), S.idxmax(axis=1), S


# ---- finch cell types and positions (snrna/integration/position_palette.py) --------------------------------------------
POSITION_COLORS = {"hvc": "#ff7f00", "lman": "#1f78b4", "ra": "#e31a1c", "nc": "#fdbf6f", "nr": "#a6cee3", "arco": "#fb9a99"}
POSITION_LABELS = {"hvc": "HVC", "lman": "LMAN", "ra": "RA", "nc": "NC", "nr": "NR", "arco": "Arco"}
POSITION_ORDER = ["hvc", "lman", "ra", "nc", "nr", "arco"]
CLUSTER_POSITION = {
    "Glut-DACH2-HVCra": "hvc", "Glut-DACH2-HVCra-Int": "hvc", "Glut-DACH2-HVCx": "hvc",
    "Glut-DACH2-LMANco": "lman", "Glut-DACH2-LMANsh": "lman", "Glut-CACNA1H-RA": "ra",
    "Glut-DACH2-1": "nc", "Glut-DACH2-2": "nc", "Glut-DACH2-3": "nc", "Glut-DACH2-4": "nc",
    "Glut-DACH2-5": "nr", "Glut-DACH2-6": "nr", "Glut-DACH2-7": "nr", "Glut-DACH2-8": "nr",
    "Glut-CACNA1H-1": "arco", "Glut-CACNA1H-2": "arco", "Glut-CACNA1H-3": "arco", "Glut-CACNA1H-4": "arco",
}
SONG_CLUSTERS = {c for c, p in CLUSTER_POSITION.items() if p in ("hvc", "lman", "ra")}
# Glut-DACH2-HVCra-Int (transcriptionally intermediate between HVCra and the surrounding nidopallium) is left out of every
# finch-chicken comparison
EXCLUDED_CLUSTERS = {"Glut-DACH2-HVCra-Int"}
SATB2_CLUSTER = "Glut-SATB2-1"
LILAC = "#cab2d6"   # Glut-SATB2-1's row colour (common_aesthetics.R's 'fieldl' slot)

# Finch row order of the heatmaps: song types, then nidopallial DACH2, Glut-SATB2-1, arcopallial CACNA1H
ROW_ORDER = [
    "Glut-DACH2-HVCra", "Glut-DACH2-HVCx", "Glut-DACH2-LMANco", "Glut-DACH2-LMANsh", "Glut-CACNA1H-RA",
    "Glut-DACH2-1", "Glut-DACH2-2", "Glut-DACH2-3", "Glut-DACH2-4",
    "Glut-DACH2-5", "Glut-DACH2-6", "Glut-DACH2-7", "Glut-DACH2-8",
    SATB2_CLUSTER,
    "Glut-CACNA1H-1", "Glut-CACNA1H-2", "Glut-CACNA1H-3", "Glut-CACNA1H-4",
]


def song_group(c):
    return "song" if c in SONG_CLUSTERS else "non-song"


def row_color(c):
    return LILAC if c == SATB2_CLUSTER else POSITION_COLORS[CLUSTER_POSITION[c]]


# Fixed chicken-column order: grouped by the marker after "Ex_" in this order, alphabetical within a group; unlisted
# prefixes (e.g. Ex_Pre_*) after, also alphabetical.
CHICKEN_PREFIX_ORDER = ["DACH2", "TSHZ2", "CACNA1H", "SATB2", "KIAA1217", "BCL6"]


def chicken_col_order(cols):
    def key(c):
        prefix = str(c).removeprefix("Ex_").split("_")[0]
        rank = CHICKEN_PREFIX_ORDER.index(prefix) if prefix in CHICKEN_PREFIX_ORDER else len(CHICKEN_PREFIX_ORDER)
        return (rank, str(c))
    return sorted(cols, key=key)
