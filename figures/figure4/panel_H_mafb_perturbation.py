#!/usr/bin/env python3
"""Figure 4, panel H -- simulated MAFB knockout (top) and overexpression (bottom) in glutamatergic neurons.

SCENIC+ in-silico perturbation on the glutamatergic-only SCENIC+ run (ra-arco-hvc-nc_glut_hybrid, config3): gradient-boosting regressors predict each
target gene from its eRegulon TFs (the top-5,000 eRegulon target genes, plus every eRegulon target of MAFB); MAFB is set to 0 (knockout) or 10 counts
(overexpression) and the effect is propagated for 5 iterations (SCENIC+ simulate_perturbation). The shift of each cell is projected on the Harmony-corrected
PCA of the eRegulon AUC and drawn as a quiver field over the cells, colored by cell type (arrows >= 2% of the longest; neighbour sampling seeded with 0).

This is the notebook grn/ra-arco-hvc-nc_glut_hybrid/perturbation/perturbation_parallel_template.ipynb (branch grn-glut-hybrid, commit 46dd519; the
version that wrote the figure files), extracted to a script: the notebook's code cells follow unchanged, preceded by the settings of the run that made the
figure (run_quiver4.log), which can be overridden with the same environment variables (see the first cell). The per-TF simulation runs in the child script
panel_H/perturbation_parallel.py. Regressors, simulations and the arrow fields are cached in the output folder and reused when present (training the
regressors takes hours; with the caches the run takes ~25 min, most of it the Harmony PCA and the gene-line plots).

Input: <SCENIC+ glut run>/results/config3/outs/scplusmdata_perturb.h5mu, made by panel_H/assemble_perturbation_h5mu.py (raw counts of the glutamatergic
AnnData + config3 AUC matrices and eRegulon tables).
Output: the notebook's figures in <results>/config3/perturbation_mafb-targets/; the two panel files are copied to FIG4_OUT as
panel_H_pca_arrows_kd_MAFB.pdf (knockout) and panel_H_pca_arrows_oe_MAFB.pdf (overexpression).

Environment: the scenicplus conda environment (Python 3.11, scenicplus, scanpy, mudata, harmonypy).
    python figures/figure4/panel_H_mafb_perturbation.py
"""
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
# settings of the run that made the figure (each can be overridden from the environment)
SCENICPLUS_JUPYTER_DIR = "/hdd/jupyter/brad/scenicplus/motor-pathway_multiome/motor-pathway_multiome_seurat_cellbender.0.05_preprocess_cr/"   # config/paths.R
for k, v in {
    "PERTURB_INPUT_DIR": os.path.join(SCENICPLUS_JUPYTER_DIR, "ra-arco-hvc-nc_glut_hybrid", "results"),
    "PERTURB_CONFIGS": "config3",
    "PERTURB_H5MU": "scplusmdata_perturb.h5mu",
    "PERTURB_KO_TFS": "MAFB",
    "PERTURB_OE_TFS": "auto",
    "PERTURB_FORCE_TARGETS_OF": "MAFB",
    "PERTURB_OUT_NAME": "perturbation_mafb-targets",
    "PERTURB_ARROW_STYLE": "quiver",
    "PERTURB_SEED": "0",
    "PERTURB_REDO_PERTURBATIONS": "0",
    "PERTURB_N_JOBS": "40",
}.items():
    os.environ.setdefault(k, v)
FIG4_OUT = os.environ.get("FIG4_OUT", os.path.join(HERE, "output"))
os.chdir(os.path.join(HERE, "panel_H"))   # the notebook finds perturbation_parallel.py in its working directory

# ======================================================================================================== notebook code cell 0
# ── User-configurable parameters ──────────────────────────────────────────────

# Directory containing config* subdirectories
import os
INPUT_DIR = os.environ.get("PERTURB_INPUT_DIR", ".")
# Only these configs (comma-separated in PERTURB_CONFIGS, e.g. "config3"); None = every config* found
CONFIGS = os.environ["PERTURB_CONFIGS"].split(",") if os.environ.get("PERTURB_CONFIGS") else None
# MuData with scRNA_counts + the AUC matrices + eRegulon tables. scplusmdata.h5mu is the full SCENIC+ file; configs run with
# build_scplus_mudata: false have none, so assemble_perturbation_h5mu.py writes scplusmdata_perturb.h5mu (real raw counts, placeholder ATAC)
H5MU_NAME = os.environ.get("PERTURB_H5MU", "scplusmdata.h5mu")

# TFs to perturb for knockout (set_value=0).  None = auto-detect from eRegulon metadata.
KO_TFS = os.environ["PERTURB_KO_TFS"].split(',') if os.environ.get("PERTURB_KO_TFS") else ['EMX2', 'MAFB', 'BCL11B', 'NR2F2','RBFOX2', 'MEIS2', 'AR', 'PRDM5' ]  # env: comma-separated; e.g. ['EMX2', 'MAFB'] or None

# TFs to perturb for overexpression (set_value=10).  None = skip OE.  'auto' = same as KO.
OE_TFS = {'none': None, 'auto': 'auto'}.get(os.environ.get('PERTURB_OE_TFS', 'auto').lower(), os.environ.get('PERTURB_OE_TFS', 'auto').split(','))  # env: none | auto | comma-separated
OE_SET_VALUE = float(os.environ.get("PERTURB_OE_SET_VALUE", 10))  # on the scale of the expression matrix: 10 counts is NOT 10 on a log1p(CP10k) scale

# "counts" (SCENIC+ default: raw counts in scRNA_counts) or "lognorm" (log1p CP10k, from assemble_perturbation_h5mu.py --lognorm). Only changes how the
# plotted fold-change is computed (means of expm1(values), not means of logs); simulation and regressors use the matrix as given.
EXPRESSION_SCALE = os.environ.get("PERTURB_EXPRESSION_SCALE", "counts")

# Number of top genes to use for regressor training
NGENES = 5000

# TFs whose eRegulon target genes ALWAYS get a regressor, on top of the NGENES top-ranked genes (env: comma-separated). Without this a TF with a
# small eRegulon next to large ones may have most of its targets outside the top NGENES, so its perturbation reaches only a fraction of them.
FORCE_TARGETS_OF = os.environ["PERTURB_FORCE_TARGETS_OF"].split(",") if os.environ.get("PERTURB_FORCE_TARGETS_OF") else []

# Output folder name under each config directory (use a different name to keep a variant's results apart)
OUT_NAME = os.environ.get("PERTURB_OUT_NAME", "perturbation")

# How the PCA flow figure is drawn. "stream" = SCENIC+'s own streamplot (unchanged default; it drops arrows shorter than 15% of the longest, so a sparse
# field, e.g. a TF expressed in few cells, can come out empty). "quiver" = the same field drawn as individual arrows at the grid points, written to
# pca_arrows_<prefix>_<TF>.pdf. "both" = both figures from one field computation.
ARROW_STYLE = os.environ.get("PERTURB_ARROW_STYLE", "stream")
# quiver only: keep arrows at least this fraction of the longest one (SCENIC+'s streamplot rule is 0.15)
ARROW_MIN_FRAC = float(os.environ.get("PERTURB_ARROW_MIN_FRAC", 0.02))
# Seed for the random neighbour sampling inside the flow-field projection (np.random.choice, unseeded in SCENIC+): None = unseeded as before
SEED = int(os.environ["PERTURB_SEED"]) if os.environ.get("PERTURB_SEED") else None
# PCA flow figures: side length in inches (small-multiple size for a manuscript panel), cell scatter size (the original used 0.5 on a 6 in figure)
FIG_SIZE = float(os.environ.get("PERTURB_FIG_SIZE", 1.8))
POINT_SIZE = float(os.environ.get("PERTURB_POINT_SIZE", 1.2))
POINT_ALPHA = float(os.environ.get("PERTURB_POINT_ALPHA", 0.7))  # semi-transparent points keep the arrows visible

# Whether to re-train regressors even if regressors.pkl exists
REDO_REGRESSORS = os.environ.get("PERTURB_REDO_REGRESSORS", "0") == "1"

# Whether to re-run perturbation simulations even if parquet files exist
REDO_PERTURBATIONS = os.environ.get("PERTURB_REDO_PERTURBATIONS", "1") == "1"  # default True as in the original; set 0 to reuse cached parquets (e.g. to only redraw plots)

# Number of perturbation iterations
N_ITER = 5

# Number of top target genes to show in line plots
N_GENES_TO_SHOW = 9

# Dimensionality reduction key to use for PCA plots (Harmony-corrected)
REDUCTION_KEY = "X_pca_harmony"

# obs column name for Harmony batch correction
HARMONY_BATCH_KEY = "replicate"

# Number of jobs for regressor training
N_JOBS_REGRESSORS = int(os.environ.get("PERTURB_N_JOBS", 20))

# obs column name for cluster labels
CLUSTER_COL = "cluster"

# ======================================================================================================== notebook code cell 1
import mudata
import os
import sys
import glob as globmod
import scanpy as sc
import anndata
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for reliable file saving
import matplotlib.pyplot as plt
import adjustText
import pickle
import joblib
import numpy as np
import pandas as pd
import h5py

from anndata._io.specs.registry import read_elem
from sklearn.ensemble import GradientBoostingRegressor
from joblib import Parallel, delayed
from tqdm import tqdm
from scenicplus.simulation import (
    simulate_perturbation,
    plot_perturbation_effect_in_embedding,
    _project_perturbation_in_embedding,  # private SCENIC+ helpers, reused so the quiver draws exactly the field the streamplot would
    _calculate_grid_arrows,
)
import matplotlib.colors
import matplotlib.patheffects as path_effects
matplotlib.rcParams['pdf.fonttype'] = 42  # embed Arial as TrueType (editable text), not Type 3

# ======================================================================================================== notebook code cell 2
def mudata_to_pandas(h5mu_fname, group_path):
    """Read an HDF5 group from a .h5mu file into a pandas DataFrame."""
    with h5py.File(h5mu_fname, 'r') as f:
        group = f[group_path]
        data = {}
        for key, item in group.items():
            if isinstance(item, h5py.Dataset):
                data[key] = read_elem(item)
    return pd.DataFrame(data)


def load_config_data(h5mu_fname, cluster_col):
    """Load all data needed from a single config's h5mu file.
    
    Returns
    -------
    dict with keys: scRNA_counts, eRegulon_gene_AUC, direct_e_regulon_metadata,
                    extended_e_regulon_metadata, gene_to_TF, gene_tf_df, celltypes
    """
    direct_gene_based_AUC = mudata.read(f"{h5mu_fname}/direct_gene_based_AUC")
    extended_gene_based_AUC = mudata.read(f"{h5mu_fname}/extended_gene_based_AUC")
    scRNA_counts = mudata.read(f"{h5mu_fname}/scRNA_counts")

    direct_e_regulon_metadata = mudata_to_pandas(h5mu_fname, "uns/direct_e_regulon_metadata")
    extended_e_regulon_metadata = mudata_to_pandas(h5mu_fname, "uns/extended_e_regulon_metadata")

    # Build eRegulon AUC
    eRegulon_gene_AUC = anndata.concat(
        [direct_gene_based_AUC, extended_gene_based_AUC], axis=1
    )
    eRegulon_gene_AUC.obs = scRNA_counts.obs

    # Gene-TF mapping
    gene_tf_df = pd.concat([
        direct_e_regulon_metadata[["Gene", "TF", "triplet_rank"]].drop_duplicates(),
        extended_e_regulon_metadata[["Gene", "TF", "triplet_rank"]].drop_duplicates(),
    ]).drop_duplicates()
    gene_to_TF = gene_tf_df.groupby("Gene")["TF"].apply(list).to_dict()

    celltypes = set(scRNA_counts.obs[cluster_col])

    # Combined metadata for plotting
    combined_metadata = pd.concat([direct_e_regulon_metadata, extended_e_regulon_metadata])
    combined_metadata['importance_x_rho_TF2G'] = (
        combined_metadata['importance_TF2G'] * combined_metadata['rho_TF2G']
    )
    combined_metadata['importance_x_abs_rho_TF2G'] = (
        combined_metadata['importance_TF2G'] * np.abs(combined_metadata['rho_TF2G'])
    )

    return {
        'scRNA_counts': scRNA_counts,
        'eRegulon_gene_AUC': eRegulon_gene_AUC,
        'direct_e_regulon_metadata': direct_e_regulon_metadata,
        'extended_e_regulon_metadata': extended_e_regulon_metadata,
        'gene_to_TF': gene_to_TF,
        'gene_tf_df': gene_tf_df,
        'combined_metadata': combined_metadata,
        'celltypes': celltypes,
    }


def make_color_dict(celltypes):
    """Create a color dictionary for the given cell types (sorted order).

    Saturated, well-separated colours: matplotlib's tab10 without its gray (a gray cluster reads as "missing") = 9 colours, then tab20
    if there are more cell types. The earlier tab20-resampled-to-n palette gave pastel neighbours that were hard to tell apart.
    """
    base = [c for i, c in enumerate(plt.get_cmap('tab10').colors) if i != 7]
    palette = base if len(celltypes) <= len(base) else list(plt.get_cmap('tab20').colors)
    return {ct: palette[i % len(palette)] for i, ct in enumerate(sorted(celltypes))}


def train_regressors(scRNA_counts, genes_to_use, gene_to_TF, out_dir, redo=False, n_jobs=32):
    """Train or load GBR regressors for gene expression prediction."""
    regressor_fname = os.path.join(out_dir, 'regressors.pkl')

    existing = {}
    if not redo and os.path.exists(regressor_fname):
        print(f"  Loading existing regressors from {regressor_fname}")
        with open(regressor_fname, 'rb') as f:
            existing = pickle.load(f)
        # top up: train only the requested genes that are not in the file yet (e.g. FORCE_TARGETS_OF added after a first run)
        missing = [g for g in genes_to_use if g in gene_to_TF and g not in existing]
        if not missing:
            return existing
        print(f"  {len(missing)} requested genes have no regressor yet; training those only")
        genes_to_use = missing

    print(f"  Training regressors ({len(genes_to_use)} genes, n_jobs={n_jobs}) ...")
    df_EXP = scRNA_counts.to_df()
    genes = [g for g in genes_to_use if g in gene_to_TF]

    def _train_one_gene(gene):
        predictor_TF = gene_to_TF[gene].copy()
        if gene in predictor_TF:
            predictor_TF.remove(gene)
        if len(predictor_TF) == 0:
            return None
        X = df_EXP[predictor_TF].to_numpy()
        if len(X.shape) == 1:
            X = X.reshape(-1, 1)
        y = df_EXP[gene].to_numpy()
        reg = GradientBoostingRegressor(
            learning_rate=0.01, n_estimators=500, max_features=0.1
        )
        reg.fit(X, y)
        return (gene, predictor_TF, reg)

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_train_one_gene)(gene) for gene in tqdm(genes)
    )

    regressors = dict(existing)
    regressors.update({
        gene: (tfs, reg)
        for r in results if r is not None
        for gene, tfs, reg in [r]
    })

    with open(regressor_fname, 'wb') as f:
        pickle.dump(regressors, f)
    print(f"  Saved {len(regressors)} regressors to {regressor_fname}")
    return regressors


def _sv_tag(set_value):
    """Return a filename-safe tag encoding the set_value, e.g. 'sv0' or 'sv10'."""
    return f"sv{set_value:g}"


PERTURBATION_SCRIPT = os.environ.get("PERTURBATION_SCRIPT", os.path.join(os.getcwd(), "perturbation_parallel.py"))  # tracked next to this notebook


def load_or_run_perturbation(out_dir, h5mu_fname, regressor_fname, tfs, set_value=0, n_iter=5, redo=False):
    """Run perturbation_parallel.py for the given TFs, or load cached results.
    
    Parquet files are named perturb_{tf}_{sv_tag}_{iter}.parquet so that
    knockout (sv0) and overexpression (sv10) results never collide.
    Set redo=True to delete cached parquets and re-run all TFs.
    """
    tag = _sv_tag(set_value)

    # Check which TFs already have parquet results
    tfs_to_run = []
    cached_results = {}
    for tf in tfs:
        parquets = globmod.glob(os.path.join(out_dir, f"perturb_{tf}_{tag}_*.parquet"))
        if parquets and not redo:
            result = {}
            for i in range(n_iter + 1):
                pq = os.path.join(out_dir, f"perturb_{tf}_{tag}_{i}.parquet")
                if os.path.exists(pq):
                    result[i] = pd.read_parquet(pq, engine="fastparquet")
            cached_results[tf] = result
            print(f"    [cached] {tf} {tag} ({len(result)} iterations)")
        else:
            if redo and parquets:
                for pq in parquets:
                    os.remove(pq)
                print(f"    [redo] {tf} {tag} — deleted {len(parquets)} cached parquets")
            tfs_to_run.append(tf)

    if tfs_to_run:
        tfs_str = ' '.join(tfs_to_run)
        cmd = (
            f"{sys.executable} {PERTURBATION_SCRIPT}"  # same interpreter as this notebook (a bare `python` may not be the scenicplus env)
            f" --tfs {tfs_str}"
            f" --h5mu {h5mu_fname}"
            f" --regressor-pkl {regressor_fname}"
            f" --out-dir {out_dir}"
            f" --set-value {set_value}"
            f" --n-iter {n_iter}"
        )
        print(f"    Running: {cmd}")
        rc = os.system(cmd)
        if rc != 0:  # the original ignored the exit status and carried on to plot nothing
            raise RuntimeError(f"perturbation_parallel.py failed (exit status {rc}); see the output above")

        # Load the newly generated results
        for tf in tfs_to_run:
            result = {}
            for i in range(n_iter + 1):
                pq = os.path.join(out_dir, f"perturb_{tf}_{tag}_{i}.parquet")
                if os.path.exists(pq):
                    result[i] = pd.read_parquet(pq, engine="fastparquet")
            if result:
                cached_results[tf] = result

    return cached_results


def get_top_tfs(gene_tf_df, n=10):
    """Get the top N most-connected TFs from the eRegulon metadata."""
    tf_counts = gene_tf_df['TF'].value_counts()
    return tf_counts.head(n).index.tolist()


def get_target_genes(combined_metadata, tf, regressor_genes, n=9):
    """Get top target genes for a TF sorted by rho_TF2G, restricted to regressor_genes."""
    ind = combined_metadata.TF.str.contains(tf, regex=False)
    df_tf = combined_metadata[ind].drop_duplicates(subset='Gene')
    df_filt = df_tf.loc[~df_tf.Gene.str.contains('LOC')]
    df_filt = df_filt.loc[~df_filt.Gene.str.contains(tf, regex=False)]
    df_filt = df_filt.loc[df_filt.Gene.isin(regressor_genes)]
    df_filt = df_filt.sort_values('rho_TF2G', ascending=False)
    return df_filt.Gene.tolist()[1:n+1]


def get_target_genes_by_rank(combined_metadata, tf, regressor_genes, n=9):
    """Get top target genes for a TF sorted by triplet_rank, restricted to regressor_genes."""
    ind = combined_metadata.TF.str.contains(tf, regex=False)
    df_tf = combined_metadata[ind]
    df_filt = df_tf.loc[~df_tf.Gene.str.contains('LOC')]
    df_filt = df_filt.loc[~df_filt.Gene.str.contains(tf, regex=False)]
    df_filt = df_filt.loc[df_filt.Gene.isin(regressor_genes)]
    df_filt = df_filt.sort_values('triplet_rank', ascending=True)
    df_filt = df_filt.drop_duplicates(subset='Gene')
    return df_filt.Gene.tolist()[1:n+1]

# ======================================================================================================== notebook code cell 3
from matplotlib.patches import FancyArrowPatch


def _compute_log2fc(perturbation_over_iters, tf, gene, cell_line, cluster_col, eRegulon_gene_AUC, n_iter):
    """Return list of log2FC values across iterations for a gene/cell_line. Returns None if flat or missing."""
    _m = (lambda x: np.expm1(x)) if EXPRESSION_SCALE == "lognorm" else (lambda x: x)
    grouped = [
        _m(perturbation_over_iters[tf][j])
        .groupby(eRegulon_gene_AUC.obs[cluster_col]).mean()
        for j in range(n_iter + 1)
    ]
    if gene not in grouped[0].columns or cell_line not in grouped[0].index:
        return None
    baseline = grouped[0].loc[cell_line, gene]
    if baseline == 0:
        return None
    vals = [np.log2(grouped[j].loc[cell_line, gene] / baseline) for j in range(1, n_iter + 1)]
    if max(abs(v) for v in vals) < 1e-6:
        return None  # flat line — gene not in regressors
    return vals


def plot_gene_lines(perturbation_over_iters, tf, genes_to_show, eRegulon_gene_AUC,
                    cell_lines, cluster_col, out_dir, prefix="kd", n_iter=5):
    """Plot log2FC of target genes across perturbation iterations (all cell types side-by-side)."""
    num_plots = len(cell_lines)
    fig, ax = plt.subplots(1, num_plots, figsize=(5 * num_plots, 5))
    if num_plots == 1:
        ax = [ax]

    for i, cell_line in enumerate(cell_lines):
        for gene in genes_to_show:
            vals = _compute_log2fc(perturbation_over_iters, tf, gene, cell_line,
                                   cluster_col, eRegulon_gene_AUC, n_iter)
            if vals is None:
                continue
            ax[i].plot(np.arange(n_iter) + 1, vals, label=gene)

        ax[i].set_xlabel("Iteration")
        ax[i].set_title(cell_line)
        ax[i].axhline(y=0, color="black")
        ax[i].grid(color="gray")
        ax[i].set_axisbelow(True)
        ax[i].set_xticks(range(1, n_iter + 1))

    # Sync y-axis limits
    all_ylims = [a.get_ylim() for a in ax]
    ymin = min(y[0] for y in all_ylims)
    ymax = max(y[1] for y in all_ylims)
    plt.setp(ax, ylim=(ymin, ymax))

    ax[0].set_ylabel("Predicted $log{_2}FC$")
    for a in ax:
        a.legend(fontsize='small')

    fig.suptitle(f"{tf} ({prefix.upper()})")
    fig.tight_layout()
    save_path = os.path.join(out_dir, f"genes_plot_{prefix}_{tf}.pdf")
    fig.savefig(save_path, transparent=True)
    plt.close(fig)
    print(f"      Saved {save_path}")


def plot_gene_lines_per_celltype(perturbation_over_iters, tf, genes_to_show,
                                  eRegulon_gene_AUC, cell_lines, cluster_col,
                                  out_dir, prefix="kd", n_iter=5):
    """Plot log2FC of target genes for each cell type in a separate figure."""
    for cell_line in cell_lines:
        fig, ax = plt.subplots(1, 1)
        for gene in genes_to_show:
            vals = _compute_log2fc(perturbation_over_iters, tf, gene, cell_line,
                                   cluster_col, eRegulon_gene_AUC, n_iter)
            if vals is None:
                continue
            ax.plot(np.arange(n_iter) + 1, vals, label=gene)

        ax.set_xlabel("Iteration")
        ax.set_title(f"{tf} — {cell_line}")
        ax.axhline(y=0, color="black")
        ax.grid(color="gray")
        ax.set_axisbelow(True)
        ax.set_xticks(range(1, n_iter + 1))
        ax.set_ylabel("Predicted $log{_2}FC$")
        ax.legend(fontsize='small')

        fig.tight_layout()
        save_path = os.path.join(out_dir, f"genes_plot_{prefix}_{tf}_{cell_line}.pdf")
        fig.savefig(save_path, transparent=True)
        plt.close(fig)
        print(f"      Saved {save_path}")


def _arrow_field(perturbed, original, embedding, n_cpu=40, seed=None, n_grid=25, cache=None):
    """Grid arrow field exactly as plot_perturbation_effect_in_embedding computes it (same private SCENIC+ functions and settings).

    With a seed and a cache path the field is stored (it takes ~10 min), so restyling a figure does not recompute it; the cache is
    reused only if the seed and a checksum of the simulated change match.
    """
    key = None
    if cache is not None and seed is not None:
        key = (seed, perturbed.shape, round(float(np.abs(perturbed.to_numpy() - original.to_numpy()).sum()), 3))
        if os.path.exists(cache):
            with open(cache, "rb") as fh:
                stored = pickle.load(fh)
            if stored.get("key") == key:
                print(f"      [cached arrow field] {os.path.basename(cache)}")
                return stored["field"]
    if seed is not None:
        np.random.seed(seed)
    delta = _project_perturbation_in_embedding(embedding=embedding, perturbed_matrix=perturbed.to_numpy(),
                                               original_matrix=original.to_numpy(), n_cpu=n_cpu)
    grid_xy, uv, mask = _calculate_grid_arrows(embedding=embedding, delta_embedding=delta, offset_frac=0.005,
                                               n_grid_cols=n_grid, n_grid_rows=n_grid, n_neighbors=25, n_cpu=n_cpu)
    dist = np.sqrt((uv ** 2).sum(1))
    scaled = (dist - dist.min()) / (dist.max() - dist.min())
    field = dict(grid_xy=grid_xy, uv=uv, mask=mask, dist=dist, scaled=scaled, n_grid=n_grid)
    if key is not None:
        with open(cache, "wb") as fh:
            pickle.dump({"key": key, "field": field}, fh)
    return field


def _draw_stream(ax, f):
    """SCENIC+'s own streamplot step (copied from plot_perturbation_effect_in_embedding)."""
    n = f["n_grid"]
    uv = f["uv"].copy()
    uv[np.logical_or(~f["mask"], f["scaled"] < 0.15)] = np.nan
    ax.streamplot(f["grid_xy"].reshape(n, n, 2)[:, :, 0], f["grid_xy"].reshape(n, n, 2)[:, :, 1],
                  uv.reshape(n, n, 2)[:, :, 0], uv.reshape(n, n, 2)[:, :, 1],
                  density=3, color=f["scaled"].reshape(n, n), cmap="Greys", zorder=10,
                  norm=matplotlib.colors.Normalize(vmin=0.15, vmax=0.5, clip=True), linewidth=0.5)


def _draw_quiver(ax, f, embedding, min_frac):
    sel = f["mask"] & (f["scaled"] >= min_frac) & (f["dist"] > 0)
    # arrow length: the longest arrow spans about one grid cell
    cell = np.ptp(embedding[:, 0]) / f["n_grid"]
    # arrows above 10% of the longest get a thin white outline so they stand out over dense points; weaker ones stay plain
    # (outlining near-zero arrows only draws white specks). minlength=0: no dot in place of a very short arrow.
    strong = sel & (f["scaled"] >= 0.10)
    for subset, halo in ((sel & ~strong, False), (strong, True)):
        if subset.any():
            q = ax.quiver(f["grid_xy"][subset, 0], f["grid_xy"][subset, 1], f["uv"][subset, 0], f["uv"][subset, 1], color="k", zorder=10,
                          angles="xy", scale_units="xy", scale=f["dist"].max() / (0.9 * cell), width=0.007, minlength=0)
            if halo:
                q.set_path_effects([path_effects.withStroke(linewidth=0.9, foreground="white")])
    return int(sel.sum())


def plot_pca_perturbation(perturbation_over_iters, eRegulon_gene_AUC, scRNA_counts,
                          color_dict, cluster_col, reduction_key, out_dir,
                          direct_metadata, extended_metadata, prefix="kd", n_iter=5):
    """Plot Harmony-corrected PCA embedding with perturbation effect arrows for each TF.

    ARROW_STYLE: "stream" (default, SCENIC+'s own streamplot), "quiver" (individual arrows -> pca_arrows_*.pdf) or "both".
    """
    tfs = list(perturbation_over_iters.keys())
    if not tfs:
        return

    # Use the embedding from eRegulon_gene_AUC (where PCA + Harmony were computed)
    reduction = eRegulon_gene_AUC.obsm[reduction_key]

    def _base_figure(tf):
        fig, ax = plt.subplots(1, 1, figsize=(FIG_SIZE, FIG_SIZE))
        texts = []
        ax.scatter(
            reduction[:, 0], reduction[:, 1],
            s=POINT_SIZE,
            color=[color_dict[line] for line in eRegulon_gene_AUC.obs[cluster_col]],
            rasterized=True,
            linewidths=0,
            alpha=POINT_ALPHA,
        )
        for line in set(eRegulon_gene_AUC.obs[cluster_col]):
            line_bc_idc = np.arange(len(eRegulon_gene_AUC.obs_names))[
                eRegulon_gene_AUC.obs[cluster_col] == line
            ]
            avg_x, avg_y = reduction[line_bc_idc, 0:2].mean(0)
            texts.append(
                ax.text(avg_x, avg_y, line, fontweight="bold",
                        fontsize=5, fontfamily="Arial")
            )
        adjustText.adjust_text(texts, ax=ax)
        return fig, ax

    for tf in tfs:
        field = None
        if ARROW_STYLE != "quiver":
            fig, ax = _base_figure(tf)
            if ARROW_STYLE == "stream":
                # unchanged SCENIC+ call (only the optional seed is new)
                if SEED is not None:
                    np.random.seed(SEED)
                plot_perturbation_effect_in_embedding(
                    perturbed_matrix=perturbation_over_iters[tf][n_iter],
                    original_matrix=perturbation_over_iters[tf][0],
                    embedding=reduction[:, 0:2],
                    AUC_kwargs={},
                    ax=ax,
                    eRegulons=pd.concat([direct_metadata, extended_metadata]),
                    n_cpu=40,
                )
            else:  # "both": one field computation, drawn as streamplot here and as quiver below
                field = _arrow_field(perturbation_over_iters[tf][n_iter], perturbation_over_iters[tf][0],
                                     reduction[:, 0:2], seed=SEED,
                                     cache=os.path.join(out_dir, f'arrow_field_{prefix}_{tf}_seed{SEED}.pkl'))
                _draw_stream(ax, field)

            # Halve arrowhead size on streamplot arrows
            for patch in ax.patches:
                if isinstance(patch, FancyArrowPatch):
                    patch.set_mutation_scale(patch.get_mutation_scale() / 2)

            ax.set_title(tf, fontsize=6, fontfamily="Arial")
            ax.set_axis_off()

            fig.tight_layout(pad=0.3)
            save_path = os.path.join(out_dir, f'pca_lines_{prefix}_{tf}.pdf')
            fig.savefig(save_path, transparent=True, dpi=600)
            plt.close(fig)
            print(f"      Saved {save_path}")

        if ARROW_STYLE in ("quiver", "both"):
            if field is None:
                field = _arrow_field(perturbation_over_iters[tf][n_iter], perturbation_over_iters[tf][0],
                                     reduction[:, 0:2], seed=SEED,
                                     cache=os.path.join(out_dir, f'arrow_field_{prefix}_{tf}_seed{SEED}.pkl'))
            fig, ax = _base_figure(tf)
            n_arrows = _draw_quiver(ax, field, reduction[:, 0:2], ARROW_MIN_FRAC)
            ax.set_title(tf, fontsize=6, fontfamily="Arial")
            print(f"      {tf}: {n_arrows} arrows >= {ARROW_MIN_FRAC:g} of the longest")
            ax.set_axis_off()
            fig.tight_layout(pad=0.3)
            save_path = os.path.join(out_dir, f'pca_arrows_{prefix}_{tf}.pdf')
            fig.savefig(save_path, transparent=True, dpi=600)
            plt.close(fig)
            print(f"      Saved {save_path}")

# ======================================================================================================== notebook code cell 4
config_dirs = sorted(globmod.glob(os.path.join(INPUT_DIR, "config*/")))
if CONFIGS is not None:
    config_dirs = [x for x in config_dirs if os.path.basename(x.rstrip("/")) in CONFIGS]

# Filter to only configs that have the required h5mu file
valid_configs = []
for d in config_dirs:
    h5mu = os.path.join(d, "outs", H5MU_NAME)
    if os.path.exists(h5mu):
        valid_configs.append(d)
    else:
        print(f"  SKIP {os.path.basename(d.rstrip('/'))} — missing {h5mu}")

print(f"\nFound {len(valid_configs)} valid config directories out of {len(config_dirs)} total:")
for d in valid_configs:
    print(f"  {os.path.basename(d.rstrip('/'))}")

# ======================================================================================================== notebook code cell 5
for config_path in valid_configs:
    config_name = os.path.basename(config_path.rstrip('/'))
    print(f"\n{'='*80}")
    print(f"Processing: {config_name}")
    print(f"{'='*80}")

    # ── Paths ──────────────────────────────────────────────────────────────────
    data_dir = os.path.join(config_path, "outs")
    h5mu_fname = os.path.join(data_dir, H5MU_NAME)
    out_dir = os.path.join(config_path, OUT_NAME)
    os.makedirs(out_dir, exist_ok=True)
    regressor_fname = os.path.join(out_dir, 'regressors.pkl')

    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"\n[1/5] Loading data from {h5mu_fname} ...")
    data = load_config_data(h5mu_fname, CLUSTER_COL)
    scRNA_counts = data['scRNA_counts']
    eRegulon_gene_AUC = data['eRegulon_gene_AUC']
    gene_to_TF = data['gene_to_TF']
    gene_tf_df = data['gene_tf_df']
    combined_metadata = data['combined_metadata']
    celltypes = data['celltypes']
    color_dict = make_color_dict(celltypes)
    cell_lines = sorted(celltypes)
    print(f"  Cells: {scRNA_counts.n_obs}, Celltypes: {cell_lines}")

    # ── PCA + Harmony batch correction ────────────────────────────────────────
    sc.pp.pca(eRegulon_gene_AUC)
    print(f"  Running Harmony batch correction (batch_key='{HARMONY_BATCH_KEY}') ...")
    sc.external.pp.harmony_integrate(eRegulon_gene_AUC, key=HARMONY_BATCH_KEY)

    # ── Select genes & train regressors ───────────────────────────────────────
    print(f"\n[2/5] Training / loading regressors ...")
    genes_to_use = gene_tf_df.sort_values("triplet_rank")["Gene"].drop_duplicates()
    genes_to_use = genes_to_use.iloc[:NGENES]
    if FORCE_TARGETS_OF:
        forced = gene_tf_df.loc[gene_tf_df["TF"].isin(FORCE_TARGETS_OF), "Gene"].drop_duplicates()
        print(f"  Forcing the {len(forced)} eRegulon targets of {FORCE_TARGETS_OF} into the regressor set "
              f"({len(set(forced) - set(genes_to_use))} beyond the top {NGENES})")
        genes_to_use = pd.concat([genes_to_use, forced]).drop_duplicates()

    regressors = train_regressors(
        scRNA_counts, genes_to_use, gene_to_TF, out_dir,
        redo=REDO_REGRESSORS, n_jobs=N_JOBS_REGRESSORS
    )
    print(f"  {len(regressors)} regressors available")
    regressor_genes = set(regressors.keys())

    # ── Determine TFs to perturb ──────────────────────────────────────────────
    ko_tfs = KO_TFS if KO_TFS is not None else get_top_tfs(gene_tf_df)
    # Filter to TFs present in expression data
    ko_tfs = [tf for tf in ko_tfs if tf in scRNA_counts.var_names]

    # ══════════════════════════════════════════════════════════════════════════
    # KNOCKOUT
    # ══════════════════════════════════════════════════════════════════════════
    print(f"\n[3/5] Knockout perturbation for TFs: {ko_tfs}")
    perturbation_kd = load_or_run_perturbation(
        out_dir, h5mu_fname, regressor_fname,
        tfs=ko_tfs, set_value=0, n_iter=N_ITER, redo=REDO_PERTURBATIONS
    )

    # ── KO plots ──────────────────────────────────────────────────────────────
    print(f"  Generating KO plots ...")
    for tf in perturbation_kd:
        genes_to_show = get_target_genes(combined_metadata, tf, regressor_genes, n=N_GENES_TO_SHOW)
        if not genes_to_show:
            print(f"    No target genes found for {tf}, skipping plots")
            continue
        print(f"    {tf}: {', '.join(genes_to_show)}")

        # All cell types side-by-side
        try:
            plot_gene_lines(
                perturbation_kd, tf, genes_to_show, eRegulon_gene_AUC,
                cell_lines, CLUSTER_COL, out_dir, prefix="kd", n_iter=N_ITER
            )
        except Exception as e:
            print(f"    WARNING: gene line plot failed for {tf}: {e}")

        # Per cell type
        try:
            genes_ranked = get_target_genes_by_rank(combined_metadata, tf, regressor_genes, n=N_GENES_TO_SHOW)
            plot_gene_lines_per_celltype(
                perturbation_kd, tf, genes_ranked, eRegulon_gene_AUC,
                cell_lines, CLUSTER_COL, out_dir, prefix="kd", n_iter=N_ITER
            )
        except Exception as e:
            print(f"    WARNING: per-celltype plot failed for {tf}: {e}")

    # PCA perturbation plots (KO)
    try:
        plot_pca_perturbation(
            perturbation_kd, eRegulon_gene_AUC, scRNA_counts,
            color_dict, CLUSTER_COL, REDUCTION_KEY, out_dir,
            data['direct_e_regulon_metadata'], data['extended_e_regulon_metadata'],
            prefix="kd", n_iter=N_ITER
        )
    except Exception as e:
        print(f"  WARNING: PCA plot failed for KO: {e}")

    # ══════════════════════════════════════════════════════════════════════════
    # OVEREXPRESSION (optional)
    # ══════════════════════════════════════════════════════════════════════════
    if OE_TFS is not None:
        oe_tfs = ko_tfs if OE_TFS == 'auto' else OE_TFS
        oe_tfs = [tf for tf in oe_tfs if tf in scRNA_counts.var_names]

        print(f"\n[4/5] Overexpression perturbation for TFs: {oe_tfs}")
        perturbation_oe = load_or_run_perturbation(
            out_dir, h5mu_fname, regressor_fname,
            tfs=oe_tfs, set_value=OE_SET_VALUE, n_iter=N_ITER, redo=REDO_PERTURBATIONS
        )

        print(f"  Generating OE plots ...")
        for tf in perturbation_oe:
            genes_to_show = get_target_genes(combined_metadata, tf, regressor_genes, n=N_GENES_TO_SHOW)
            if not genes_to_show:
                continue
            print(f"    {tf}: {', '.join(genes_to_show)}")

            try:
                plot_gene_lines(
                    perturbation_oe, tf, genes_to_show, eRegulon_gene_AUC,
                    cell_lines, CLUSTER_COL, out_dir, prefix="oe", n_iter=N_ITER
                )
            except Exception as e:
                print(f"    WARNING: gene line plot failed for {tf}: {e}")

            try:
                genes_ranked = get_target_genes_by_rank(combined_metadata, tf, regressor_genes, n=N_GENES_TO_SHOW)
                plot_gene_lines_per_celltype(
                    perturbation_oe, tf, genes_ranked, eRegulon_gene_AUC,
                    cell_lines, CLUSTER_COL, out_dir, prefix="oe", n_iter=N_ITER
                )
            except Exception as e:
                print(f"    WARNING: per-celltype plot failed for {tf}: {e}")

        try:
            plot_pca_perturbation(
                perturbation_oe, eRegulon_gene_AUC, scRNA_counts,
                color_dict, CLUSTER_COL, REDUCTION_KEY, out_dir,
                data['direct_e_regulon_metadata'], data['extended_e_regulon_metadata'],
                prefix="oe", n_iter=N_ITER
            )
        except Exception as e:
            print(f"  WARNING: PCA plot failed for OE: {e}")
    else:
        print(f"\n[4/5] Overexpression: skipped (OE_TFS=None)")

    print(f"\n[5/5] Done with {config_name}. Output: {out_dir}")

print(f"\n{'='*80}")
print(f"All {len(valid_configs)} configs processed.")

# ============================================================================================================== panel files
os.makedirs(FIG4_OUT, exist_ok=True)
for prefix in ("kd", "oe"):
    src = os.path.join(os.environ["PERTURB_INPUT_DIR"], "config3", os.environ["PERTURB_OUT_NAME"], f"pca_arrows_{prefix}_MAFB.pdf")
    shutil.copy(src, os.path.join(FIG4_OUT, f"panel_H_pca_arrows_{prefix}_MAFB.pdf"))
    print("saved", os.path.join(FIG4_OUT, f"panel_H_pca_arrows_{prefix}_MAFB.pdf"))
