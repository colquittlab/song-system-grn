#!/bin/bash
''''exec micromamba run -n scenicplus python -- "$0" "$@" # '''
"""
Efficient parallel perturbation simulation for specified transcription factors.

Uses multiprocessing.Pool with an initializer so that the large expression
DataFrame and regressors dict are loaded once per worker (shared via
copy-on-write on Linux fork), instead of being serialized per task.

Usage examples:
  # Knockout (set_value=0, default)
  python perturbation_parallel.py --tfs LHX9 UNCX KLF5 FOXO1 HLF NR2F2

  # Overexpression
  python perturbation_parallel.py --tfs LHX9 UNCX --set-value 10

  # Custom paths & parallelism
  python perturbation_parallel.py --tfs LHX9 --n-jobs 4 --n-iter 5 \
      --h5mu /path/to/scplusmdata.h5mu --regressor-pkl /path/to/regressors.pkl \
      --out-dir /path/to/output
"""

import argparse
import glob
import multiprocessing as mp
import os
import pickle
import sys
import time

import h5py
import mudata
import numpy as np
import pandas as pd
from anndata._io.specs.registry import read_elem
from scenicplus.simulation import simulate_perturbation

# ---------------------------------------------------------------------------
# Worker state (set once per worker via Pool initializer, shared copy-on-write)
# ---------------------------------------------------------------------------
_worker_expr_df = None
_worker_regressors = None


def _init_worker(expr_df, regressors):
    """Initializer called once per worker process."""
    global _worker_expr_df, _worker_regressors
    _worker_expr_df = expr_df
    _worker_regressors = regressors


def _run_perturbation(args):
    """Worker function: run perturbation for a single TF.

    Reads shared state from module-level globals set by _init_worker,
    avoiding repeated serialization of large objects.
    """
    tf_name, out_dir, set_value, n_iter = args
    sv_tag = f"sv{set_value:g}"
    fname_prefix = os.path.join(out_dir, f"perturb_{tf_name}_{sv_tag}")
    existing = glob.glob(f"{fname_prefix}_*.parquet")

    if existing:
        # Load cached results
        result = {}
        for i in range(n_iter + 1):
            pq = f"{fname_prefix}_{i}.parquet"
            if os.path.exists(pq):
                result[i] = pd.read_parquet(pq, engine="fastparquet")
        print(f"  [cached]  {tf_name}  ({len(result)} iterations loaded)")
        return tf_name, result

    print(f"  [running] {tf_name}  (set_value={set_value}, n_iter={n_iter})")
    t0 = time.time()

    result = simulate_perturbation(
        df_EXP=_worker_expr_df,
        perturbation={tf_name: set_value},
        keep_intermediate=True,
        n_iter=n_iter,
        regressors=_worker_regressors,
    )

    # Save each iteration as parquet
    for i in result:
        result[i].to_parquet(f"{fname_prefix}_{i}.parquet")

    elapsed = time.time() - t0
    print(f"  [done]    {tf_name}  ({elapsed:.1f}s)")
    return tf_name, result


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------
def mudata_to_pandas(h5mu_fname, group_path):
    """Read an HDF5 group from a .h5mu file into a pandas DataFrame."""
    with h5py.File(h5mu_fname, "r") as f:
        group = f[group_path]
        data = {}
        for key, item in group.items():
            if isinstance(item, h5py.Dataset):
                data[key] = read_elem(item)
    return pd.DataFrame(data)


def load_data(h5mu_fname, regressor_pkl):
    """Load expression data and pre-trained regressors.

    Returns
    -------
    expr_df : pd.DataFrame
        Cells × genes expression matrix.
    regressors : dict
        Pre-trained gene expression regressors.
    """
    print(f"Loading expression data from {h5mu_fname} ...")
    scRNA_counts = mudata.read(f"{h5mu_fname}/scRNA_counts")
    expr_df = scRNA_counts.to_df()
    print(f"  expression matrix: {expr_df.shape[0]} cells × {expr_df.shape[1]} genes")

    print(f"Loading regressors from {regressor_pkl} ...")
    with open(regressor_pkl, "rb") as fh:
        regressors = pickle.load(fh)
    print(f"  {len(regressors)} gene regressors loaded")

    return expr_df, regressors


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Run SCENIC+ perturbation simulations in parallel.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--tfs",
        nargs="+",
        required=True,
        help="Transcription factor names to perturb (e.g. LHX9 UNCX KLF5).",
    )
    parser.add_argument(
        "--set-value",
        type=float,
        default=0,
        help="Expression value to set for the perturbed TF. "
        "0 = knockout (default), >0 = overexpression.",
    )
    parser.add_argument(
        "--n-iter",
        type=int,
        default=5,
        help="Number of simulation iterations (default: 5).",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="Number of parallel workers. Defaults to min(len(tfs), cpu_count).",
    )
    parser.add_argument(
        "--h5mu",
        required=True,
        help="Path to the .h5mu file.",
    )
    parser.add_argument(
        "--regressor-pkl",
        default=None,
        help="Path to regressors pickle. Defaults to <out-dir>/regressors.pkl.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory for parquet results. "
        "Defaults to <h5mu parent>/perturbation.",
    )

    args = parser.parse_args()

    # Resolve paths
    if args.out_dir is None:
        args.out_dir = os.path.join(os.path.dirname(args.h5mu), "perturbation")
    os.makedirs(args.out_dir, exist_ok=True)

    if args.regressor_pkl is None:
        args.regressor_pkl = os.path.join(args.out_dir, "regressors.pkl")

    # Load shared data
    expr_df, regressors = load_data(args.h5mu, args.regressor_pkl)

    # Validate TFs
    tfs_valid = [tf for tf in args.tfs if tf in expr_df.columns]
    tfs_missing = set(args.tfs) - set(tfs_valid)
    if tfs_missing:
        print(f"WARNING: TFs not found in expression columns, skipping: {tfs_missing}")
    if not tfs_valid:
        print("ERROR: No valid TFs to process.")
        sys.exit(1)

    print(f"\nTFs to process: {tfs_valid}")
    print(f"  set_value={args.set_value}  n_iter={args.n_iter}")

    # Determine worker count
    n_jobs = args.n_jobs or min(len(tfs_valid), mp.cpu_count())
    n_jobs = min(n_jobs, len(tfs_valid))  # no point exceeding TF count
    print(f"  n_jobs={n_jobs}\n")

    # Build task list
    tasks = [
        (tf, args.out_dir, args.set_value, args.n_iter) for tf in tfs_valid
    ]

    # Run in parallel using Pool + initializer
    # The initializer passes expr_df and regressors to each worker ONCE,
    # avoiding per-task serialization overhead (major speedup).
    t_start = time.time()

    if n_jobs == 1:
        # Skip pool overhead for single TF
        _init_worker(expr_df, regressors)
        results = [_run_perturbation(t) for t in tasks]
    else:
        with mp.Pool(
            processes=n_jobs,
            initializer=_init_worker,
            initargs=(expr_df, regressors),
        ) as pool:
            results = pool.map(_run_perturbation, tasks)

    # Collect into dict
    perturbation_over_iters = {tf: res for tf, res in results}

    elapsed = time.time() - t_start
    print(f"\nAll done. {len(perturbation_over_iters)} TFs processed in {elapsed:.1f}s.")
    print(f"Results saved to: {args.out_dir}")

    return perturbation_over_iters


if __name__ == "__main__":
    main()
