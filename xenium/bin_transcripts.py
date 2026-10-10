"""Bin Xenium transcripts into a 2D density grid without loading the file.

nanoparquet segfaults on these 2 GB transcript-metadata files, so stream with
pyarrow instead: iterate row-group batches, keep three columns, accumulate a
fixed-extent 2D histogram. Peak memory is one batch, not the file.

Uses proseg's transcript-metadata.parquet rather than the raw Xenium
transcripts.parquet so the coordinates are already in the same frame as the
cell centroids, and reads observed_x/observed_y (the detections as measured)
rather than x/y (proseg's model-repositioned coordinates) -- the point of a
cell-free view is that no segmentation model has touched it.
"""
import sys, numpy as np, pyarrow.parquet as pq

path, out_csv, bin_um, qv_min = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
pf = pq.ParquetFile(path)
cols = ["observed_x", "observed_y", "qv"]

# pass 1: extents (cheap -- min/max only)
xmin = ymin = np.inf; xmax = ymax = -np.inf; n_tot = 0
for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
    x = b.column("observed_x").to_numpy(zero_copy_only=False)
    y = b.column("observed_y").to_numpy(zero_copy_only=False)
    xmin = min(xmin, x.min()); xmax = max(xmax, x.max())
    ymin = min(ymin, y.min()); ymax = max(ymax, y.max())
    n_tot += len(x)
nx = int(np.ceil((xmax - xmin) / bin_um)); ny = int(np.ceil((ymax - ymin) / bin_um))
print(f"transcripts: {n_tot:,}   extent {xmax-xmin:.0f} x {ymax-ymin:.0f} um   grid {nx} x {ny}", flush=True)

# pass 2: accumulate
H = np.zeros((nx, ny), dtype=np.int64); n_kept = 0
for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
    x = b.column("observed_x").to_numpy(zero_copy_only=False)
    y = b.column("observed_y").to_numpy(zero_copy_only=False)
    q = b.column("qv").to_numpy(zero_copy_only=False)
    m = q >= qv_min
    x, y = x[m], y[m]
    ix = np.clip(((x - xmin) / bin_um).astype(np.int64), 0, nx - 1)
    iy = np.clip(((y - ymin) / bin_um).astype(np.int64), 0, ny - 1)
    np.add.at(H, (ix, iy), 1)
    n_kept += len(x)
print(f"kept qv>={qv_min}: {n_kept:,} ({100*n_kept/n_tot:.1f}%)", flush=True)

ix, iy = np.nonzero(H)
with open(out_csv, "w") as f:
    f.write("x,y,n\n")
    for i, j in zip(ix, iy):
        f.write(f"{xmin + (i + 0.5) * bin_um:.1f},{ymin + (j + 0.5) * bin_um:.1f},{H[i, j]}\n")
print(f"wrote {len(ix):,} occupied bins -> {out_csv}", flush=True)
