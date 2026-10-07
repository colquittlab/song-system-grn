# cisTarget database for ra-arco-hvc-nc_hybrid, built on prism

The database is built on the 499,348 consensus regions of the current cisTopic object (see
`../create_cistarget_db_hybrid.sh` for why the old database does not match them). Scoring ~10,249 motifs over those regions is about
19 h on 40 cores, longer than a single local background job may run, so it is split into parts that run in parallel:

1. `score_part.sh` -- one part of the motif list (`create_cistarget_motif_databases.py -p PART NPARTS`); an array job, 6 parts of ~3 h.
2. `combine_and_rank.sh` -- combines the partial score databases, makes the rankings, removes the partial files.

`submit_cistarget_build.sh` submits both with a dependency. Tested locally on 200 regions x 6 motifs in 3 parts: the combined scores
are identical to a single-pass run.

## What lives where

The database directory holds the inputs and receives the outputs, so nothing needs copying back:

```
/private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid/     (= CISTARGET_HYBRID_DIR in make_configs.py)
    lonStrDom2_1kb_bg_padding.fa    padded sequences of the consensus regions (1.3 GB)
    consensus_regions.bed           the regions the database is built on
    motifs.txt                      the 10,249 Cluster-Buster motif files to score
    singletons/                     those motif files (42 MB)
    tools/                          create_cisTarget_databases scripts + a static cbust binary (3 MB)
    prism/                          these scripts
```

Built from the local copy at `/hdd/jupyter/brad/scenicplus/cistarget/ra-arco-hvc-nc_hybrid/`.

## Run

On prism (the `scenicplus` conda env needs only numpy, pandas and pyarrow; `cbust` is bundled):

```
rclone copy -vP --exclude build.log lark:/hdd/jupyter/brad/scenicplus/cistarget/ra-arco-hvc-nc_hybrid \
    /private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid
cd /private/groups/colquittlab/scenicplus/cistarget/ra-arco-hvc-nc_hybrid
DRY_RUN=1 ./prism/submit_cistarget_build.sh     # check the commands
./prism/submit_cistarget_build.sh               # 6 parts; pass another number to change it
squeue -u $USER
```

Output: `ra-arco-hvc-nc_hybrid.regions_vs_motifs.rankings.feather` and `.regions_vs_motifs.scores.feather` (read by SCENIC+),
plus `.motifs_vs_regions.scores.feather`. Check at the end that the rankings file has 10,249 rows and 499,348 columns.
A failed or requeued scoring task resumes: finished parts are skipped.
