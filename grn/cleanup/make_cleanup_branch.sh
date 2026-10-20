#!/usr/bin/env bash
# Build (locally) the branch that removes development material from main. Dry run by default; never pushes.
#
#   ./make_cleanup_branch.sh                       # dry run: print what would be removed
#   APPLY=1 ./make_cleanup_branch.sh               # create $BRANCH from $BASE, `git rm` the files, commit (local only)
#   CATS=DEVELOPMENT,DECIDE APPLY=1 ./make_cleanup_branch.sh   # also remove files you have not yet resolved (see DECISIONS in README.md)
#
# What is removed is whatever inventory.csv lists in the chosen categories. Edit the CSV first: change a row's `category` to
# PUBLICATION to keep that file, to DEVELOPMENT to remove it. Paths that no longer exist on $BASE are skipped.
#
# Env:  BASE (default origin/main)   BRANCH (default main-cleanup)   CATS (default DEVELOPMENT)   APPLY (1 to act)
#
# After applying: review `git show --stat`, push the branch yourself and open the PR against main. See README.md for what to do with
# grn-dev afterwards (do NOT merge main into grn-dev wholesale: it would delete these files there too; use the revert recipe).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
inv="$here/inventory.csv"
BASE="${BASE:-origin/main}"
BRANCH="${BRANCH:-main-cleanup}"
CATS="${CATS:-DEVELOPMENT}"
[ -s "$inv" ] || { echo "missing $inv" >&2; exit 1; }

list="$(python3 - "$inv" "$CATS" <<'PY'
import csv, sys
cats = set(sys.argv[2].split(","))
for r in csv.DictReader(open(sys.argv[1])):
    if r["category"] in cats:
        print(r["path"])
PY
)"
n="$(printf '%s\n' "$list" | grep -c . || true)"
echo "categories to remove: $CATS -> $n files (from $inv, snapshot of the ref it was built on)"
printf '%s\n' "$list" | sed 's#^\(grn/[^/]*/[^/]*\).*#\1#' | sort | uniq -c | sort -rn | head -25

if [ "${APPLY:-0}" != "1" ]; then
    echo; echo "dry run only. APPLY=1 to create branch '$BRANCH' from '$BASE' and remove these files locally (nothing is pushed)."
    exit 0
fi

toplevel="$(git rev-parse --show-toplevel)"
cd "$toplevel"
git diff --quiet && git diff --cached --quiet || { echo "working tree has uncommitted changes: commit or stash first" >&2; exit 1; }
git rev-parse --verify "$BASE" >/dev/null
git rev-parse --verify "$BRANCH" >/dev/null 2>&1 && { echo "branch $BRANCH already exists" >&2; exit 1; }
git switch -c "$BRANCH" "$BASE"
printf '%s\n' "$list" | xargs -r git rm -q -r --ignore-unmatch --
msg="Move GRN development material to grn-dev ($CATS)

Removes $n files listed in grn/cleanup/inventory.csv (categories: $CATS). Everything removed remains on grn-dev and in history."
git commit -q -m "$msg"
echo; echo "created $BRANCH:"; git show --stat --format='%h %s' HEAD | tail -3
echo "next: review, then push $BRANCH and open the PR against main (this script does not push)."
