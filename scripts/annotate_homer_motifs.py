#!/usr/bin/env python3
"""
Copy HOMER de novo motif files to a new directory, renaming the name field
with manual annotations from a CSV produced by parse_homer_results.py.

Rows where manual_annotation is blank or "Exclude" are not copied.

The name field is the second tab-separated field on the motif header line, e.g.:
  >ATGACTCA   1-ATGACTCA,BestGuess:BNC2/MA1928.2/Jaspar(0.976)   5.83 ...
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                replaced with manual_annotation value

Output files are written into OUT_DIR with the same filename as the source (motif{rank}.motif).

Usage:
    python annotate_homer_motifs.py CSV_FILE BASE_DIR OUT_DIR [--dry-run]

Arguments:
    CSV_FILE    CSV with comparison, rank, manual_annotation columns.
    BASE_DIR    Directory containing the motif files directly
                (e.g., peaks_all_homer/homer_dars/ra_pos_vs_neg/homerResults/).
    OUT_DIR     Directory to write annotated motif files into.

Options:
    --dry-run   Print what would be written without creating any files.

Example:
    python annotate_homer_motifs.py homer_results.csv ra_pos_vs_neg/homerResults/ annotated_motifs/
    python annotate_homer_motifs.py homer_results.csv ra_pos_vs_neg/homerResults/ annotated_motifs/ --dry-run
"""

import argparse
import csv
import sys
from pathlib import Path


EXCLUDE_VALUES = {"exclude"}


def annotate_and_write(motif_path, annotation, out_path, dry_run=False):
    lines = motif_path.read_text().splitlines(keepends=True)
    if not lines:
        print(f"WARNING: empty file {motif_path}", file=sys.stderr)
        return

    fields = lines[0].rstrip("\n").split("\t")
    if len(fields) < 2:
        print(f"WARNING: unexpected header in {motif_path}", file=sys.stderr)
        return

    fields[1] = annotation
    lines[0] = "\t".join(fields) + "\n"

    if dry_run:
        print(f"  {motif_path} -> {out_path}  [{annotation}]")
    else:
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text("".join(lines))


def main():
    parser = argparse.ArgumentParser(
        description="Copy annotated HOMER motif files, skipping excluded/blank rows."
    )
    parser.add_argument("csv_file", help="CSV with comparison, rank, manual_annotation columns")
    parser.add_argument("base_dir", help="Directory containing comparison subdirectories")
    parser.add_argument("out_dir", help="Output directory for annotated motif files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be written without creating files")
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    n_written = 0
    n_skipped = 0

    with open(args.csv_file, newline="") as f:
        reader = csv.DictReader(f)
        if "manual_annotation" not in reader.fieldnames:
            print("ERROR: CSV has no 'manual_annotation' column", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            annotation = row["manual_annotation"].strip()
            if not annotation or annotation.lower() in EXCLUDE_VALUES:
                n_skipped += 1
                continue

            motif_path = base_dir / f"motif{row['rank']}.motif"

            if not motif_path.exists():
                print(f"WARNING: not found: {motif_path}", file=sys.stderr)
                n_skipped += 1
                continue

            out_path = out_dir / f"motif{row['rank']}.motif"

            annotate_and_write(motif_path, annotation, out_path, dry_run=args.dry_run)
            n_written += 1

    action = "Would write" if args.dry_run else "Wrote"
    print(f"{action} {n_written} motif file(s), skipped {n_skipped}.", file=sys.stderr)


if __name__ == "__main__":
    main()
