#!/usr/bin/env python3
"""
Extract Best Match fields and known motif matches from homerResults.html files.

Usage:
    # Single file
    python parse_homer_results.py homerResults.html

    # Multiple files
    python parse_homer_results.py dir1/homerResults.html dir2/homerResults.html

    # Recursively find all homerResults.html under a directory
    python parse_homer_results.py --dir /path/to/homer/results/

    # Write to file instead of stdout
    python parse_homer_results.py --dir . -o results.csv

    # Limit number of known motif match columns (default: 10)
    python parse_homer_results.py --dir . --n-matches 5 -o results.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from bs4 import BeautifulSoup

MAX_KNOWN_MATCHES = 10


def parse_info_html(info_path):
    """
    Parse a motifN.info.html file and return list of known motif matches.

    Returns list of dicts with keys: name, score, orientation (up to MAX_KNOWN_MATCHES).
    """
    if not info_path.exists():
        return []

    with open(info_path) as f:
        raw = f.read()

    # HOMER writes malformed </TD</TR> (missing >) — fix before parsing
    raw = raw.replace("</TD</TR>", "</TD></TR>")
    soup = BeautifulSoup(raw, "html.parser")

    matches = []
    # Each known motif entry is headed by an <H4> tag with the motif name.
    for h4 in soup.find_all("h4"):
        name = h4.get_text(strip=True)
        score = ""
        orientation = ""
        # The details table immediately follows inside the same outer <TD>
        detail_table = h4.find_next("table")
        if detail_table:
            for row in detail_table.find_all("tr"):
                cells = row.find_all("td")
                if len(cells) < 2:
                    continue
                label = cells[0].get_text(strip=True).rstrip(":")
                value = cells[1].get_text(strip=True)
                if label == "Score":
                    score = value
                elif label == "Orientation":
                    orientation = value
        matches.append({"name": name, "score": score, "orientation": orientation})

    return matches


def parse_homer_html(html_path, n_matches):
    """Parse a homerResults.html file and return list of row dicts."""
    html_path = Path(html_path)
    comparison = html_path.parent.name
    info_dir = html_path.parent / "homerResults"

    with open(html_path) as f:
        soup = BeautifulSoup(f, "html.parser")

    table = soup.find("table")
    if table is None:
        return []

    rows = table.find_all("tr")
    records = []
    for row in rows[1:]:  # skip header
        cols = row.find_all("td")
        if len(cols) < 8:
            continue

        rank = cols[0].get_text(strip=True)
        pvalue = cols[2].get_text(strip=True)
        log_pvalue = cols[3].get_text(strip=True)
        pct_targets = cols[4].get_text(strip=True)
        pct_bg = cols[5].get_text(strip=True)
        std = cols[6].get_text(strip=True)

        # Best Match/Details: text before the <BR/> tag
        details_td = cols[7]
        br = details_td.find("br")
        if br:
            best_match = "".join(
                str(s) for s in br.previous_siblings
                if isinstance(s, str)
            ).strip()
        else:
            best_match = details_td.get_text(strip=True)

        record = {
            "comparison": comparison,
            "rank": rank,
            "best_match": best_match,
            "pvalue": pvalue,
            "log_pvalue": log_pvalue,
            "pct_targets": pct_targets,
            "pct_background": pct_bg,
            "std_bg_std": std,
        }

        # Parse the linked motifN.info.html
        info_path = info_dir / f"motif{rank}.info.html"
        known = parse_info_html(info_path)
        for i in range(1, n_matches + 1):
            if i <= len(known):
                record[f"known_match_{i}"] = known[i - 1]["name"]
                record[f"known_match_{i}_score"] = known[i - 1]["score"]
                record[f"known_match_{i}_orientation"] = known[i - 1]["orientation"]
            else:
                record[f"known_match_{i}"] = ""
                record[f"known_match_{i}_score"] = ""
                record[f"known_match_{i}_orientation"] = ""

        records.append(record)

    return records


def build_fieldnames(n_matches):
    base = ["comparison", "rank", "best_match", "pvalue", "log_pvalue",
            "pct_targets", "pct_background", "std_bg_std"]
    for i in range(1, n_matches + 1):
        base += [f"known_match_{i}", f"known_match_{i}_score", f"known_match_{i}_orientation"]
    return base


def find_html_files(directory):
    return sorted(Path(directory).resolve().rglob("homerResults.html"))


def main():
    parser = argparse.ArgumentParser(
        description="Extract best match and known motif info from homerResults.html to CSV."
    )
    parser.add_argument("files", nargs="*", help="homerResults.html file(s)")
    parser.add_argument("--dir", help="Directory to search recursively for homerResults.html")
    parser.add_argument("-o", "--output", help="Output CSV file (default: stdout)")
    parser.add_argument("--n-matches", type=int, default=MAX_KNOWN_MATCHES,
                        help=f"Number of known motif match columns to include (default: {MAX_KNOWN_MATCHES})")
    args = parser.parse_args()

    html_files = []
    if args.dir:
        html_files = find_html_files(args.dir)
        if not html_files:
            print(f"No homerResults.html found under {args.dir}", file=sys.stderr)
            sys.exit(1)
    elif args.files:
        html_files = [Path(f).resolve() for f in args.files]
    else:
        parser.print_help()
        sys.exit(1)

    all_records = []
    for f in html_files:
        records = parse_homer_html(f, args.n_matches)
        if not records:
            print(f"Warning: no data parsed from {f}", file=sys.stderr)
        all_records.extend(records)

    fieldnames = build_fieldnames(args.n_matches)

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_records)
    finally:
        if args.output:
            out.close()
            print(f"Wrote {len(all_records)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
