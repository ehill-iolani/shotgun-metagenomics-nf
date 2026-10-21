#!/usr/bin/env python3
"""
Combine Deacon filter summaries (`deacon filter -s <sample>.deacon.json`)
into one row per sample: final_report/host_removal.tsv, columns in
docs/CONTRACT.md. The sample name is the JSON's filename up to
`.deacon.json`.

Pure standard library, like read_stats.py.

  host_removal_report.py --out host_removal.tsv *.deacon.json
"""
import argparse
import csv
import json
import os

COLS = [
    "sample", "reads_in", "reads_removed", "pct_reads_removed",
    "bases_in", "bases_removed", "pct_bases_removed", "host_index",
]
SUFFIX = ".deacon.json"


def pct(part, whole):
    return "NA" if not whole else f"{part / whole * 100:.2f}"


def row(path):
    name = os.path.basename(path)
    if not name.endswith(SUFFIX):
        raise SystemExit(f"{path}: expected a <sample>{SUFFIX} file")
    with open(path) as fh:
        s = json.load(fh)
    return {
        "sample": name[: -len(SUFFIX)],
        "reads_in": s["seqs_in"],
        "reads_removed": s["seqs_removed"],
        "pct_reads_removed": pct(s["seqs_removed"], s["seqs_in"]),
        "bases_in": s["bp_in"],
        "bases_removed": s["bp_removed"],
        "pct_bases_removed": pct(s["bp_removed"], s["bp_in"]),
        "host_index": os.path.basename(s["index"]),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("summaries", nargs="+", help="per-sample Deacon summary JSONs")
    p.add_argument("--out", required=True, help="output TSV")
    args = p.parse_args()

    rows = sorted((row(path) for path in args.summaries), key=lambda r: r["sample"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
