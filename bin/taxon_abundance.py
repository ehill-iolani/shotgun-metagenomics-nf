#!/usr/bin/env python3
"""
Turn KRAKEN2's per-read calls for all samples into
final_report/taxon_abundance.tsv (columns in docs/CONTRACT.md) and a
Kraken2-style report per sample (<sample>.kraken2.report.txt).

KRAKEN2 classifies every sample in one run, so its input is one file of
per-read calls (C/U, read ID, taxid, length) whose read IDs are prefixed
with the sample's position in --samples: `<i>|<read id>`. Kraken2's own
--report would cover all samples at once, so the per-sample reports are
written here instead, in the same format: percent of the sample's reads,
clade reads, direct reads, rank code, taxid, indented name.

Names, ranks and lineages come from the NCBI taxonomy dump the database
was built with (nodes.dmp, names.dmp). Pure standard library, like
read_stats.py.

  taxon_abundance.py --reads kraken2.reads.tsv.gz --samples a b \\
      --nodes db/nodes.dmp --names db/names.dmp \\
      --min-reads 10 --min-rel-abundance 0.0001 \\
      --out taxon_abundance.tsv --report-dir .
"""
import argparse
import csv
import gzip
import os

LINEAGE_RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
# NCBI renamed `superkingdom` to `domain` in 2025, and gave Viruses the
# rank `acellular root`: the superkingdom column takes any of them, so
# viruses still have one
RANK_COLUMN = {r: r for r in LINEAGE_RANKS}
RANK_COLUMN.update({"domain": "superkingdom", "acellular root": "superkingdom"})
COLS = ["sample", "taxid", "name", "rank", "reads", "bases", "rel_abundance", "below_threshold"] + LINEAGE_RANKS

# Kraken2's report rank codes. Any other rank gets the code above it plus
# its depth below that (root's children are R1, a strain under S is S1),
# the same as Kraken2 2.17 writes, `acellular root` included
RANK_CODES = {
    "superkingdom": "D", "domain": "D", "kingdom": "K", "phylum": "P", "class": "C",
    "order": "O", "family": "F", "genus": "G", "species": "S",
}


def read_calls(path, n_samples):
    """Per sample: {taxid: [reads, bases]}; unclassified is taxid 0."""
    counts = [{0: [0, 0]} for _ in range(n_samples)]
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, read_id, taxid, length = line.rstrip("\n").split("\t")
            i = int(read_id.split("|", 1)[0])
            c = counts[i].setdefault(int(taxid), [0, 0])
            c[0] += 1
            c[1] += int(length)
    return counts


def read_nodes(path):
    """{taxid: (parent, rank)} for the whole dump."""
    nodes = {}
    ranks = {}   # one string object per distinct rank: the dump has millions of rows
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t|\t", 3)
            rank = parts[2].split("\t|")[0]
            nodes[int(parts[0])] = (int(parts[1]), ranks.setdefault(rank, rank))
    return nodes


def read_names(path, wanted):
    names = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if parts[3].startswith("scientific name"):
                taxid = int(parts[0])
                if taxid in wanted:
                    names[taxid] = parts[1]
    return names


def lineage(taxid, nodes):
    """[taxid, parent, ..., 1]"""
    path = [taxid]
    while taxid != 1:
        if taxid not in nodes:
            raise SystemExit(f"taxid {taxid} isn't in nodes.dmp: is it the taxonomy the database was built with?")
        taxid = nodes[taxid][0]
        path.append(taxid)
    return path


def write_report(path, counts, nodes, names, lineages):
    total = sum(r for r, _ in counts.values())
    clade = {}
    for taxid, (reads, _) in counts.items():
        if taxid and reads:
            for t in lineages[taxid]:
                clade[t] = clade.get(t, 0) + reads
    children = {}
    for t in clade:
        if t != 1:
            children.setdefault(nodes[t][0], []).append(t)

    with open(path, "w") as fh:
        def line(taxid, code, depth, name):
            n = counts[0][0] if taxid == 0 else clade[taxid]
            direct = counts.get(taxid, [0])[0]
            fh.write(f"{n / total * 100:6.2f}\t{n}\t{direct}\t{code}\t{taxid}\t{'  ' * depth}{name}\n")

        # nothing at all for an empty sample, as Kraken2 does
        if counts[0][0]:
            line(0, "U", 0, "unclassified")
        if 1 not in clade:
            return
        stack = [(1, "R", 0)]
        while stack:
            taxid, code, depth = stack.pop()
            line(taxid, code, depth, names[taxid])
            kids = sorted(children.get(taxid, []), key=lambda t: (-clade[t], t))
            for kid in reversed(kids):
                rank = nodes[kid][1]
                if rank in RANK_CODES:
                    kid_code = RANK_CODES[rank]
                else:
                    kid_code = code[0] + str(int(code[1:] or 0) + 1)
                stack.append((kid, kid_code, depth + 1))


def rows(sample, counts, nodes, names, lineages, min_reads, min_rel):
    total = sum(r for r, _ in counts.values())
    out = []
    for taxid, (reads, bases) in counts.items():
        rel = reads / total if total else 0
        row = {
            "sample": sample, "taxid": taxid, "reads": reads, "bases": bases,
            "rel_abundance": f"{rel:.6g}",
            "below_threshold": "true" if reads < min_reads or rel < min_rel else "false",
        }
        row.update({r: "" for r in LINEAGE_RANKS})
        if taxid == 0:
            row.update(name="Unclassified", rank="unclassified")
        else:
            row.update(name=names[taxid], rank=nodes[taxid][1])
            # nearest first: where two ranks share a column, the lower one wins
            for t in lineages[taxid]:
                col = RANK_COLUMN.get(nodes[t][1])
                if col and not row[col]:
                    row[col] = names[t]
        out.append(row)
    out.sort(key=lambda r: (-r["reads"], r["taxid"]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reads", required=True, help="KRAKEN2's per-read calls, read IDs prefixed <i>|")
    p.add_argument("--samples", required=True, nargs="+", help="sample names, in prefix order")
    p.add_argument("--nodes", required=True)
    p.add_argument("--names", required=True)
    p.add_argument("--min-reads", type=int, default=0)
    p.add_argument("--min-rel-abundance", type=float, default=0)
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--report-dir", required=True, help="where <sample>.kraken2.report.txt go")
    args = p.parse_args()

    counts = read_calls(args.reads, len(args.samples))
    nodes = read_nodes(args.nodes)
    lineages = {t: lineage(t, nodes) for c in counts for t in c if t}
    names = read_names(args.names, {t for path in lineages.values() for t in path})

    table = []
    for sample, c in zip(args.samples, counts):
        write_report(os.path.join(args.report_dir, f"{sample}.kraken2.report.txt"), c, nodes, names, lineages)
        table += rows(sample, c, nodes, names, lineages, args.min_reads, args.min_rel_abundance)

    table.sort(key=lambda r: r["sample"])   # stable: each sample's rows keep their order
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(table)


if __name__ == "__main__":
    main()
