#!/usr/bin/env python3
"""
Turn each sample's AMR gene hits into final_report/amr_genes.tsv and
<sample>.amr_reads.tsv, naming the organism each gene was found in from
the same reads' Kraken2 calls.

Hits are minimap2 PAF (--format paf) or ABRicate's table (--format
abricate, run with --noheader), one file per sample named
<sample>.amr_hits.<ext>, against the sequences AMR_DB laid out: PAF target
names are amr_genes.tsv's `id`, ABRicate's are matched back to it by gene
and accession. A hit passes at --min-identity and --min-coverage (percent
of the gene's length). Where passing hits overlap on a read by more than
half the shorter one, only the best is kept (identity x coverage), so a
read counts once per gene however many near-identical alleles it matches.

KRAKEN2's per-read calls give each sample's total reads and bases (the
normalising denominators: every host-removed read is in them, classified
or not) and the taxon each AMR read was called as. Read IDs there are
prefixed with the sample's position in --samples (`<i>|<read id>`).
Pure standard library, like taxon_abundance.py.

  amr_report.py --hits a.amr_hits.paf b.amr_hits.paf --format paf \\
      --genes amr_genes.tsv --samples a b --kraken2 kraken2.reads.tsv.gz \\
      --nodes db/nodes.dmp --names db/names.dmp --intrinsic intrinsic_genes.tsv \\
      --min-identity 90 --min-coverage 80 --out amr_genes.tsv --reads-dir .
"""
import argparse
import collections
import csv
import gzip
import os
import re

GENE_COLS = ["sample", "gene_family", "closest_allele", "drug_class", "drug_subclass", "reads",
             "reads_per_million", "copies_per_gb", "call", "host_genera", "intrinsic"]
READ_COLS = ["read_id", "allele", "gene_family", "identity", "gene_coverage", "aligned_bases",
             "kraken2_taxid", "kraken2_name"]
Hit = collections.namedtuple("Hit", "read gene_id read_start read_end identity coverage aligned_bases")


def read_table(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def read_paf(path, genes):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        # primary alignments only: secondaries are other alleles of the
        # same stretch of read
        if "tp:A:P" not in f[12:]:
            continue
        gene_len, g_start, g_end = int(f[6]), int(f[7]), int(f[8])
        yield Hit(f[0], f[5], int(f[2]), int(f[3]),
                  int(f[9]) / int(f[10]) * 100, (g_end - g_start) / gene_len * 100, g_end - g_start)


def read_abricate(path, genes):
    ids = {(g["gene"], g["accession"]): g["id"] for g in genes.values()}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        # FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS
        # %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE
        g_start, g_end = re.match(r"(\d+)-(\d+)/", f[6]).groups()
        yield Hit(f[1], ids[(f[5], f[12])], int(f[2]) - 1, int(f[3]),
                  float(f[10]), float(f[9]), int(g_end) - int(g_start) + 1)


def best_hits(hits, min_identity, min_coverage):
    """{read: [Hit]}, passing and non-overlapping, best first."""
    by_read = collections.defaultdict(list)
    for h in hits:
        if h.identity >= min_identity and h.coverage >= min_coverage:
            by_read[h.read].append(h)
    kept = {}
    for read, hs in by_read.items():
        hs.sort(key=lambda h: (-h.identity * h.coverage, int(h.gene_id)))
        keep = []
        for h in hs:
            if all(min(h.read_end, k.read_end) - max(h.read_start, k.read_start)
                   <= min(h.read_end - h.read_start, k.read_end - k.read_start) / 2 for k in keep):
                keep.append(h)
        kept[read] = keep
    return kept


def read_kraken2(path, wanted):
    """Per-sample [reads, bases], and {prefixed read ID: taxid} for wanted reads."""
    totals = collections.defaultdict(lambda: [0, 0])
    calls = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            _, read_id, taxid, length = line.rstrip("\n").split("\t")
            t = totals[int(read_id.split("|", 1)[0])]
            t[0] += 1
            t[1] += int(length)
            if read_id in wanted:
                calls[read_id] = int(taxid)
    return totals, calls


def read_nodes(path):
    nodes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t|\t", 3)
            nodes[int(parts[0])] = (int(parts[1]), parts[2].split("\t|")[0])
    return nodes


def read_names(path, wanted):
    names = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if parts[3].startswith("scientific name") and int(parts[0]) in wanted:
                names[int(parts[0])] = parts[1]
    return names


def genus_of(taxid, nodes):
    """taxid of the genus at or above taxid; None if the call is above genus."""
    while taxid not in (0, 1):
        if taxid not in nodes:
            raise SystemExit(f"taxid {taxid} isn't in nodes.dmp: is it the taxonomy the database was built with?")
        parent, rank = nodes[taxid]
        if rank == "genus":
            return taxid
        taxid = parent
    return None


def fmt(x):
    return f"{x:.6g}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hits", required=True, nargs="+", help="<sample>.amr_hits.<ext>, one per sample")
    p.add_argument("--format", required=True, choices=["paf", "abricate"])
    p.add_argument("--genes", required=True, help="AMR_DB's annotation table")
    p.add_argument("--samples", required=True, nargs="+", help="sample names, in KRAKEN2's prefix order")
    p.add_argument("--kraken2", required=True, help="KRAKEN2's per-read calls, read IDs prefixed <i>|")
    p.add_argument("--nodes", required=True)
    p.add_argument("--names", required=True)
    p.add_argument("--intrinsic", required=True, help="gene, genus, note")
    p.add_argument("--min-identity", type=float, required=True, help="percent")
    p.add_argument("--min-coverage", type=float, required=True, help="percent of the gene's length")
    p.add_argument("--out", required=True)
    p.add_argument("--reads-dir", required=True, help="where <sample>.amr_reads.tsv go")
    args = p.parse_args()

    genes = {g["id"]: g for g in read_table(args.genes)}
    intrinsic = collections.defaultdict(set)
    for r in read_table(args.intrinsic):
        intrinsic[r["gene"]].add(r["genus"])
    index = {s: i for i, s in enumerate(args.samples)}
    parse = read_paf if args.format == "paf" else read_abricate

    per_sample = {}
    for path in args.hits:
        sample = re.sub(r"\.amr_hits\.\w+$", "", os.path.basename(path))
        if sample not in index:
            raise SystemExit(f"{path}: sample '{sample}' isn't one of --samples")
        per_sample[sample] = best_hits(parse(path, genes), args.min_identity, args.min_coverage)

    wanted = {f"{index[s]}|{read}" for s, kept in per_sample.items() for read in kept}
    totals, calls = read_kraken2(args.kraken2, wanted)
    if wanted - calls.keys():
        raise SystemExit(f"{len(wanted - calls.keys())} AMR reads aren't in {args.kraken2}, e.g. {sorted(wanted - calls.keys())[0]}")
    nodes = read_nodes(args.nodes)
    genus = {t: genus_of(t, nodes) for t in set(calls.values())}
    names = read_names(args.names, set(calls.values()) | {g for g in genus.values() if g})

    table = []
    for sample in args.samples:
        kept = per_sample.get(sample, {})
        n_reads, n_bases = totals[index[sample]]
        read_rows = []
        families = collections.defaultdict(lambda: {"reads": set(), "alleles": collections.Counter(), "annotation": {}, "copies": 0.0, "genera": collections.Counter()})
        for read, hs in kept.items():
            taxid = calls[f"{index[sample]}|{read}"]
            for h in hs:
                g = genes[h.gene_id]
                fam = families[g["gene_family"]]
                if read not in fam["reads"]:
                    fam["reads"].add(read)
                    fam["alleles"][g["allele"]] += 1
                    fam["annotation"][g["allele"]] = g
                    if genus[taxid]:
                        fam["genera"][names[genus[taxid]]] += 1
                fam["copies"] += h.aligned_bases / int(g["length"])
                read_rows.append({
                    "read_id": read, "allele": g["allele"], "gene_family": g["gene_family"],
                    "identity": f"{h.identity:.2f}", "gene_coverage": f"{h.coverage:.2f}",
                    "aligned_bases": h.aligned_bases, "kraken2_taxid": taxid,
                    "kraken2_name": names[taxid] if taxid else "Unclassified",
                })

        read_rows.sort(key=lambda r: (r["gene_family"], r["read_id"], r["allele"]))
        with open(os.path.join(args.reads_dir, f"{sample}.amr_reads.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=READ_COLS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(read_rows)

        rows = []
        for family, fam in families.items():
            allele = min(fam["alleles"], key=lambda a: (-fam["alleles"][a], a))
            hosts = [g for g, n in sorted(fam["genera"].items(), key=lambda x: (-x[1], x[0])) if n >= 2]
            flagged = intrinsic[family] | intrinsic[allele]
            reads = len(fam["reads"])
            rows.append({
                "sample": sample, "gene_family": family, "closest_allele": allele,
                "drug_class": fam["annotation"][allele]["drug_class"],
                "drug_subclass": fam["annotation"][allele]["drug_subclass"],
                "reads": reads,
                "reads_per_million": fmt(reads / n_reads * 1e6),
                "copies_per_gb": fmt(fam["copies"] / (n_bases / 1e9)),
                "call": "detected" if reads >= 2 else "possible",
                "host_genera": ";".join(hosts),
                "intrinsic": ";".join(g for g in hosts if g in flagged),
            })
        rows.sort(key=lambda r: (-r["reads"], r["gene_family"]))
        table += rows

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=GENE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(table)


if __name__ == "__main__":
    main()
