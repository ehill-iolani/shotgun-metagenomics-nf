#!/usr/bin/env python3
"""
Writes the synthetic fixtures -profile test runs on (tests/data/), so they
can be regenerated, or extended for a new stage, reproducibly.

Random sequence, not real genomes: the test profile checks pipeline
mechanics, not biology. A "host" genome stands in for human (the test
profile builds its Deacon index from it, instead of downloading the 3.3 GB
panhuman one), and three "microbe" genomes are what the reads are meant to
be classified as.

Reads are simulated long reads (600-5000 bp, 1% substitutions, Q20), plus a
few too-short and low-quality ones that chopper should drop. Output is
deterministic (fixed seed, gzip mtime 0), so rerunning without changes
leaves git with nothing to commit.

The microbes get a made-up NCBI-style taxonomy (kraken2_taxonomy/), which
the test profile builds its Kraken2 database from instead of downloading
PlusPF. It uses NCBI's current rank names (`cellular root`, `domain`) and
leaves gaps -- no kingdom anywhere, no class above microbe_c, a strain
below a species -- so the lineage columns of taxon_abundance.tsv get
exercised.

  python3 tests/make_test_data.py
"""
import gzip
import os
import random

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
rng = random.Random(20261021)

GENOMES = {"host": 200_000, "microbe_a": 100_000, "microbe_b": 100_000, "microbe_c": 100_000}
# reads per source genome; the two samples differ so abundances differ
SAMPLES = {
    "sample_a": {"host": 30, "microbe_a": 40, "microbe_b": 15, "microbe_c": 5},
    "sample_b": {"host": 5, "microbe_a": 5, "microbe_b": 20, "microbe_c": 30},
}
# (taxid, parent, rank, name); genome names in GENOME_TAXIDS
TAXONOMY = [
    (1, 1, "no rank", "root"),
    (131567, 1, "cellular root", "cellular organisms"),
    (100, 131567, "domain", "Fakebacteria"),
    (200, 100, "phylum", "Mockota"),
    (300, 200, "class", "Mockia"),
    (400, 300, "order", "Mockales"),
    (500, 400, "family", "Mockaceae"),
    (600, 500, "genus", "Mockus"),
    (1001, 600, "species", "Mockus alpha"),
    (1002, 600, "species", "Mockus beta"),
    (210, 100, "phylum", "Testota"),
    (410, 210, "order", "Testales"),
    (510, 410, "family", "Testaceae"),
    (610, 510, "genus", "Testella"),
    (1003, 610, "species", "Testella gamma"),
    (2003, 1003, "strain", "Testella gamma C1"),
]
GENOME_TAXIDS = {"microbe_a": 1001, "microbe_b": 1002, "microbe_c": 2003}
N_SHORT = 5    # < min_len 500
N_LOW_Q = 5    # mean Q7, < min_qual 10
GOOD_Q, LOW_Q = 20, 7


def write_gz(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as fh:
        fh.write(text.encode())


def random_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate=0.01):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def sample_read(genome, lo, hi):
    n = rng.randint(lo, hi)
    start = rng.randint(0, len(genome) - n)
    read = mutate(genome[start:start + n])
    return revcomp(read) if rng.random() < 0.5 else read


def fastq_record(name, seq, q):
    return f"@{name}\n{seq}\n+\n{chr(q + 33) * len(seq)}\n"


def main():
    genomes = {name: random_seq(n) for name, n in GENOMES.items()}
    write_gz(os.path.join(OUT, "host.fasta.gz"), f">host\n{genomes['host']}\n")
    for name in GENOMES:
        if name != "host":
            write_gz(os.path.join(OUT, "genomes", f"{name}.fasta.gz"), f">{name}\n{genomes[name]}\n")

    for sample, counts in SAMPLES.items():
        records = []
        for source, n in counts.items():
            for i in range(n):
                records.append(fastq_record(f"{sample}_{source}_{i}", sample_read(genomes[source], 600, 5000), GOOD_Q))
        for i in range(N_SHORT):
            records.append(fastq_record(f"{sample}_short_{i}", sample_read(genomes["microbe_a"], 100, 400), GOOD_Q))
        for i in range(N_LOW_Q):
            records.append(fastq_record(f"{sample}_lowq_{i}", sample_read(genomes["microbe_a"], 600, 5000), LOW_Q))
        rng.shuffle(records)
        # sample_a split into two part-files, like a barcode folder from MinKNOW
        parts = [records[: len(records) // 2], records[len(records) // 2:]] if sample == "sample_a" else [records]
        for i, part in enumerate(parts):
            write_gz(os.path.join(OUT, "reads", sample, f"part_{i}.fastq.gz"), "".join(part))

    with open(os.path.join(OUT, "samplesheet.csv"), "w") as fh:
        fh.write("sample,fastq\n")
        for sample in SAMPLES:
            fh.write(f"{sample},tests/data/reads/{sample}/*.fastq.gz\n")

    # NCBI dump format: fields joined by "\t|\t", each line ending "\t|"
    tax_dir = os.path.join(OUT, "kraken2_taxonomy")
    os.makedirs(tax_dir, exist_ok=True)
    with open(os.path.join(tax_dir, "nodes.dmp"), "w") as fh:
        for taxid, parent, rank, _ in TAXONOMY:
            fh.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\n")
    with open(os.path.join(tax_dir, "names.dmp"), "w") as fh:
        for taxid, _, _, name in TAXONOMY:
            fh.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    with open(os.path.join(tax_dir, "seqid2taxid.map"), "w") as fh:
        for genome, taxid in GENOME_TAXIDS.items():
            fh.write(f"{genome}\t{taxid}\n")


if __name__ == "__main__":
    main()
