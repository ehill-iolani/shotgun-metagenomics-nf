#!/usr/bin/env python3
"""
Writes assets/amr/ncbi_reference_genes.tsv, the gene family and drug class
of every gene in ABRicate's `ncbi` database, which AMR_REPORT names its
results by.

ABRicate 1.4.0's `ncbi` database (in the pinned abricate image, built
2026-04-03) is the core AMR genes of NCBI's Reference Gene Catalog release
2026-03-24.1: every one of its 8,232 sequences is a core/AMR/AMR row of that
release, matched on RefSeq accession. Its FASTA headers carry only the
allele name and subclass, so this pulls the family and class from the
catalog of that same release (pinned by md5 below).

One change from the catalog: NCBI files all 1,381 OXA alleles under one
family, `blaOXA`, which would put an acquired carbapenemase (OXA-48) in the
same row as a gene every Acinetobacter baumannii carries (OXA-51). They're
split here by the OXA-n family their product names give ("OXA-48 family
carbapenem-hydrolyzing class D beta-lactamase OXA-181" -> blaOXA-48-like).
The namesake allele joins its family (blaOXA-48 -> blaOXA-48-like), and so
do lettered variants of one number (blaOXA-114a to -114x ->
blaOXA-114-like); any other OXA allele stays a family of its own.

  python3 scripts/make_amr_catalog.py

Output is deterministic, so rerunning without changes leaves git with
nothing to commit.
"""
import csv
import hashlib
import io
import os
import re
import urllib.request

RELEASE = "2026-03-24.1"
URL = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/{RELEASE}/ReferenceGeneCatalog.txt"
MD5 = "57e1d18b5ca5e218c97194fe9b981c03"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "amr", "ncbi_reference_genes.tsv")
COLS = ["accession", "gene", "gene_family", "drug_class", "drug_subclass"]

OXA_FAMILY = re.compile(r"\bOXA-(\S+) family\b")


def main():
    raw = urllib.request.urlopen(URL).read()
    if hashlib.md5(raw).hexdigest() != MD5:
        raise SystemExit(f"md5 mismatch for {URL}")
    # ABRicate's `ncbi` keeps the same subset (abricate-get_db, get_ncbi)
    rows = [r for r in csv.DictReader(io.StringIO(raw.decode()), delimiter="\t")
            if (r["scope"], r["type"], r["subtype"]) == ("core", "AMR", "AMR")]

    oxa_families = set()
    for r in rows:
        m = OXA_FAMILY.search(r["product_name"])
        if r["gene_family"] == "blaOXA" and m:
            oxa_families.add(m.group(1))

    out = []
    for r in rows:
        # ABRicate names a sequence by its allele, or its hierarchy node
        # when the catalog gives no allele (e.g. BcII, whose family is bla2)
        gene = r["allele"] or r["hierarchy_node"]
        family = r["gene_family"]
        if family == "blaOXA":
            m = OXA_FAMILY.search(r["product_name"])
            number = gene.removeprefix("blaOXA-")
            base = re.sub(r"(?<=\d)[a-z]$", "", number)
            if m:
                family = f"blaOXA-{m.group(1)}-like"
            elif base in oxa_families or base != number:
                family = f"blaOXA-{base}-like"
            else:
                family = gene
        out.append({
            "accession": r["refseq_nucleotide_accession"], "gene": gene, "gene_family": family,
            "drug_class": r["class"], "drug_subclass": r["subclass"],
        })
    out.sort(key=lambda r: r["accession"])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        fh.write(f"# NCBI Reference Gene Catalog {RELEASE}, core AMR genes; written by scripts/make_amr_catalog.py\n")
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} genes, {len({r['gene_family'] for r in out})} families -> {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()
