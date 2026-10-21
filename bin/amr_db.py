#!/usr/bin/env python3
"""
Lay out one of ABRicate's bundled AMR databases for AMR_MINIMAP2 and
AMR_REPORT: its sequences as a FASTA (renamed 1, 2, 3...), and a table
giving each one's name in the database (`gene`, `accession`: how
ABRicate's hits are matched back), allele name for reports, gene family
and drug class.

Families and classes come from NCBI's Reference Gene Catalog
(assets/amr/ncbi_reference_genes.tsv), so results read the same whichever
database found them. A gene from another database gets the NCBI
annotation of the NCBI allele with the identical sequence (either strand),
found in ABRicate's own `ncbi` database. Most genes have one (73% of
megares to 91% of card); the rest are annotated from their own database's
header, as far as it allows, and left empty otherwise:

  resfinder  >resfinder~~~blaCTX-M-15_1~~~AY044436~~~Amoxicillin;Ampicillin;...
             family = name without its _<n> variant number; antibiotics -> drug_subclass
  argannot   >argannot~~~(Bla)blaTEM-1D~~~...~~~ (Bla)blaTEM-1D
             family = name without its (class) prefix; that class -> drug_class
  card       >card~~~TEM-1~~~...~~~monobactam;cephalosporin;... description
             family = name; drug classes -> drug_class
  megares    >megares~~~TEM~~~MEG_7181~~~ Drugs:betalactams:Class_A_betalactamases:TEM
             family = name (MEGARes's group); class -> drug_class

The allele is the database's own name, minus resfinder's _<n> and
argannot's (class). MEGARes names only the group, so its allele is the
identical NCBI allele's name where there is one, else `<group>
(<MEG id>)`.

Drug classes are written as NCBI writes them: upper case, '/'-joined.
Pure standard library: it runs in the abricate image.

  amr_db.py --datadir /usr/local/db --db resfinder \\
      --catalog ncbi_reference_genes.tsv --fasta amr_genes.fasta --table amr_genes.tsv
"""
import argparse
import csv
import os
import re

COLS = ["id", "gene", "accession", "allele", "length", "gene_family", "drug_class", "drug_subclass"]

# ARG-ANNOT's class codes, as NCBI names those classes; codes not listed
# (Mca, Tcm) leave drug_class empty
ARGANNOT_CLASS = {
    "AGly": "AMINOGLYCOSIDE", "AGly_Flqn": "AMINOGLYCOSIDE/QUINOLONE", "Bla": "BETA-LACTAM",
    "Col": "COLISTIN", "Fcd": "FUSIDIC ACID", "Fcyn": "FOSFOMYCIN", "Flq": "QUINOLONE",
    "Gly": "GLYCOPEPTIDE", "MLS": "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN", "Ntmdz": "NITROIMIDAZOLE",
    "Oxzln": "OXAZOLIDINONE", "Phe": "PHENICOL", "Rif": "RIFAMYCIN", "Sul": "SULFONAMIDE",
    "Tet": "TETRACYCLINE", "Tmt": "TRIMETHOPRIM",
}


def read_fasta(path):
    """[(header, sequence)], sequences upper case."""
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append([line[1:], []])
            elif line:
                records[-1][1].append(line.strip())
    return [(h, "".join(s).upper()) for h, s in records]


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def ncbi_style(names):
    """'cephalosporin;penicillin_beta-lactam' -> 'CEPHALOSPORIN/PENICILLIN BETA-LACTAM'"""
    return "/".join(n.strip().replace("_", " ").upper() for n in names if n.strip())


def allele_name(db, gene, accession, ncbi):
    if db == "resfinder":
        return re.sub(r"_\d+$", "", gene)
    if db == "argannot":
        return re.sub(r"^\(\w+\)", "", gene)
    if db == "megares":
        return ncbi["gene"] if ncbi else f"{gene} ({accession})"
    return gene


def own_annotation(db, gene, description):
    """(gene_family, drug_class, drug_subclass) from the database's own header."""
    if db == "resfinder":
        # description: "<antibiotic>;<antibiotic>... <gene>", or " <gene>" with none
        drugs = description.split(" ")[0] if not description.startswith(" ") else ""
        return re.sub(r"_\d+$", "", gene), "", ncbi_style(drugs.split(";"))
    if db == "argannot":
        m = re.match(r"\((\w+)\)(.*)", gene)
        return (m.group(2), ARGANNOT_CLASS.get(m.group(1), ""), "") if m else (gene, "", "")
    if db == "card":
        return gene, ncbi_style(description.split(" ")[0].split(";")), ""
    if db == "megares":
        # Drugs:<class>:<mechanism>:<group>, or Multi-compound:<class>:...
        parts = description.strip().split(":")
        return gene, ncbi_style(parts[1:2]), ""
    return gene, "", ""


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datadir", required=True, help="ABRicate's database folder (abricate --help shows it)")
    p.add_argument("--db", required=True, choices=["ncbi", "resfinder", "argannot", "card", "megares"])
    p.add_argument("--catalog", required=True, help="assets/amr/ncbi_reference_genes.tsv")
    p.add_argument("--fasta", required=True, help="output FASTA, sequences named by id")
    p.add_argument("--table", required=True, help="output annotation table")
    args = p.parse_args()

    with open(args.catalog) as fh:
        catalog = {r["accession"]: r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}

    # NCBI annotation of each NCBI sequence, either strand
    by_seq = {}
    for header, seq in read_fasta(os.path.join(args.datadir, "ncbi", "sequences")):
        _, gene, accession, _ = header.split("~~~", 3)
        if accession not in catalog:
            raise SystemExit(f"{accession} ({gene}) isn't in {args.catalog}: was it made from the catalog this ABRicate database was built from?")
        by_seq[seq] = by_seq[revcomp(seq)] = catalog[accession]

    rows = []
    with open(args.fasta, "w") as fa:
        for i, (header, seq) in enumerate(read_fasta(os.path.join(args.datadir, args.db, "sequences")), 1):
            # >db~~~gene~~~accession~~~description
            _, gene, accession, description = header.split("~~~", 3)
            ncbi = by_seq.get(seq)
            if ncbi:
                family, drug_class, drug_subclass = ncbi["gene_family"], ncbi["drug_class"], ncbi["drug_subclass"]
            else:
                family, drug_class, drug_subclass = own_annotation(args.db, gene, description)
            rows.append({
                "id": i, "gene": gene, "accession": accession,
                "allele": allele_name(args.db, gene, accession, ncbi), "length": len(seq),
                "gene_family": family, "drug_class": drug_class, "drug_subclass": drug_subclass,
            })
            fa.write(f">{i}\n{seq}\n")

    with open(args.table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
