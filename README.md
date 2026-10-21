# shotgun-metagenomics-nf

Taxonomic identification from unamplified Oxford Nanopore reads: samples
prepared with the Rapid Barcoding Kit (SQK-RBK114.24/.96), R10.4.1, one
MinION flow cell per run. Built for the iolani-bioinformatics platform,
following the stages of EPI2ME's
[wf-metagenomics](https://github.com/epi2me-labs/wf-metagenomics), but
written from scratch.

**Status:** read QC, human read removal and Kraken2 classification are
in. AMR gene detection comes next; see [docs/CONTRACT.md](docs/CONTRACT.md)
for the planned outputs.

## Steps

1. **Merge** each sample's part-files (a MinKNOW barcode folder) into one.
2. **Filter** with [chopper](https://github.com/wdecoster/chopper): at least
   500 bp and mean Q10 by default.
3. **Remove human reads** with [Deacon](https://github.com/bede/deacon) and
   its panhuman-1 index (47 human haplotypes + CHM13 + GRCh38, downloaded
   once and cached). Reads are only ever published after this step.
4. **Read QC** at each stage: length and Q-score summaries.
5. **Classify** with [Kraken2](https://github.com/DerrickWood/kraken2)
   against the PlusPF-16 database (RefSeq bacteria, archaea, viruses,
   plasmids, human, protozoa and fungi), with `--confidence 0.05
   --minimum-hit-groups 3`. All samples go through one Kraken2 run, so the
   database loads once. No Bracken: its read-length model doesn't fit long
   reads. Taxa under `--min_reads` (10) or `--min_rel_abundance` (0.0001)
   are flagged `below_threshold`, not dropped.

## Usage

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

`samplesheet.csv` has two columns, `sample,fastq`. `fastq` may be a glob,
e.g. `fastq_pass/barcode01/*.fastq.gz`, which resolves relative to the
launch directory.

The first run downloads the 3.3 GB human index and the 11 GB PlusPF-16
archive (15 GB unpacked) into `assets/db_cache/` (`--db_cache`); later runs
reuse them. Kraken2 needs about 20 GB of memory for PlusPF-16. For a
quicker test on real data, or a smaller machine, `--kraken2_db pluspf8`
uses PlusPF-8 instead (5.5 GB download, about 10 GB of memory), which misses
more reads.

### Test

```bash
nextflow run main.nf -profile test,docker
```

Runs offline on the synthetic data in `tests/data/` (regenerate with
`python3 tests/make_test_data.py`). It builds its host index from a small
synthetic "host" genome instead of downloading the human one, and its
Kraken2 database from the synthetic microbes and a made-up taxonomy
(`tests/data/kraken2_taxonomy/`) instead of downloading PlusPF.

## Outputs

See [docs/CONTRACT.md](docs/CONTRACT.md).
