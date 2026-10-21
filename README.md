# shotgun-metagenomics-nf

Taxonomic identification from unamplified Oxford Nanopore reads: samples
prepared with the Rapid Barcoding Kit (SQK-RBK114.24/.96), R10.4.1, one
MinION flow cell per run. Built for the iolani-bioinformatics platform,
following the stages of EPI2ME's
[wf-metagenomics](https://github.com/epi2me-labs/wf-metagenomics), but
written from scratch.

**Status:** read QC and human read removal are in. Kraken2 classification
and AMR gene detection come next; see [docs/CONTRACT.md](docs/CONTRACT.md)
for the planned outputs.

## Steps

1. **Merge** each sample's part-files (a MinKNOW barcode folder) into one.
2. **Filter** with [chopper](https://github.com/wdecoster/chopper): at least
   500 bp and mean Q10 by default.
3. **Remove human reads** with [Deacon](https://github.com/bede/deacon) and
   its panhuman-1 index (47 human haplotypes + CHM13 + GRCh38, downloaded
   once and cached). Reads are only ever published after this step.
4. **Read QC** at each stage: length and Q-score summaries.

## Usage

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

`samplesheet.csv` has two columns, `sample,fastq`. `fastq` may be a glob,
e.g. `fastq_pass/barcode01/*.fastq.gz`, which resolves relative to the
launch directory.

The first run downloads the 3.3 GB human index into `assets/db_cache/`
(`--db_cache`); later runs reuse it.

### Test

```bash
nextflow run main.nf -profile test,docker
```

Runs offline on the synthetic data in `tests/data/` (regenerate with
`python3 tests/make_test_data.py`). It builds its host index from a small
synthetic "host" genome instead of downloading the human one.

## Outputs

See [docs/CONTRACT.md](docs/CONTRACT.md).
