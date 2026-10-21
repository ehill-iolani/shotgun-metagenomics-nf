#!/usr/bin/env nextflow
/*
 * shotgun-metagenomics-nf: taxonomic identification (and, later, AMR gene
 * detection) from unamplified ONT reads -- Rapid Barcoding Kit
 * (SQK-RBK114), R10.4.1, one MinION flow cell per run.
 *
 * Follows EPI2ME wf-metagenomics' stages (ingest -> filter -> host removal
 * -> classify -> abundance report) but is written from scratch: that
 * workflow's code is under ONT's research-only licence. Its defaults are
 * also loose for shotgun data (no quality filter, Kraken2 confidence 0,
 * Bracken modelled on 300 bp reads), so this one filters by quality,
 * always removes human reads, and reports Kraken2 calls directly.
 *
 * Output files and their columns: docs/CONTRACT.md.
 */

nextflow.enable.dsl = 2

include { WF_SHOTGUN         } from './workflows/shotgun.nf'
include { FETCH_DEACON_INDEX } from './modules/fetch_deacon_index.nf'
include { BUILD_DEACON_INDEX } from './modules/build_deacon_index.nf'

// ---- top-level params (override via -params-file or --flag) ----
// filtering/host-removal defaults live in nextflow.config, not here --
// anything read directly inside an included module or workflow script
// (rather than only inside this file's own `workflow` block) must be set
// there to be reliably visible by the time that module/workflow script runs
params.input  = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir = "results"
params.help   = false

def helpMessage() {
    log.info """
    shotgun-metagenomics-nf: ONT shotgun metagenomics (Rapid Barcoding Kit)
    ------------------------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Required:
      --input       CSV: sample,fastq_path

    Key optional:
      --min_len / --min_qual   chopper filtering thresholds (default ${params.min_len} bp / Q${params.min_qual})
      --host_fasta  build the host-removal index from this FASTA instead of
                    downloading the panhuman-1 human pangenome index
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

// Nothing validates params against nextflow_schema.json (nf-schema needs a
// newer Nextflow than the platform runs), so check the ones that would
// otherwise fail deep inside a task, or not fail at all
def checkParams() {
    def problems = []
    if (!(params.min_len instanceof Integer) || params.min_len < 1) {
        problems << "--min_len must be a positive whole number (got '${params.min_len}')"
    }
    if (!(params.min_qual instanceof Number) || params.min_qual < 0) {
        problems << "--min_qual must be a number >= 0 (got '${params.min_qual}')"
    }
    if (problems) {
        problems.each { log.error it }
        exit 1
    }
}

workflow {
    if (params.help || !params.input) {
        helpMessage()
        exit 0
    }
    checkParams()

    if (params.host_fasta) {
        BUILD_DEACON_INDEX(Channel.fromPath(params.host_fasta, checkIfExists: true))
        host_index_ch = BUILD_DEACON_INDEX.out.index
    } else {
        FETCH_DEACON_INDEX()
        host_index_ch = FETCH_DEACON_INDEX.out.index
    }

    WF_SHOTGUN(samplesheetToChannel(params.input), host_index_ch)
}
