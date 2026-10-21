#!/usr/bin/env nextflow
/*
 * shotgun-metagenomics-nf: taxonomic identification and AMR gene
 * detection from unamplified ONT reads -- Rapid Barcoding Kit
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
include { FETCH_KRAKEN2_DB; kraken2Dbs } from './modules/fetch_kraken2_db.nf'
include { BUILD_KRAKEN2_DB   } from './modules/build_kraken2_db.nf'

// ---- top-level params (override via -params-file or --flag) ----
// filtering/host-removal/classification defaults live in nextflow.config, not here --
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
      --kraken2_db  Kraken2 database, downloaded once and cached (default ${params.kraken2_db}):
                      pluspf16   PlusPF-16: RefSeq bacteria, archaea, viruses, plasmids,
                                 human, protozoa and fungi, capped at 16 GB (needs ~20 GB RAM)
                      pluspf8    PlusPF-8: the same, capped at 8 GB (needs ~10 GB RAM) --
                                 less sensitive, for testing or smaller machines
      --min_reads / --min_rel_abundance
                    taxa with fewer reads, or a smaller fraction of the sample's
                    reads, get below_threshold=true in taxon_abundance.tsv, not
                    dropped (default ${params.min_reads} / ${params.min_rel_abundance}, 0 = off)
      --amr_tool    how reads are matched to AMR genes (default ${params.amr_tool}):
                      minimap2   fast; the default
                      abricate   ABRicate (BLAST): the same calls in testing, ~25x slower
      --amr_db      AMR gene database, from those bundled with ABRicate 1.4.0 (default ${params.amr_db}):
                      ncbi       NCBI Reference Gene Catalog (AMRFinderPlus), acquired + intrinsic AMR genes
                      resfinder  ResFinder, acquired genes
                      argannot   ARG-ANNOT
                      card       CARD -- advanced: also lists housekeeping efflux pumps and their
                                 regulators, so its results aren't comparable with the others
                      megares    MEGARes -- advanced, for the same reason as card
      --amr_min_identity / --amr_min_coverage
                    a read's hit to a gene counts at this percent identity, over at
                    least this percent of the gene's length (default ${params.amr_min_identity} / ${params.amr_min_coverage})
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
    if (!(params.kraken2_db in kraken2Dbs())) {
        problems << "--kraken2_db must be one of ${kraken2Dbs().keySet().join(', ')} (got '${params.kraken2_db}')"
    }
    if (!(params.min_reads instanceof Integer) || params.min_reads < 0) {
        problems << "--min_reads must be a whole number >= 0 (got '${params.min_reads}')"
    }
    if (!(params.min_rel_abundance instanceof Number) || params.min_rel_abundance < 0 || params.min_rel_abundance > 1) {
        problems << "--min_rel_abundance must be a fraction from 0 to 1 (got '${params.min_rel_abundance}')"
    }
    if (!(params.amr_tool in ['minimap2', 'abricate'])) {
        problems << "--amr_tool must be minimap2 or abricate (got '${params.amr_tool}')"
    }
    if (!(params.amr_db in ['ncbi', 'resfinder', 'argannot', 'card', 'megares'])) {
        problems << "--amr_db must be one of ncbi, resfinder, argannot, card, megares (got '${params.amr_db}')"
    }
    if (!(params.amr_min_identity instanceof Number) || params.amr_min_identity < 0 || params.amr_min_identity > 100) {
        problems << "--amr_min_identity must be a percentage from 0 to 100 (got '${params.amr_min_identity}')"
    }
    if (!(params.amr_min_coverage instanceof Number) || params.amr_min_coverage < 0 || params.amr_min_coverage > 100) {
        problems << "--amr_min_coverage must be a percentage from 0 to 100 (got '${params.amr_min_coverage}')"
    }
    if ((params.kraken2_build_fasta as boolean) != (params.kraken2_build_taxonomy as boolean)) {
        problems << "--kraken2_build_fasta and --kraken2_build_taxonomy go together: set both or neither"
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

    // same idea for the Kraken2 database: the test profile builds a tiny
    // one from its synthetic genomes instead of downloading PlusPF
    if (params.kraken2_build_fasta) {
        BUILD_KRAKEN2_DB(
            Channel.fromPath(params.kraken2_build_fasta, checkIfExists: true).collect(),
            Channel.fromPath(params.kraken2_build_taxonomy, type: 'dir', checkIfExists: true)
        )
        kraken2_db_ch = BUILD_KRAKEN2_DB.out.db
    } else {
        FETCH_KRAKEN2_DB(params.kraken2_db)
        kraken2_db_ch = FETCH_KRAKEN2_DB.out.db
    }

    WF_SHOTGUN(samplesheetToChannel(params.input), host_index_ch, kraken2_db_ch)
}
