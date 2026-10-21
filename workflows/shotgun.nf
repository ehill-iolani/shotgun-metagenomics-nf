include { MERGE_FASTQ    } from '../modules/merge_fastq.nf'
include { CHOPPER        } from '../modules/chopper.nf'
include { DEACON_DEPLETE; HOST_REMOVAL_REPORT } from '../modules/deacon_deplete.nf'
include { READ_STATS; READ_STATS_REPORT       } from '../modules/read_stats.nf'

workflow WF_SHOTGUN {

    take:
    reads_ch     // tuple(sample, fastq)
    host_index   // Deacon index of the host genome(s) to remove

    main:
    // 1. merge multi-part fastq(.gz) files per barcode into one file per sample
    MERGE_FASTQ(reads_ch)

    // 2. length/quality filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 3. human (host) read removal. Runs before anything is published, so
    // no published file holds a read from the person the sample came from.
    // .first() makes the single-emission index a value channel, reused for
    // every sample instead of being consumed by the first one
    DEACON_DEPLETE(CHOPPER.out.filtered, host_index.first())
    HOST_REMOVAL_REPORT(DEACON_DEPLETE.out.summary.map { sample, json -> json }.collect())

    // 3b. read length / Q-score summary at each stage -- a side branch,
    // nothing downstream depends on it. Publishes final_report/read_stats.tsv,
    // read_length_qscore_hist.tsv and read_qc_summary.html, which the
    // platform's Read QC tab charts.
    READ_STATS(
        MERGE_FASTQ.out.merged
            .join(CHOPPER.out.filtered)
            .join(DEACON_DEPLETE.out.reads)
    )
    READ_STATS_REPORT(
        READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
        READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
    )

    emit:
    reads = DEACON_DEPLETE.out.reads
}
