include { MERGE_FASTQ    } from '../modules/merge_fastq.nf'
include { CHOPPER        } from '../modules/chopper.nf'
include { DEACON_DEPLETE; HOST_REMOVAL_REPORT } from '../modules/deacon_deplete.nf'
include { READ_STATS; READ_STATS_REPORT       } from '../modules/read_stats.nf'
include { KRAKEN2; TAXON_ABUNDANCE            } from '../modules/kraken2.nf'
include { AMR_DB; AMR_MINIMAP2; AMR_ABRICATE; AMR_REPORT } from '../modules/amr.nf'

workflow WF_SHOTGUN {

    take:
    reads_ch     // tuple(sample, fastq)
    host_index   // Deacon index of the host genome(s) to remove
    kraken2_db   // Kraken2 database dir (.k2d files + nodes.dmp/names.dmp)

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

    // 4. taxonomic classification: every sample in one Kraken2 task, so the
    // database loads once. Sorted by sample so the task's inputs, and so
    // -resume's cache key, don't depend on which sample finished first
    KRAKEN2(
        DEACON_DEPLETE.out.reads
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { it.transpose() },   // [[samples], [fastqs]]
        kraken2_db
    )
    // publishes final_report/taxon_abundance.tsv and
    // <sample>/03_kraken2/<sample>.kraken2.report.txt
    TAXON_ABUNDANCE(KRAKEN2.out.reads, kraken2_db)

    // 5. AMR genes in each sample's host-removed reads, with minimap2 or
    // ABRicate, then one report across samples that names each gene's
    // host genus from step 4's calls for the same reads
    AMR_DB(params.amr_db, file("${projectDir}/assets/amr/ncbi_reference_genes.tsv", checkIfExists: true))
    if (params.amr_tool == 'abricate') {
        AMR_ABRICATE(DEACON_DEPLETE.out.reads, params.amr_db)
        amr_hits = AMR_ABRICATE.out.hits
    } else {
        AMR_MINIMAP2(DEACON_DEPLETE.out.reads, AMR_DB.out.fasta)
        amr_hits = AMR_MINIMAP2.out.hits
    }
    // publishes final_report/amr_genes.tsv and
    // <sample>/04_amr/<sample>.amr_reads.tsv. Sorted for the same reason
    // as KRAKEN2's input
    AMR_REPORT(
        KRAKEN2.out.reads,
        amr_hits.toSortedList { a, b -> a[0] <=> b[0] }.map { it.collect { sample, hits -> hits } },
        AMR_DB.out.table,
        kraken2_db,
        file(params.amr_intrinsic, checkIfExists: true)
    )

    emit:
    reads   = DEACON_DEPLETE.out.reads
    kraken2 = KRAKEN2.out.reads   // tuple([samples], per-read calls); read IDs prefixed with the sample's index
}
