// AMR gene detection in host-removed reads, against one of the AMR
// databases bundled in the pinned ABRicate image (--amr_db, all dated
// 2026-04-03; `ncbi` is NCBI's Reference Gene Catalog 2026-03-24.1, core
// AMR genes). Reads are matched with minimap2 (--amr_tool minimap2, the
// default) or ABRicate itself (--amr_tool abricate). Both found the same
// genes on simulated R10 reads; ABRicate took ~25x longer.
//
// Every read is searched, not an assembly: with one MinION flow cell a
// sample's minor members are rarely deep enough to assemble, and a long
// read usually spans a whole gene (~1 kb) plus enough flank for Kraken2 to
// say whose it is.

// The chosen database as a FASTA for minimap2, and a table of each
// sequence's gene family and drug class (bin/amr_db.py). Taken from the
// abricate image itself, so both tools search exactly the same sequences
process AMR_DB {
    label 'process_low'
    container "${params.container_registry}/biocontainers/abricate:1.4.0--h05cac1d_0"

    input:
    val db          // ncbi, resfinder, argannot, card or megares
    path catalog    // assets/amr/ncbi_reference_genes.tsv

    output:
    path "amr_db.fasta", emit: fasta
    path "amr_db.tsv",   emit: table

    script:
    // ABRicate's database folder differs between its container and a
    // conda install; its --help gives the one in use
    """
    datadir=\$(abricate --help | sed -n 's/.*--datadir.*\\[\\(.*\\)\\].*/\\1/p')
    amr_db.py --datadir "\$datadir" --db ${db} --catalog ${catalog} \\
        --fasta amr_db.fasta --table amr_db.tsv
    """
}

// -c: base-level alignment, so the identity (matching bases / alignment
// columns) is exact rather than minimap2's chaining estimate.
// --secondary=no: only each stretch of read's best allele; near-identical
// alleles would otherwise each add a line per read
process AMR_MINIMAP2 {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/minimap2:2.31--h118bc1c_0"

    input:
    tuple val(sample), path(fastq)
    path genes

    output:
    tuple val(sample), path("${sample}.amr_hits.paf"), emit: hits

    script:
    """
    minimap2 -c -x map-ont --secondary=no -t ${task.cpus} ${genes} ${fastq} > ${sample}.amr_hits.paf
    """
}

// ABRicate takes FASTA, and its any2fasta step holds a whole input file in
// memory (~4-9 bytes per base, so ~40-90 GB for a flow cell), so reads go
// through it 50,000 at a time: ~100-250 Mb, 1-2 GB, for rapid-kit reads.
// Its BLAST threading gains little past 4, hence process_medium
process AMR_ABRICATE {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/abricate:1.4.0--h05cac1d_0"

    input:
    tuple val(sample), path(fastq)
    val db

    output:
    tuple val(sample), path("${sample}.amr_hits.tsv"), emit: hits

    script:
    def abricate = "abricate --db ${db} --minid ${params.amr_min_identity} --mincov ${params.amr_min_coverage} --threads ${task.cpus} --noheader --nopath"
    """
    set -o pipefail
    : > ${sample}.amr_hits.tsv

    zcat -f ${fastq} | awk -v n=50000 -v cmd="${abricate}" -v out=${sample}.amr_hits.tsv '
        function flush() {
            close("chunk.fa")
            if (system(cmd " chunk.fa >> " out) != 0) { failed = 1; exit 1 }
            reads = 0
        }
        NR % 4 == 1 { if (reads == n) flush(); print ">" substr(\$1, 2) > "chunk.fa"; reads++; next }
        NR % 4 == 2 { print > "chunk.fa" }
        END { if (failed) exit 1; if (reads) flush() }'
    rm -f chunk.fa
    """
}

// final_report/amr_genes.tsv and <sample>/04_amr/<sample>.amr_reads.tsv,
// with each gene's host genus from KRAKEN2's calls for the same reads
process AMR_REPORT {
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
        saveAs: { fn ->
            fn.endsWith('.amr_reads.tsv')
                ? "${fn - '.amr_reads.tsv'}/04_amr/${fn}"
                : "final_report/${fn}"
        }
    )

    input:
    tuple val(samples), path(kraken2_reads)
    path hits        // <sample>.amr_hits.{paf,tsv}, one per sample
    path genes       // AMR_DB's table
    path kraken2_db
    path intrinsic   // assets/amr/intrinsic_genes.tsv

    output:
    path "amr_genes.tsv",    emit: table
    path "*.amr_reads.tsv",  emit: reads

    script:
    """
    amr_report.py \\
        --hits ${hits} \\
        --format ${params.amr_tool == 'abricate' ? 'abricate' : 'paf'} \\
        --genes ${genes} \\
        --samples ${samples.join(' ')} \\
        --kraken2 ${kraken2_reads} \\
        --nodes ${kraken2_db}/nodes.dmp \\
        --names ${kraken2_db}/names.dmp \\
        --intrinsic ${intrinsic} \\
        --min-identity ${params.amr_min_identity} \\
        --min-coverage ${params.amr_min_coverage} \\
        --out amr_genes.tsv \\
        --reads-dir .
    """
}
