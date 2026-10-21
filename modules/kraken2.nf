// Taxonomic classification of every sample's host-removed reads in one
// Kraken2 run, so its database (15 GB for PlusPF-16) is read into memory
// once per pipeline run rather than once per sample. Each read ID gets its
// sample's index as a prefix (`<i>|<read id>`), which TAXON_ABUNDANCE uses
// to split the calls back out by sample.
//
// --confidence 0.05: a read is only called at a taxon when >= 5% of its
// k-mers fall in that taxon or below it, else the call moves up the tree
// until they do. Kraken2's default (0) takes any call at all, so a few
// stray k-mers in a multi-kb read can make a species call.
// --minimum-hit-groups 3 (default 2) asks for one more independent hit
// before any call.
process KRAKEN2 {
    // memory follows the database instead, set in nextflow.config
    label 'process_high'
    container "${params.container_registry}/biocontainers/kraken2:2.17.2--pl5321h3be2455_0"

    input:
    tuple val(samples), path(fastqs)
    path db

    output:
    tuple val(samples), path("kraken2.reads.tsv.gz"), emit: reads

    script:
    def files = fastqs instanceof List ? fastqs : [fastqs]
    // only the first four columns of Kraken2's per-read output (C/U, read
    // ID, taxid, length): the fifth, its k-mer-by-k-mer LCA list, runs to
    // hundreds of entries per long read and nothing reads it
    """
    set -o pipefail

    i=0
    for fq in ${files.join(' ')}; do
        zcat -f "\$fq" | awk -v i="\$i" 'NR % 4 == 1 { print "@" i "|" substr(\$0, 2); next } { print }'
        i=\$((i + 1))
    done \\
      | kraken2 --db ${db} --threads ${task.cpus} \\
            --confidence 0.05 --minimum-hit-groups 3 \\
            /dev/stdin \\
      | cut -f 1-4 \\
      | gzip > kraken2.reads.tsv.gz
    """
}

process TAXON_ABUNDANCE {
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    // per-sample reports into <sample>/03_kraken2/, the table into final_report/
    publishDir(
        path: "${params.outdir}",
        mode: 'copy',
        saveAs: { fn ->
            fn.endsWith('.kraken2.report.txt')
                ? "${fn - '.kraken2.report.txt'}/03_kraken2/${fn}"
                : "final_report/${fn}"
        }
    )

    input:
    tuple val(samples), path(reads)
    path db

    output:
    path "taxon_abundance.tsv",   emit: table
    path "*.kraken2.report.txt",  emit: reports

    script:
    """
    taxon_abundance.py \\
        --reads ${reads} \\
        --samples ${samples.join(' ')} \\
        --nodes ${db}/nodes.dmp \\
        --names ${db}/names.dmp \\
        --min-reads ${params.min_reads} \\
        --min-rel-abundance ${params.min_rel_abundance} \\
        --out taxon_abundance.tsv \\
        --report-dir .
    """
}
