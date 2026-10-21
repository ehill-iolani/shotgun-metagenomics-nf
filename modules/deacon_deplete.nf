// Human (host) read removal: drop every read sharing enough minimizers with
// the host index (Deacon's defaults, >= 2 and >= 1% of the read's
// minimizers). Deacon had the best balanced accuracy for long-read panhuman
// depletion in its benchmark, removing ~20x fewer bacterial reads than
// Hostile (https://doi.org/10.1101/2025.06.09.658732).
process DEACON_DEPLETE {
    tag "$sample"
    // the panhuman-1 index is 3.3 GB on disk, ~5 GB in memory
    label 'process_medium'
    container "${params.container_registry}/biocontainers/deacon:0.18.0--h3edb6b3_0"
    publishDir(path: { "${params.outdir}/${sample}/02_host_removed" }, mode: 'copy', pattern: '*.fastq.gz')

    input:
    tuple val(sample), path(fastq)
    path host_index

    output:
    tuple val(sample), path("${sample}.host_removed.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.deacon.json"),           emit: summary

    script:
    """
    deacon filter -d -t ${task.cpus} \\
        ${host_index} ${fastq} \\
        -o ${sample}.host_removed.fastq.gz \\
        -s ${sample}.deacon.json
    """
}

process HOST_REMOVAL_REPORT {
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path summaries

    output:
    path "host_removal.tsv", emit: table

    script:
    """
    host_removal_report.py --out host_removal.tsv ${summaries}
    """
}
