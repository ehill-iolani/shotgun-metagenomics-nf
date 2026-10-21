process CHOPPER {
    tag "$sample"
    label 'process_low'
    // same image as 16S-nf, so it's already in the Artifact Registry mirror
    container "${params.container_registry}/biocontainers/chopper:0.7.0--hdcf5f25_0"
    // not published, same reason as MERGE_FASTQ

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: filtered

    script:
    """
    set -o pipefail

    zcat -f ${fastq} \\
      | chopper -t ${task.cpus} -q ${params.min_qual} -l ${params.min_len} \\
      | gzip > ${sample}.filtered.fastq.gz
    """
}
