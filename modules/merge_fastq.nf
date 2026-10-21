process MERGE_FASTQ {
    tag "$sample"
    label 'process_low'
    // not published: these reads still include any human reads (see
    // docs/CONTRACT.md), so they stay in the work dir only

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: merged

    script:
    // ONT barcode folders are usually split into many part-files
    // (fastq_pass/barcodeNN/*.fastq.gz); cat them into one file per sample
    // so downstream steps don't need to care how many parts came in.
    def infiles = fastq instanceof List ? fastq : [fastq]
    """
    cat ${infiles.join(' ')} > ${sample}.merged.fastq.gz
    """
}
