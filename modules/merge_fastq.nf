process MERGE_FASTQ {
    tag "$sample"
    label 'process_low'
    // not published: these reads still include any human reads (see
    // docs/CONTRACT.md), so they stay in the work dir only

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: merged
    tuple val(sample), env(model),                        emit: model

    script:
    // ONT barcode folders are usually split into many part-files
    // (fastq_pass/barcodeNN/*.fastq.gz); cat them into one file per sample
    // so downstream steps don't need to care how many parts came in.
    def infiles = fastq instanceof List ? fastq : [fastq]
    // The basecall model named in the first read's header, for WF_SHOTGUN's
    // fast-model warning: MinKNOW writes basecall_model_version_id=<model>,
    // Dorado an RG tag ending _<model>. Empty when there's neither
    """
    cat ${infiles.join(' ')} > ${sample}.merged.fastq.gz
    model=\$(zcat -f ${sample}.merged.fastq.gz | head -n 1 | grep -oE '[A-Za-z0-9_.]+_(fast|hac|sup)@v[0-9.]+' | head -n 1 || true)
    """
}
