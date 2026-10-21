// Host-removal index from a FASTA of host genome(s), for --host_fasta --
// the test profile's synthetic host, or a non-human host. Default k/w
// (31/15), the same as the prebuilt panhuman-1 index.
process BUILD_DEACON_INDEX {
    label 'process_medium'
    container "${params.container_registry}/biocontainers/deacon:0.18.0--h3edb6b3_0"

    input:
    path host_fasta

    output:
    path "host.idx", emit: index

    script:
    """
    deacon index build -t ${task.cpus} ${host_fasta} > host.idx
    """
}
