process FETCH_DEACON_INDEX {
    label 'process_low'
    // needs /bin/bash to exist (the Nextflow docker executor execs it
    // directly, regardless of the `shell` directive) -- ubuntu has it plus
    // md5sum already installed, only curl needs adding
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    storeDir "${params.db_cache}/deacon"

    output:
    path "panhuman-1.k31w15.idx", emit: index

    script:
    // panhuman-1: minimizers of 47 HPRC Year 1 haplotypes + CHM13v2.0 +
    // GRCh38.p14, minus any shared with FDA-ARGOS bacteria or RefSeq viruses
    // (https://github.com/bede/deacon#prebuilt-indexes). Index format 3,
    // read by Deacon >= 0.11.0. The checksum is the Zenodo record's
    // (https://zenodo.org/records/17288185), which serves the same file.
    // --retry-all-errors because plain --retry skips a connection dropped
    // mid-transfer (curl exit 56, seen 2 GB into this 3.3 GB file), and
    // -C - so each retry carries on from the bytes already saved
    """
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    curl -fsSL --retry 5 --retry-all-errors --retry-delay 10 -C - -o panhuman-1.k31w15.idx \\
        "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w15.idx"

    echo "4bd0223e0619c26f1d6dc968896bc0df  panhuman-1.k31w15.idx" | md5sum -c -
    """
}
