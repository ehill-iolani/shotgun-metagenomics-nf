// Kraken2 database built from local genomes, for --kraken2_build_fasta --
// the test profile's synthetic microbes and their made-up taxonomy, so it
// runs offline instead of downloading PlusPF. The output has the same
// files as FETCH_KRAKEN2_DB's, so the steps after it can't tell them apart.
process BUILD_KRAKEN2_DB {
    label 'process_low'
    container "${params.container_registry}/biocontainers/kraken2:2.17.2--pl5321h3be2455_0"

    input:
    path fastas
    path taxonomy   // dir: nodes.dmp, names.dmp, seqid2taxid.map (FASTA ID -> taxid)

    output:
    path "kraken2_db", emit: db

    script:
    // kraken2-build reads each sequence's taxid from its header
    // (>id|kraken:taxid|N), so add it there from seqid2taxid.map; a
    // sequence missing from the map fails the build rather than being left
    // out of it
    """
    set -o pipefail

    zcat -f ${fastas} \\
      | awk 'NR == FNR { taxid[\$1] = \$2; next }
             /^>/ {
                 id = substr(\$1, 2)
                 if (!(id in taxid)) { print "no taxid in seqid2taxid.map for " id > "/dev/stderr"; exit 1 }
                 print ">" id "|kraken:taxid|" taxid[id]
                 next
             }
             { print }' ${taxonomy}/seqid2taxid.map - \\
      > library.fa

    mkdir -p build/taxonomy kraken2_db
    cp ${taxonomy}/nodes.dmp ${taxonomy}/names.dmp build/taxonomy/
    kraken2-build --db build --add-to-library library.fa
    kraken2-build --db build --build --threads ${task.cpus}

    mv build/*.k2d kraken2_db/
    cp ${taxonomy}/nodes.dmp ${taxonomy}/names.dmp kraken2_db/
    """
}
