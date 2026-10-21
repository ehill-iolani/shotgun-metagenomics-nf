// Prebuilt Kraken2 databases from the Langmead lab
// (https://benlangmead.github.io/aws-indexes/k2), June 2026 build: RefSeq
// archaea, bacteria, viral, plasmid, human, UniVec_Core, protozoa and
// fungi. The -16/-8 builds are capped at 16/8 GB by keeping a subset of
// minimizers, so they're smaller and less sensitive than the full PlusPF
// (110 GB); PlusPF-8 runs on an ordinary laptop, for testing.
// md5s are the build's own, for the members kept (e.g.
// https://genome-idx.s3.amazonaws.com/kraken/pluspf_16_GB_20260626/pluspf_16_GB.md5).
def kraken2Dbs() {
    [
        pluspf16: [
            name: 'k2_pluspf_16_GB_20260626',
            md5: [
                'hash.k2d':  '429dce699dc0fdc72828597ab1eb0141',
                'opts.k2d':  '9c98c500a8f28f84f831b698437f4016',
                'taxo.k2d':  '626f4e8cd33cfaa0e1420987a17c08a0',
                'nodes.dmp': '88a1963036ad02403f63098884d88b21',
                'names.dmp': 'f416ac36abeb8098fc982f216972eb8d',
            ],
        ],
        pluspf8: [
            name: 'k2_pluspf_08_GB_20260626',
            md5: [
                'hash.k2d':  '52230b381ee3631fe73c1595c8169ad9',
                'opts.k2d':  '0bdf0ee1621b0c0783bf98c9cce28ae8',
                'taxo.k2d':  '626f4e8cd33cfaa0e1420987a17c08a0',
                'nodes.dmp': '88a1963036ad02403f63098884d88b21',
                'names.dmp': 'f416ac36abeb8098fc982f216972eb8d',
            ],
        ],
    ]
}

process FETCH_KRAKEN2_DB {
    label 'process_low'
    // same image and reason as FETCH_DEACON_INDEX
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    storeDir "${params.db_cache}/kraken2"

    input:
    val db   // a key of kraken2Dbs()

    output:
    path "${kraken2Dbs()[db].name}", emit: db

    script:
    def info = kraken2Dbs()[db]
    // Streamed straight into tar rather than saved first, so the archive
    // (11 GB for PlusPF-16) never sits on disk next to what's extracted
    // from it. Only the members classification needs are kept: the three
    // .k2d files, plus NCBI's taxonomy dump for TAXON_ABUNDANCE's ranks and
    // lineages (the Bracken .kmer_distrib files etc. are left out). The
    // checksums are of the extracted files, which also catches a download
    // cut short or corrupted mid-stream.
    """
    set -o pipefail
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    mkdir ${info.name}
    curl -fsSL "https://genome-idx.s3.amazonaws.com/kraken/${info.name}.tar.gz" \\
        | tar -xzf - -C ${info.name} ${info.md5.keySet().join(' ')}

    cd ${info.name}
    printf '%s\\n' ${info.md5.collect { f, sum -> "'${sum}  ${f}'" }.join(' ')} | md5sum -c -
    """
}
