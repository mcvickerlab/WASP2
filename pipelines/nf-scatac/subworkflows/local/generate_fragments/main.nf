/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE_FRAGMENTS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generates 10x-compatible fragments.tsv.gz from scATAC-seq BAM files using sinto.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SINTO_FRAGMENTS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sinto=0.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sinto:0.9.0--pyhdfd78af_0' :
        'biocontainers/sinto:0.9.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.fragments.tsv.gz"), path("*.fragments.tsv.gz.tbi"), emit: fragments
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_tag = meta.cell_barcode_tag ?: 'CB'
    """
    set -o pipefail

    sinto fragments \\
        -b ${bam} \\
        -f ${prefix}.fragments.tsv \\
        --barcodetag ${barcode_tag} \\
        -p ${task.cpus} \\
        ${args}

    # Validate sinto output
    if [[ ! -s ${prefix}.fragments.tsv ]]; then
        echo "ERROR: sinto produced empty fragments file. Check BAM has '${barcode_tag}' tag." >&2
        exit 1
    fi

    sort -k1,1 -k2,2n ${prefix}.fragments.tsv > ${prefix}.fragments.sorted.tsv
    bgzip -c ${prefix}.fragments.sorted.tsv > ${prefix}.fragments.tsv.gz
    tabix -p bed ${prefix}.fragments.tsv.gz

    # Validate final output before cleanup
    if [[ ! -s ${prefix}.fragments.tsv.gz ]]; then
        echo "ERROR: bgzip output is empty" >&2
        exit 1
    fi

    rm ${prefix}.fragments.tsv ${prefix}.fragments.sorted.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version 2>&1 | sed 's/sinto //' || echo "unavailable")
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //' || echo "unavailable")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "chr1\\t100\\t500\\tAAACGAACAAGTCAGT-1\\t1" | bgzip > ${prefix}.fragments.tsv.gz
    tabix -p bed ${prefix}.fragments.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: 0.9.0
        tabix: 1.17
    END_VERSIONS
    """
}

workflow GENERATE_FRAGMENTS {
    take:
    ch_bam  // channel: [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()

    SINTO_FRAGMENTS ( ch_bam )
    ch_versions = ch_versions.mix(SINTO_FRAGMENTS.out.versions.first())

    emit:
    fragments = SINTO_FRAGMENTS.out.fragments  // channel: [ val(meta), path(fragments.tsv.gz), path(fragments.tbi) ]
    versions  = ch_versions                    // channel: [ path(versions.yml) ]
}
