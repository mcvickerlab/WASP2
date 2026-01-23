process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::macs2=2.2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs2:2.2.9.1--py39hf95cd2a_0' :
        'biocontainers/macs2:2.2.9.1--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)
    val   gsize

    output:
    tuple val(meta), path("*.narrowPeak"),  emit: peak
    tuple val(meta), path("*.xls"),         emit: xls
    tuple val(meta), path("*.summits.bed"), emit: summits, optional: true
    tuple val(meta), path("*.bdg"),         emit: bedgraph, optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = meta.single_end ? 'BAM' : 'BAMPE'
    """
    macs2 callpeak \\
        $args \\
        -g $gsize \\
        -f $format \\
        -t $bam \\
        -n $prefix \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version 2>&1 | sed 's/macs2 //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_peaks.narrowPeak
    touch ${prefix}_peaks.xls
    touch ${prefix}_summits.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: 2.2.9.1
    END_VERSIONS
    """
}
