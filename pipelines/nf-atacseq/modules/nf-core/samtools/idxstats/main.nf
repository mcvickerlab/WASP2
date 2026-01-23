process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19--h50ea8bc_0' :
        'biocontainers/samtools:1.19--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools idxstats \\
        $args \\
        $bam \\
        > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.19
    END_VERSIONS
    """
}
