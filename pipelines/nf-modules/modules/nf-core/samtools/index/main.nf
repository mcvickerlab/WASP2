process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19--h50ea8bc_0' :
        'biocontainers/samtools:1.19--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools index $args -@ $task.cpus $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.19
    END_VERSIONS
    """
}
