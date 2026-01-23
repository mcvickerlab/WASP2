process BWA_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa=0.7.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_0' :
        'biocontainers/bwa:0.7.18--he4a0461_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwa
    bwa index $args -p bwa/${fasta.baseName} $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -o 'Version: [0-9.]*' | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa
    touch bwa/${fasta.baseName}.amb
    touch bwa/${fasta.baseName}.ann
    touch bwa/${fasta.baseName}.bwt
    touch bwa/${fasta.baseName}.pac
    touch bwa/${fasta.baseName}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.18
    END_VERSIONS
    """
}
