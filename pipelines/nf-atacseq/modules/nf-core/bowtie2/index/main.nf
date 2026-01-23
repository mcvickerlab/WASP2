process BOWTIE2_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bowtie2=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.2--py39h6fed5c7_0' :
        'biocontainers/bowtie2:2.5.2--py39h6fed5c7_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bowtie2"), emit: index
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p bowtie2
    bowtie2-build \\
        $args \\
        --threads $task.cpus \\
        $fasta \\
        bowtie2/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = fasta.baseName
    """
    mkdir -p bowtie2
    touch bowtie2/${prefix}.{1,2,3,4}.bt2 bowtie2/${prefix}.rev.{1,2}.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.2
    END_VERSIONS
    """
}
