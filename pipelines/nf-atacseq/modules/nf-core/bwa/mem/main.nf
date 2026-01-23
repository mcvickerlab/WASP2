process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3571' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3571' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  fasta
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"

    def samtools_command = sort_bam ? "samtools sort -@ ${task.cpus} -o ${prefix}.bam -" : "samtools view -@ ${task.cpus} $args2 -o ${prefix}.bam -"

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        -R "$read_group" \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | $samtools_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -o 'Version: [0-9.]*' | sed 's/Version: //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.18
        samtools: 1.19
    END_VERSIONS
    """
}
