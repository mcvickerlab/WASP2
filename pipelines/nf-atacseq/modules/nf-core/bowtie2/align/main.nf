process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bowtie2=2.5.2 bioconda::samtools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6' :
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  fasta
    val   save_unaligned
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"),     emit: aligned
    tuple val(meta), path("*.log"),     emit: log
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_inputs = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command = sort_bam ? "samtools sort -@ ${task.cpus} -o ${prefix}.bam -" : "samtools view -@ ${task.cpus} -bS -o ${prefix}.bam -"

    """
    INDEX=`find -L ./ -name "*.1.bt2" | sed 's/\\.1.bt2\$//'`

    bowtie2 \\
        $args \\
        --threads $task.cpus \\
        -x \$INDEX \\
        $read_inputs \\
        2> ${prefix}.bowtie2.log \\
        | $samtools_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n1 | sed 's/.*version //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bowtie2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.2
        samtools: 1.19
    END_VERSIONS
    """
}
