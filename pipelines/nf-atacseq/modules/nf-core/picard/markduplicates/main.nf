process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.markdup.bam"),    emit: bam
    tuple val(meta), path("*.markdup.bam.bai"), emit: bai
    tuple val(meta), path("*.metrics.txt"),    emit: metrics
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.markdup.bam \\
        --METRICS_FILE ${prefix}.metrics.txt \\
        --CREATE_INDEX true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:[0-9.]*' | sed 's/Version://')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.markdup.bam
    touch ${prefix}.markdup.bam.bai
    touch ${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.1.1
    END_VERSIONS
    """
}
