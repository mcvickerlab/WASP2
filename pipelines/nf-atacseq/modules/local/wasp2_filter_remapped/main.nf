process WASP2_FILTER_REMAPPED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(remapped_bam), path(remapped_bai), path(to_remap_bam), path(keep_bam), path(wasp_json)

    output:
    tuple val(meta), path("*_wasp_filt.bam"), path("*_wasp_filt.bam.bai"), emit: bam
    tuple val(meta), path("*_wasp_stats.txt"),                             emit: stats
    path "versions.yml",                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Use JSON mode - WASP2 reads intermediate file paths from JSON
    wasp2-map filter-remapped \\
        ${remapped_bam} \\
        --json ${wasp_json} \\
        --threads ${task.cpus} \\
        ${args} \\
        --out_bam ${prefix}_wasp_filt.bam

    # Index the output BAM
    samtools index ${prefix}_wasp_filt.bam || { echo "ERROR: Failed to index BAM file" >&2; exit 1; }

    # Generate statistics
    samtools flagstat ${prefix}_wasp_filt.bam > ${prefix}_wasp_stats.txt || { echo "ERROR: Failed to generate flagstat" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-map --version 2>&1 | head -n1 || { echo "WARNING: Could not determine wasp2 version" >&2; echo "unknown"; })
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //' || { echo "WARNING: Could not determine samtools version" >&2; echo "unknown"; })
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_wasp_filt.bam
    touch ${prefix}_wasp_filt.bam.bai
    echo "WASP filtered reads" > ${prefix}_wasp_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
        samtools: 1.17
    END_VERSIONS
    """
}
