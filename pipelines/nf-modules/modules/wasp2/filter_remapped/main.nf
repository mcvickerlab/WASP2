process WASP2_FILTER_REMAPPED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(remapped_bam), path(remapped_bai)
    tuple val(_meta2), path(to_remap_bam)
    tuple val(_meta3), path(keep_bam)
    tuple val(_meta4), path(wasp_json)

    output:
    tuple val(meta), path("*_wasp_filt.bam"), path("*_wasp_filt.bam.bai"), emit: bam
    tuple val(meta), path("*.filter_stats.txt")                          , emit: stats, optional: true
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 4
    """
    # Filter remapped reads using WASP algorithm
    # Note: use positional args for BAMs, not --json (which overrides positional args)
    wasp2-map filter-remapped \\
        ${remapped_bam} \\
        ${to_remap_bam} \\
        ${keep_bam} \\
        --out_bam ${prefix}_remapped_filt.bam \\
        --threads ${threads} \\
        --use-rust \\
        ${args}

    # Merge filtered remapped reads with keep reads
    samtools merge \\
        -@ ${threads} \\
        -f \\
        ${prefix}_wasp_filt.bam \\
        ${keep_bam} \\
        ${prefix}_remapped_filt.bam

    # Sort and index final BAM
    samtools sort -@ ${threads} -o ${prefix}_wasp_filt_sorted.bam ${prefix}_wasp_filt.bam
    mv ${prefix}_wasp_filt_sorted.bam ${prefix}_wasp_filt.bam
    samtools index -@ ${threads} ${prefix}_wasp_filt.bam

    # Generate filter statistics
    echo "Sample: ${prefix}" > ${prefix}.filter_stats.txt
    echo "Total reads in remapped BAM: \$(samtools view -c ${remapped_bam})" >> ${prefix}.filter_stats.txt
    echo "Total reads in keep BAM: \$(samtools view -c ${keep_bam})" >> ${prefix}.filter_stats.txt
    echo "Total reads after WASP filter: \$(samtools view -c ${prefix}_wasp_filt.bam)" >> ${prefix}.filter_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_wasp_filt.bam
    touch ${prefix}_wasp_filt.bam.bai
    echo "stub" > ${prefix}.filter_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: dev
        samtools: 1.18
    END_VERSIONS
    """
}
