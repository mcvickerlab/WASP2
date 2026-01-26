process WASP2_UNIFIED_MAKE_READS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_remap_r1.fq.gz"), path("*_remap_r2.fq.gz"), emit: remap_fastq
    tuple val(meta), path("*_to_remap.bam")                            , emit: to_remap_bam
    tuple val(meta), path("*_keep.bam")                                , emit: keep_bam
    tuple val(meta), path("*_wasp_data.json")                          , emit: wasp_json
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus ?: 4
    def sample_arg = meta.sample ? "--samples ${meta.sample}" : "--samples ${meta.id}"
    """
    # Run WASP2 make-reads to generate swapped allele FASTQs
    wasp2-map make-reads \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        --out_dir ./ \\
        --out_json ${prefix}_wasp_data.json \\
        --threads ${threads} \\
        ${args}

    # Rename outputs to standardized names with prefix
    for suffix in remap_r1.fq.gz remap_r2.fq.gz to_remap.bam keep.bam; do
        for f in *_\${suffix} \${suffix}; do
            [ -f "\$f" ] && mv "\$f" ${prefix}_\${suffix} && break
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "@read1" | gzip > ${prefix}_remap_r1.fq.gz
    echo "@read1" | gzip > ${prefix}_remap_r2.fq.gz
    touch ${prefix}_to_remap.bam
    touch ${prefix}_keep.bam
    echo '{}' > ${prefix}_wasp_data.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """
}
