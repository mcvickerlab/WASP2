process WASP2_MAKE_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path vcf

    output:
    tuple val(meta), path("*_remap_r1.fq.gz"), path("*_remap_r2.fq.gz"), emit: fastq
    tuple val(meta), path("*_to_remap.bam"),                             emit: to_remap_bam
    tuple val(meta), path("*_keep.bam"),                                 emit: keep_bam
    tuple val(meta), path("*_wasp_data_files.json"),                     emit: json
    path "versions.yml",                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = meta.sample_name ? "-s ${meta.sample_name}" : ""
    """
    wasp2-map make-reads \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        ${args} \\
        --out_dir . \\
        --out_json ${prefix}_wasp_data_files.json \\
        --threads ${task.cpus}

    # Rename outputs to include sample prefix (with validation)
    rename_single_file() {
        local pattern="\$1" target="\$2"
        local matches=(\$pattern)
        if [ \${#matches[@]} -eq 0 ] || [ ! -f "\${matches[0]}" ]; then
            echo "ERROR: No files matching pattern '\$pattern'" >&2
            exit 1
        fi
        if [ \${#matches[@]} -gt 1 ]; then
            echo "ERROR: Multiple files match '\$pattern': \${matches[*]}" >&2
            exit 1
        fi
        if [ "\${matches[0]}" != "\$target" ]; then
            mv "\${matches[0]}" "\$target" || exit 1
        fi
    }
    rename_single_file "*_remap_r1.fq.gz" "${prefix}_remap_r1.fq.gz"
    rename_single_file "*_remap_r2.fq.gz" "${prefix}_remap_r2.fq.gz"
    rename_single_file "*_to_remap.bam"   "${prefix}_to_remap.bam"
    rename_single_file "*_keep.bam"       "${prefix}_keep.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-map --version 2>&1 | head -n1 || { echo "WARNING: Could not determine wasp2 version" >&2; echo "unknown"; })
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_remap_r1.fq.gz
    touch ${prefix}_remap_r2.fq.gz
    touch ${prefix}_to_remap.bam
    touch ${prefix}_keep.bam
    echo '{}' > ${prefix}_wasp_data_files.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
    END_VERSIONS
    """
}
