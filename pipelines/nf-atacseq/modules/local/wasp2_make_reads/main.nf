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

    # Rename R1 FASTQ (may be uncompressed .fq or compressed .fq.gz)
    r1_found=false
    for f in *_swapped_alleles_r1.fq.gz; do
        [ -f "\$f" ] && mv "\$f" ${prefix}_remap_r1.fq.gz && r1_found=true && break
    done
    if [ "\$r1_found" = "false" ]; then
        for f in *_swapped_alleles_r1.fq; do
            [ -f "\$f" ] && gzip -c "\$f" > ${prefix}_remap_r1.fq.gz && r1_found=true && break
        done
    fi

    # Rename R2 FASTQ
    r2_found=false
    for f in *_swapped_alleles_r2.fq.gz; do
        [ -f "\$f" ] && mv "\$f" ${prefix}_remap_r2.fq.gz && r2_found=true && break
    done
    if [ "\$r2_found" = "false" ]; then
        for f in *_swapped_alleles_r2.fq; do
            [ -f "\$f" ] && gzip -c "\$f" > ${prefix}_remap_r2.fq.gz && r2_found=true && break
        done
    fi

    # Rename to_remap and keep BAMs
    for f in *_to_remap.bam; do
        [ -f "\$f" ] && [ "\$f" != "${prefix}_to_remap.bam" ] && mv "\$f" ${prefix}_to_remap.bam && break
    done
    for f in *_keep.bam; do
        [ -f "\$f" ] && [ "\$f" != "${prefix}_keep.bam" ] && mv "\$f" ${prefix}_keep.bam && break
    done

    # Validate outputs
    for expected in ${prefix}_remap_r1.fq.gz ${prefix}_remap_r2.fq.gz ${prefix}_to_remap.bam ${prefix}_keep.bam ${prefix}_wasp_data_files.json; do
        if [ ! -f "\$expected" ]; then
            echo "ERROR: Expected output \$expected not found" >&2
            ls -la >&2
            exit 1
        fi
    done

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
