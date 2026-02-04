process WASP2_MAP_MAKE_READS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.remap_r1.fq.gz"), path("*.remap_r2.fq.gz"), emit: fastq
    tuple val(meta), path("*.to_remap.bam")                            , emit: to_remap_bam
    tuple val(meta), path("*.keep.bam")                                , emit: keep_bam
    tuple val(meta), path("*.wasp_data.json")                          , emit: wasp_json
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Sanitize inputs to prevent shell injection
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    def sample = meta.sample?.toString()?.replaceAll(/[^a-zA-Z0-9._-]/, '_') ?: ''
    def sample_arg = sample ? "-s ${sample}" : ''
    """
    set -euo pipefail

    wasp2-map \\
        make-reads \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        -o . \\
        -j ${prefix}.wasp_data.json \\
        --threads ${task.cpus} \\
        --phased \\
        ${args}

    # Rename R1 FASTQ - explicit error checking
    r1_found=false
    for f in *_swapped_alleles_r1.fq.gz; do
        if [ -f "\$f" ]; then
            mv "\$f" ${prefix}.remap_r1.fq.gz
            r1_found=true
            break
        fi
    done
    if [ "\$r1_found" = "false" ]; then
        for f in *_swapped_alleles_r1.fq; do
            if [ -f "\$f" ]; then
                gzip -c "\$f" > ${prefix}.remap_r1.fq.gz
                r1_found=true
                break
            fi
        done
    fi

    # Rename R2 FASTQ - explicit error checking
    r2_found=false
    for f in *_swapped_alleles_r2.fq.gz; do
        if [ -f "\$f" ]; then
            mv "\$f" ${prefix}.remap_r2.fq.gz
            r2_found=true
            break
        fi
    done
    if [ "\$r2_found" = "false" ]; then
        for f in *_swapped_alleles_r2.fq; do
            if [ -f "\$f" ]; then
                gzip -c "\$f" > ${prefix}.remap_r2.fq.gz
                r2_found=true
                break
            fi
        done
    fi

    # Rename to_remap BAM
    remap_found=false
    for f in *_to_remap.bam; do
        if [ -f "\$f" ]; then
            mv "\$f" ${prefix}.to_remap.bam
            remap_found=true
            break
        fi
    done

    # Rename keep BAM
    keep_found=false
    for f in *_keep.bam; do
        if [ -f "\$f" ]; then
            mv "\$f" ${prefix}.keep.bam
            keep_found=true
            break
        fi
    done

    # Validate ALL required outputs exist
    missing_files=""
    if [ ! -f "${prefix}.remap_r1.fq.gz" ]; then
        missing_files="\${missing_files} remap_r1.fq.gz"
    fi
    if [ ! -f "${prefix}.remap_r2.fq.gz" ]; then
        missing_files="\${missing_files} remap_r2.fq.gz"
    fi
    if [ ! -f "${prefix}.to_remap.bam" ]; then
        missing_files="\${missing_files} to_remap.bam"
    fi
    if [ ! -f "${prefix}.keep.bam" ]; then
        missing_files="\${missing_files} keep.bam"
    fi
    if [ ! -f "${prefix}.wasp_data.json" ]; then
        missing_files="\${missing_files} wasp_data.json"
    fi

    if [ -n "\$missing_files" ]; then
        echo "ERROR: Required output files not found:\$missing_files" >&2
        exit 1
    fi

    # Get version - fail if not detectable
    WASP2_VERSION=\$(wasp2-map --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || true)
    if [ -z "\$WASP2_VERSION" ]; then
        echo "ERROR: Could not determine wasp2-map version" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \${WASP2_VERSION}
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        samtools: \$(samtools --version 2>&1 | head -1 | grep -oE '[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    echo "" | gzip > ${prefix}.remap_r1.fq.gz
    echo "" | gzip > ${prefix}.remap_r2.fq.gz
    touch ${prefix}.to_remap.bam
    touch ${prefix}.keep.bam
    echo '{"bam_file": "test.bam", "variant_file": "test.vcf", "read_mappings": {}}' > ${prefix}.wasp_data.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
        python: 3.11.0
        samtools: 1.17
    END_VERSIONS
    """
}

process WASP2_MAP_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(remapped_bam), path(remapped_bai)
    tuple val(meta2), path(to_remap_bam)
    tuple val(meta3), path(keep_bam)
    tuple val(meta4), path(wasp_json)

    output:
    tuple val(meta), path("*.wasp_filt.bam"), path("*.wasp_filt.bam.bai"), emit: bam
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Sanitize inputs to prevent shell injection
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    def use_rust = task.ext.use_rust != false ? '--use-rust' : '--no-rust'
    """
    set -euo pipefail

    wasp2-map \\
        filter-remapped \\
        ${remapped_bam} \\
        -j ${wasp_json} \\
        -o ${prefix}.wasp_filt.bam \\
        --threads ${task.cpus} \\
        ${use_rust} \\
        ${args}

    # Validate output BAM was created
    if [ ! -f "${prefix}.wasp_filt.bam" ]; then
        echo "ERROR: Output BAM ${prefix}.wasp_filt.bam was not created" >&2
        exit 1
    fi

    # Index the output BAM
    if [ ! -f "${prefix}.wasp_filt.bam.bai" ]; then
        samtools index ${prefix}.wasp_filt.bam || {
            echo "ERROR: Failed to index ${prefix}.wasp_filt.bam" >&2
            exit 1
        }
    fi

    # Validate index was created
    if [ ! -f "${prefix}.wasp_filt.bam.bai" ]; then
        echo "ERROR: BAM index ${prefix}.wasp_filt.bam.bai was not created" >&2
        exit 1
    fi

    # Get version - fail if not detectable
    WASP2_VERSION=\$(wasp2-map --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || true)
    if [ -z "\$WASP2_VERSION" ]; then
        echo "ERROR: Could not determine wasp2-map version" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \${WASP2_VERSION}
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        samtools: \$(samtools --version 2>&1 | head -1 | grep -oE '[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    touch ${prefix}.wasp_filt.bam
    touch ${prefix}.wasp_filt.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
        python: 3.11.0
        samtools: 1.17
    END_VERSIONS
    """
}
