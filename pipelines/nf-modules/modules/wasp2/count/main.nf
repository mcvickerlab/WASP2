process WASP2_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_index)
    path regions  // Optional: BED, GTF, GFF3, or narrowPeak/broadPeak file

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Sanitize inputs to prevent shell injection
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    def sample = meta.sample?.toString()?.replaceAll(/[^a-zA-Z0-9._-]/, '_') ?: ''
    def sample_arg = sample ? "-s ${sample}" : ''
    def region_arg = regions ? "-r ${regions}" : ''
    def use_rust = task.ext.use_rust != false ? '--use-rust' : '--no-rust'
    """
    set -euo pipefail

    wasp2-count \\
        count-variants \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        ${region_arg} \\
        -o ${prefix}.counts.tsv \\
        ${use_rust} \\
        ${args}

    # Validate output was created and has content
    if [ ! -f "${prefix}.counts.tsv" ]; then
        echo "ERROR: Output file ${prefix}.counts.tsv was not created" >&2
        exit 1
    fi

    # Get version - fail if not detectable
    WASP2_VERSION=\$(wasp2-count --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || true)
    if [ -z "\$WASP2_VERSION" ]; then
        echo "ERROR: Could not determine wasp2-count version. Tool may not be installed correctly." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \${WASP2_VERSION}
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    cat <<-END_HEADER > ${prefix}.counts.tsv
    chrom	pos	ref	alt	ref_count	alt_count	other_count
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
        python: 3.11.0
    END_VERSIONS
    """
}
