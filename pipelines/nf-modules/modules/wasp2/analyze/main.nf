process WASP2_ANALYZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.stats.tsv"), emit: stats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Sanitize inputs to prevent shell injection
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    def min_count = task.ext.min_count ?: 10
    def pseudocount = task.ext.pseudocount ?: 1
    def phased_arg = meta.phased ? '--phased' : ''
    """
    set -euo pipefail

    wasp2-analyze \\
        find-imbalance \\
        ${counts} \\
        --min ${min_count} \\
        --pseudocount ${pseudocount} \\
        -o ${prefix}.stats.tsv \\
        ${phased_arg} \\
        ${args}

    # Validate output was created
    if [ ! -f "${prefix}.stats.tsv" ]; then
        echo "ERROR: Output file ${prefix}.stats.tsv was not created" >&2
        exit 1
    fi

    # Get version - fail if not detectable
    WASP2_VERSION=\$(wasp2-analyze --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || true)
    if [ -z "\$WASP2_VERSION" ]; then
        echo "ERROR: Could not determine wasp2-analyze version. Tool may not be installed correctly." >&2
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
    cat <<-END_HEADER > ${prefix}.stats.tsv
    region	ref_count	alt_count	ratio	log2fc	pval	fdr_pval	dispersion
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
        python: 3.11.0
    END_VERSIONS
    """
}
