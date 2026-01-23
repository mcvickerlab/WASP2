process WASP2_ANALYZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/jaureguy760/wasp2:latest' :
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_count = task.ext.min_count ?: 10
    def pseudocount = task.ext.pseudocount ?: 1
    def phased_arg = meta.phased ? '--phased' : ''
    """
    wasp2-analyze \\
        find-imbalance \\
        ${counts} \\
        --min ${min_count} \\
        --pseudocount ${pseudocount} \\
        -o ${prefix}.stats.tsv \\
        ${phased_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-analyze --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || echo 'unknown')
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
