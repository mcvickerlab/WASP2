process WASP2_ANALYZE_IMBALANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jaureguy760/wasp2:latest' :
        'jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*_ai_results.tsv"), emit: results
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_count = params.min_count ?: 10
    def pseudocount = params.pseudocount ?: 1
    """
    wasp2-analyze find-imbalance \\
        ${counts} \\
        --out_file ${prefix}_ai_results.tsv \\
        --min ${min_count} \\
        --pseudocount ${pseudocount} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "region\\tsnp_count\\tref_sum\\talt_sum\\tmu\\tnull_ll\\talt_ll\\tLRT\\tpvalue\\tfdr" > ${prefix}_ai_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: dev
    END_VERSIONS
    """
}
