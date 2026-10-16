process WASP2_FIND_IMBALANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mcvickerlab/wasp2:1.5.0' :
        'ghcr.io/mcvickerlab/wasp2:1.5.0' }"

    input:
    tuple val(meta), path(counts)
    val min_count
    val pseudocount

    output:
    tuple val(meta), path("*_ai_results.tsv"), emit: results
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wasp2-analyze find-imbalance \\
        ${counts} \\
        --min ${min_count} \\
        --pseudocount ${pseudocount} \\
        ${args} \\
        --out_file ${prefix}_ai_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-analyze --version 2>&1 | head -n1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "region\\tref_count\\talt_count\\tN\\tsnp_count\\tnull_ll\\talt_ll\\tmu\\tlrt\\tpval\\tfdr_pval" > ${prefix}_ai_results.tsv
    echo -e "peak_1\\t10\\t8\\t18\\t1\\t-2.5\\t-2.4\\t0.55\\t0.2\\t0.655\\t0.655" >> ${prefix}_ai_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.5.0
    END_VERSIONS
    """
}
