process WASP2_FIND_IMBALANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

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
    # WASP2 v1.4.0 workaround: normalize count-variants output to
    # the 9-column format expected by find-imbalance Rust backend:
    #   chrom  pos0  pos  ref  alt  GT  ref_count  alt_count  other_count
    # count-variants outputs (with --region):
    #   chrom  pos  ref  alt  GT  region  ref_count  alt_count  other_count
    # We strip 'region' and insert 'pos0' (= pos - 1).
    if head -1 ${counts} | grep -q 'region'; then
        # Strip region column
        REGION_COL=\$(head -1 ${counts} | tr '\\t' '\\n' | grep -n 'region' | cut -d: -f1)
        cut -f1-\$((REGION_COL-1)),\$((REGION_COL+1))- ${counts} > ${prefix}_counts_stripped.tsv
    else
        cp ${counts} ${prefix}_counts_stripped.tsv
    fi

    # Insert pos0 column if missing (count-variants lacks it)
    if ! head -1 ${prefix}_counts_stripped.tsv | grep -q 'pos0'; then
        awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{print \$1,"pos0",\$0; next} {print \$1,\$2-1,\$0}' ${prefix}_counts_stripped.tsv | cut -f1,2,4- > ${prefix}_counts_fixed.tsv
        COUNTS_INPUT=${prefix}_counts_fixed.tsv
    else
        COUNTS_INPUT=${prefix}_counts_stripped.tsv
    fi

    wasp2-analyze find-imbalance \\
        \$COUNTS_INPUT \\
        --min ${min_count} \\
        --pseudocount ${pseudocount} \\
        ${args} \\
        --out_file ${prefix}_ai_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-analyze --version 2>&1 | head -n1 || { echo "WARNING: Could not determine wasp2 version" >&2; echo "unknown"; })
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "region\\tref_count\\talt_count\\tpval\\tfdr_pval\\tlog2_ratio" > ${prefix}_ai_results.tsv
    echo -e "peak_1\\t10\\t8\\t0.05\\t0.1\\t0.322" >> ${prefix}_ai_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
    END_VERSIONS
    """
}
