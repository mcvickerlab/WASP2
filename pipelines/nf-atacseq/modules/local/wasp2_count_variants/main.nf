process WASP2_COUNT_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path vcf
    path peaks

    output:
    tuple val(meta), path("*_counts.tsv"), emit: counts
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = meta.sample_name ? "-s ${meta.sample_name}" : ""
    """
    wasp2-count count-variants \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        --region ${peaks} \\
        ${args} \\
        --out_file ${prefix}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-count --version 2>&1 | head -n1 || { echo "WARNING: Could not determine wasp2 version" >&2; echo "unknown"; })
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "chrom\\tpos\\tref\\talt\\tref_count\\talt_count\\tregion" > ${prefix}_counts.tsv
    echo -e "chr1\\t12345\\tA\\tG\\t10\\t8\\tpeak_1" >> ${prefix}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.2.1
    END_VERSIONS
    """
}
