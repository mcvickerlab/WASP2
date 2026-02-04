process WASP2_COUNT_ALLELES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_index)
    path gtf

    output:
    tuple val(meta), path("*_counts.tsv"), emit: counts
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = meta.sample ? "--samples ${meta.sample}" : "--samples ${meta.id}"
    def region_arg = gtf ? "--region ${gtf}" : ''
    def gene_feature = gtf ? "--gene_feature exon" : ''
    def gene_attr = gtf ? "--gene_attribute gene_id" : ''
    """
    wasp2-count count-variants \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        ${region_arg} \\
        ${gene_feature} \\
        ${gene_attr} \\
        --out_file ${prefix}_counts.tsv \\
        --use-rust \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(python -c "import wasp2; print(wasp2.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "chrom\\tpos\\tref\\talt\\tregion\\tref_count\\talt_count\\tother_count\\tN" > ${prefix}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: dev
    END_VERSIONS
    """
}
