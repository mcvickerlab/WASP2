process WASP2_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_arg = meta.sample ? "-s ${meta.sample}" : ''
    def region_arg = regions ? "-r ${regions}" : ''
    def use_rust = task.ext.use_rust != false ? '--use-rust' : '--no-rust'
    """
    wasp2-count \\
        count-variants \\
        ${bam} \\
        ${vcf} \\
        ${sample_arg} \\
        ${region_arg} \\
        -o ${prefix}.counts.tsv \\
        ${use_rust} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-count --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' || echo 'unknown')
        python: \$(python --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
