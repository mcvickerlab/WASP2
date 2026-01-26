process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.11a bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721f691a6113293be77c7a1dff72e6a-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721f691a6113293be77c7a1dff72e6a-0' }"

    input:
    tuple val(meta), path(reads)
    path star_index
    path gtf

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), path("*.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("*.Log.final.out")                                                             , emit: log_final
    tuple val(meta), path("*.Log.out")                                                                   , emit: log_out
    tuple val(meta), path("*.SJ.out.tab")                                                                , emit: sj_tab
    path "versions.yml"                                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_arg = gtf ? "--sjdbGTFfile ${gtf}" : ''
    def read_files = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${read_files} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMattributes NH HI AS nM NM MD \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 999 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --twopassMode Basic \\
        ${gtf_arg} \\
        ${args}

    # Index BAM
    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam.bai
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.SJ.out.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11a
        samtools: 1.18
    END_VERSIONS
    """
}
