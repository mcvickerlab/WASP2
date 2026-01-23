process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path  multiqc_config
    path  extra_multiqc_config
    path  multiqc_logo

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "*_plots",              emit: plots, optional: true
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ""
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ""
    def logo = multiqc_logo ? "--cl-config 'custom_logo: $multiqc_logo'" : ""
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed 's/multiqc, version //' )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.19
    END_VERSIONS
    """
}
