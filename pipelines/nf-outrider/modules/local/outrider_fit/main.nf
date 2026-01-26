process OUTRIDER_FIT {
    tag "outrider"
    label 'process_high'
    label 'process_high_memory'

    conda "bioconda::bioconductor-outrider=1.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/gagneurlab/outrider:latest' :
        'ghcr.io/gagneurlab/outrider:latest' }"

    input:
    path count_matrix
    val padj_cutoff
    val zscore_cutoff
    val encoding_dim
    val max_iterations
    val convergence

    output:
    path "outrider_model.rds"   , emit: model
    path "outrider_results.tsv" , emit: results
    path "outrider_summary.html", emit: summary, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    Rscript ${projectDir}/bin/runOutrider.R \\
        --counts ${count_matrix} \\
        --output_model outrider_model.rds \\
        --output_results outrider_results.tsv \\
        --padj ${padj_cutoff} \\
        --zscore ${zscore_cutoff} \\
        ${encoding_dim ? "--q ${encoding_dim}" : ""} \\
        --iterations ${max_iterations} \\
        --convergence ${convergence} \\
        --threads ${task.cpus}

    # Validate output was created
    if [ ! -f "outrider_results.tsv" ]; then
        echo "ERROR: OUTRIDER results file was not created" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //' | cut -d' ' -f1)
        OUTRIDER: \$(Rscript -e "cat(as.character(packageVersion('OUTRIDER')))")
    END_VERSIONS
    """

    stub:
    """
    touch outrider_model.rds
    cat <<-END_HEADER > outrider_results.tsv
    geneID	sampleID	pValue	padjust	zScore	l2fc	rawcounts	normcounts	meanCorrected	theta	aberrant
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.0
        OUTRIDER: 1.16.0
    END_VERSIONS
    """
}
