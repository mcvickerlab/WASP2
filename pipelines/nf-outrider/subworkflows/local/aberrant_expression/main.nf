/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ABERRANT_EXPRESSION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTRIDER-based aberrant expression detection with WASP2 integration

    This subworkflow handles:
    1. Count matrix preparation
    2. OUTRIDER model fitting
    3. Outlier calling with multiple testing correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERGE_COUNTS   } from '../../modules/local/merge_counts'
include { OUTRIDER_FIT   } from '../../modules/local/outrider_fit'

workflow ABERRANT_EXPRESSION {
    take:
    ch_gene_counts   // channel: [ val(meta), path(gene_counts) ]
    padj_cutoff      // val: adjusted p-value cutoff
    zscore_cutoff    // val: z-score cutoff
    encoding_dim     // val: encoding dimension (null for auto)
    max_iterations   // val: max OUTRIDER iterations
    convergence      // val: convergence threshold

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Merge individual sample counts into matrix
    //
    MERGE_COUNTS(
        ch_gene_counts.map { meta, counts -> counts }.collect()
    )
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    //
    // MODULE: Fit OUTRIDER autoencoder and detect outliers
    //
    OUTRIDER_FIT(
        MERGE_COUNTS.out.count_matrix,
        padj_cutoff,
        zscore_cutoff,
        encoding_dim,
        max_iterations,
        convergence
    )
    ch_versions = ch_versions.mix(OUTRIDER_FIT.out.versions)

    emit:
    count_matrix = MERGE_COUNTS.out.count_matrix // channel: path(count_matrix)
    model        = OUTRIDER_FIT.out.model        // channel: path(model.rds)
    results      = OUTRIDER_FIT.out.results      // channel: path(results.tsv)
    summary      = OUTRIDER_FIT.out.summary      // channel: path(summary.html)
    versions     = ch_versions                    // channel: path(versions.yml)
}
