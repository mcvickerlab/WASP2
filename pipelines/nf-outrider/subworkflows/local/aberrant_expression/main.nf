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
    // Parameter validation
    //
    if (padj_cutoff <= 0 || padj_cutoff >= 1) {
        error "ERROR: padj_cutoff must be between 0 and 1 (exclusive), got: ${padj_cutoff}"
    }
    if (zscore_cutoff <= 0) {
        error "ERROR: zscore_cutoff must be positive, got: ${zscore_cutoff}"
    }
    if (max_iterations <= 0) {
        error "ERROR: max_iterations must be positive, got: ${max_iterations}"
    }
    if (convergence <= 0) {
        error "ERROR: convergence threshold must be positive, got: ${convergence}"
    }

    //
    // Validate minimum sample count for OUTRIDER
    //
    ch_gene_counts
        .count()
        .map { n ->
            if (n < 15) {
                log.warn "WARNING: OUTRIDER requires >= 15 samples for reliable results. Found ${n} samples."
            }
        }

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
