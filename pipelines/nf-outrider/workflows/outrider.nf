/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTRIDER WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP2 + OUTRIDER for Aberrant Expression and MAE Detection

    Workflow:
    RNA-seq BAMs → WASP2 Filter → WASP2 Count → Gene Aggregation → OUTRIDER → Outlier Calls
                        ↓                              ↓
                  Bias-corrected             Autoencoder-based
                      reads                  outlier detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WASP2 modules from shared nf-modules
include { WASP2_COUNT             } from '../../nf-modules/modules/wasp2/count/main'
include { WASP2_ML_OUTPUT         } from '../../nf-modules/modules/wasp2/ml_output/main'

// Local modules
include { AGGREGATE_COUNTS        } from '../modules/local/aggregate_counts'
include { MERGE_COUNTS            } from '../modules/local/merge_counts'
include { OUTRIDER_FIT            } from '../modules/local/outrider_fit'
include { MAE_DETECT              } from '../modules/local/mae_detect'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OUTRIDER {

    take:
    ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions    = Channel.empty()
    ch_counts      = Channel.empty()
    ch_gene_counts = Channel.empty()
    ch_outliers    = Channel.empty()
    ch_mae_results = Channel.empty()
    ch_ml_zarr     = Channel.empty()
    ch_ml_parquet  = Channel.empty()
    ch_ml_anndata  = Channel.empty()

    //
    // Prepare reference files
    //
    ch_vcf = params.vcf ? Channel.fromPath(params.vcf, checkIfExists: true).collect() : Channel.empty()
    ch_gtf = params.gtf ? Channel.fromPath(params.gtf, checkIfExists: true).collect() : Channel.empty()

    // Validate required inputs
    if (!params.vcf) {
        error "ERROR: --vcf is required for WASP2 allele counting"
    }
    if (!params.gtf) {
        error "ERROR: --gtf is required for gene aggregation"
    }

    //
    // STEP 1: WASP2 Allele Counting
    //
    // Count alleles at heterozygous variants using WASP2's Rust-accelerated counter
    // This replaces GATK ASEReadCounter with 61× faster processing
    //

    // Prepare VCF channel with meta for WASP2_COUNT
    ch_vcf_with_meta = ch_vcf.map { vcf ->
        def vcf_index = file("${vcf}.tbi").exists() ? file("${vcf}.tbi") :
                        (file("${vcf}.csi").exists() ? file("${vcf}.csi") : null)
        if (!vcf_index) {
            error "VCF index not found for ${vcf}. Please provide .tbi or .csi index."
        }
        [[id: 'variants'], vcf, vcf_index]
    }

    WASP2_COUNT(
        ch_samplesheet,           // [ meta, bam, bai ]
        ch_vcf_with_meta.first(), // [ meta, vcf, vcf_index ]
        ch_gtf.first()            // GTF for region filtering
    )
    ch_counts = WASP2_COUNT.out.counts
    ch_versions = ch_versions.mix(WASP2_COUNT.out.versions.first())

    //
    // STEP 2: Gene Aggregation
    //
    // Aggregate variant-level counts to gene level for OUTRIDER input
    //

    AGGREGATE_COUNTS(
        ch_counts,
        ch_gtf.first(),
        params.aggregation_method
    )
    ch_gene_counts = AGGREGATE_COUNTS.out.gene_counts
    ch_versions = ch_versions.mix(AGGREGATE_COUNTS.out.versions.first())

    //
    // STEP 3: Merge Sample Counts
    //
    // Create gene x sample count matrix for OUTRIDER
    //

    MERGE_COUNTS(
        ch_gene_counts.map { meta, counts -> counts }.collect()
    )
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    //
    // STEP 4: OUTRIDER Aberrant Expression Detection
    //
    // Fit autoencoder model and detect expression outliers
    //

    OUTRIDER_FIT(
        MERGE_COUNTS.out.count_matrix,
        params.outrider_padj,
        params.outrider_zScore,
        params.outrider_q,
        params.outrider_iterations,
        params.outrider_convergence
    )
    ch_outliers = OUTRIDER_FIT.out.results
    ch_versions = ch_versions.mix(OUTRIDER_FIT.out.versions)

    //
    // STEP 5: MAE Detection (Optional)
    //
    // Detect mono-allelic expression using binomial test with FDR correction
    //

    if (!params.skip_mae) {
        MAE_DETECT(
            ch_counts,
            params.mae_min_count,
            params.mae_padj,
            params.mae_alt_ratio
        )
        ch_mae_results = MAE_DETECT.out.mae_results
        ch_versions = ch_versions.mix(MAE_DETECT.out.versions.first())
    }

    //
    // STEP 6: ML Output Formats (Optional)
    //
    // Convert counts to ML-ready formats: Zarr, Parquet, AnnData
    //

    if (params.output_format) {
        WASP2_ML_OUTPUT(
            ch_counts,
            params.output_format
        )
        ch_versions = ch_versions.mix(WASP2_ML_OUTPUT.out.versions.first())
        ch_ml_zarr = WASP2_ML_OUTPUT.out.zarr
        ch_ml_parquet = WASP2_ML_OUTPUT.out.parquet
        ch_ml_anndata = WASP2_ML_OUTPUT.out.anndata
    }

    emit:
    counts      = ch_counts                                    // channel: [ val(meta), path(counts) ]
    gene_counts = ch_gene_counts                               // channel: [ val(meta), path(gene_counts) ]
    outliers    = ch_outliers                                  // channel: path(outliers.tsv)
    mae_results = params.skip_mae ? Channel.empty() : ch_mae_results // channel: [ val(meta), path(mae_results) ]
    ml_zarr     = ch_ml_zarr                                   // channel: [ val(meta), path(*.zarr) ]
    ml_parquet  = ch_ml_parquet                                // channel: [ val(meta), path(*.parquet) ]
    ml_anndata  = ch_ml_anndata                                // channel: [ val(meta), path(*.h5ad) ]
    versions    = ch_versions                                  // channel: path(versions.yml)
}
