/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Local modules
include { STAR_ALIGN as STAR_ALIGN_INITIAL } from '../modules/local/star_align/main'
include { WASP2_COUNT_ALLELES              } from '../modules/local/wasp2_count_alleles/main'
include { WASP2_ANALYZE_IMBALANCE          } from '../modules/local/wasp2_analyze_imbalance/main'
include { WASP2_ML_OUTPUT                  } from '../modules/local/wasp2_ml_output/main'

// Local subworkflows
include { WASP_RNASEQ_MAPPING              } from '../subworkflows/local/wasp_rnaseq_mapping/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ_ASE {

    take:
    ch_samplesheet // channel: [ val(meta), [ fastq ] ]
    ch_vcf         // channel: [ val(meta), path(vcf), path(vcf_index) ]

    main:
    // Channel for version tracking
    ch_versions = Channel.empty()

    // Initialize output channels for conditional steps
    ch_ai_results = Channel.empty()
    ch_ml_zarr = Channel.empty()
    ch_ml_parquet = Channel.empty()
    ch_ml_anndata = Channel.empty()

    //
    // Load reference files
    //
    ch_star_index = file(params.star_index)
    ch_gtf = params.gtf ? file(params.gtf) : []

    //
    // STEP 1: Initial STAR alignment
    //
    STAR_ALIGN_INITIAL(
        ch_samplesheet,
        ch_star_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_INITIAL.out.versions)

    //
    // STEP 2-4: WASP mapping bias correction
    // Includes: make_reads -> remap -> filter
    //
    WASP_RNASEQ_MAPPING(
        STAR_ALIGN_INITIAL.out.bam,
        ch_vcf,
        ch_star_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(WASP_RNASEQ_MAPPING.out.versions)

    //
    // STEP 5: Count alleles at heterozygous SNPs
    //
    WASP2_COUNT_ALLELES(
        WASP_RNASEQ_MAPPING.out.bam,
        ch_vcf.first(),
        ch_gtf
    )
    ch_versions = ch_versions.mix(WASP2_COUNT_ALLELES.out.versions)

    //
    // STEP 6: Statistical testing for allelic imbalance (optional)
    // Skip if params.skip_analysis is true
    //
    if (!params.skip_analysis) {
        WASP2_ANALYZE_IMBALANCE(
            WASP2_COUNT_ALLELES.out.counts
        )
        ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions)
        ch_ai_results = WASP2_ANALYZE_IMBALANCE.out.results
    }

    //
    // STEP 7: Convert to ML output formats (optional)
    // Run if params.output_format is specified
    //
    if (params.output_format) {
        WASP2_ML_OUTPUT(
            WASP2_COUNT_ALLELES.out.counts,
            params.output_format
        )
        ch_versions = ch_versions.mix(WASP2_ML_OUTPUT.out.versions)
        ch_ml_zarr = WASP2_ML_OUTPUT.out.zarr
        ch_ml_parquet = WASP2_ML_OUTPUT.out.parquet
        ch_ml_anndata = WASP2_ML_OUTPUT.out.anndata
    }

    emit:
    wasp_bam    = WASP_RNASEQ_MAPPING.out.bam     // channel: [ val(meta), path(bam), path(bai) ]
    counts      = WASP2_COUNT_ALLELES.out.counts  // channel: [ val(meta), path(counts) ]
    results     = ch_ai_results                    // channel: [ val(meta), path(results) ]
    ml_zarr     = ch_ml_zarr                       // channel: [ val(meta), path(zarr) ]
    ml_parquet  = ch_ml_parquet                    // channel: [ val(meta), path(parquet) ]
    ml_anndata  = ch_ml_anndata                    // channel: [ val(meta), path(anndata) ]
    versions    = ch_versions                      // channel: path(versions.yml)
}
