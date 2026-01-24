/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCATAC WORKFLOW - Single-Cell ATAC-seq Allelic Imbalance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Allelic imbalance analysis on scATAC-seq data using WASP2.
    Accepts 10x-format fragments files and VCF with heterozygous SNPs.

    Features:
    - Cell barcode filtering (optional)
    - Peak region filtering (optional)
    - AnnData H5AD output for scverse ecosystem
    - Zarr output for GenVarLoader (optional)
    - Pseudo-bulk aggregation for statistical power
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASP_ALLELIC_SC } from '../subworkflows/local/wasp_allelic_sc/main'
include { WASP2_ML_OUTPUT } from '../../nf-modules/modules/wasp2/ml_output/main'

workflow SCATAC {
    take:
    samplesheet  // channel: [ val(meta), path(fragments.tsv.gz), path(fragments.tbi), path(barcodes), path(peaks) ]

    main:
    ch_versions = Channel.empty()
    ch_ml_zarr = Channel.empty()
    ch_ml_parquet = Channel.empty()
    ch_ml_anndata = Channel.empty()

    // Validate required VCF parameter
    if (!params.vcf) {
        error "ERROR: --vcf parameter is required. Provide path to indexed VCF/BCF with heterozygous SNPs."
    }

    // Validate VCF input and resolve index
    ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            def csi = file("${vcf}.csi")
            def idx = tbi.exists() ? tbi : (csi.exists() ? csi : null)
            if (!idx) {
                error "VCF index not found: expected ${vcf}.tbi or ${vcf}.csi. Run: tabix -p vcf ${vcf}"
            }
            [ [id: 'variants'], vcf, idx ]
        }

    // Run WASP2 single-cell allelic imbalance analysis
    // Includes: counting, AnnData creation, pseudo-bulk aggregation, statistical analysis
    WASP_ALLELIC_SC ( samplesheet, ch_vcf )
    ch_versions = ch_versions.mix(WASP_ALLELIC_SC.out.versions)

    // Convert to ML output formats (optional)
    if (params.output_format) {
        WASP2_ML_OUTPUT(
            WASP_ALLELIC_SC.out.cell_counts,
            params.output_format
        )
        ch_versions = ch_versions.mix(WASP2_ML_OUTPUT.out.versions)
        ch_ml_zarr = WASP2_ML_OUTPUT.out.zarr
        ch_ml_parquet = WASP2_ML_OUTPUT.out.parquet
        ch_ml_anndata = WASP2_ML_OUTPUT.out.anndata
    }

    emit:
    cell_counts = WASP_ALLELIC_SC.out.cell_counts  // channel: [ val(meta), path(counts.tsv) ]
    count_stats = WASP_ALLELIC_SC.out.count_stats  // channel: [ val(meta), path(stats.tsv) ]
    anndata     = WASP_ALLELIC_SC.out.anndata      // channel: [ val(meta), path(*.h5ad) ]
    zarr        = WASP_ALLELIC_SC.out.zarr         // channel: [ val(meta), path(*.zarr) ]
    cell_qc     = WASP_ALLELIC_SC.out.cell_qc      // channel: [ val(meta), path(cell_qc.tsv) ]
    pseudobulk  = WASP_ALLELIC_SC.out.pseudobulk   // channel: [ val(meta), path(pseudobulk.tsv) ]
    imbalance   = WASP_ALLELIC_SC.out.imbalance    // channel: [ val(meta), path(results.tsv) ]
    ml_zarr     = ch_ml_zarr                        // channel: [ val(meta), path(*.zarr) ] (ML format)
    ml_parquet  = ch_ml_parquet                     // channel: [ val(meta), path(*.parquet) ]
    ml_anndata  = ch_ml_anndata                     // channel: [ val(meta), path(*.h5ad) ] (ML format)
    versions    = ch_versions                       // channel: [ path(versions.yml) ]
}
