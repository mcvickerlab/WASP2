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

workflow SCATAC {
    take:
    samplesheet  // channel: [ val(meta), path(fragments.tsv.gz), path(fragments.tbi), path(barcodes), path(peaks) ]

    main:
    ch_versions = Channel.empty()

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

    emit:
    cell_counts = WASP_ALLELIC_SC.out.cell_counts  // channel: [ val(meta), path(counts.tsv) ]
    count_stats = WASP_ALLELIC_SC.out.count_stats  // channel: [ val(meta), path(stats.tsv) ]
    anndata     = WASP_ALLELIC_SC.out.anndata      // channel: [ val(meta), path(*.h5ad) ]
    zarr        = WASP_ALLELIC_SC.out.zarr         // channel: [ val(meta), path(*.zarr) ]
    cell_qc     = WASP_ALLELIC_SC.out.cell_qc      // channel: [ val(meta), path(cell_qc.tsv) ]
    pseudobulk  = WASP_ALLELIC_SC.out.pseudobulk   // channel: [ val(meta), path(pseudobulk.tsv) ]
    imbalance   = WASP_ALLELIC_SC.out.imbalance    // channel: [ val(meta), path(results.tsv) ]
    versions    = ch_versions                       // channel: [ path(versions.yml) ]
}
