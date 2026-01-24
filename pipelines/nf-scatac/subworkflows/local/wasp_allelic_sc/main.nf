/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP_ALLELIC_SC SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP2 single-cell allelic imbalance analysis for scATAC-seq data.
    Includes per-cell counting, AnnData output, and optional pseudo-bulk analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASP2_VCF_TO_BED        } from '../../../../nf-modules/modules/wasp2/vcf_to_bed/main'
include { SCATAC_COUNT_ALLELES    } from '../../../modules/local/scatac_count_alleles/main'
include { SCATAC_CREATE_ANNDATA   } from '../../../modules/local/scatac_create_anndata/main'
include { SCATAC_PSEUDOBULK       } from '../../../modules/local/scatac_pseudobulk/main'
include { WASP2_ANALYZE_IMBALANCE } from '../../../../nf-modules/modules/wasp2/analyze_imbalance/main'

workflow WASP_ALLELIC_SC {
    take:
    ch_fragments   // channel: [ val(meta), path(fragments.tsv.gz), path(fragments.tbi), path(barcodes), path(peaks) ]
    ch_vcf         // channel: [ val(meta), path(vcf), path(tbi) ]

    main:
    ch_versions = Channel.empty()

    // Convert VCF to BED for heterozygous SNP positions
    WASP2_VCF_TO_BED ( ch_vcf, '' )
    ch_versions = ch_versions.mix(WASP2_VCF_TO_BED.out.versions)

    // Prepare input for counting: combine fragments with SNP BED
    // Extract barcodes and peaks from the fragment channel
    ch_count_input = ch_fragments
        .combine(WASP2_VCF_TO_BED.out.bed)
        .map { meta, fragments, fragments_tbi, barcodes, peaks, var_meta, bed ->
            [ meta, fragments, fragments_tbi, bed, barcodes, peaks ]
        }

    // Count fragment overlaps per cell at SNP positions
    // Supports optional barcode filtering and peak region filtering
    SCATAC_COUNT_ALLELES (
        ch_count_input.map { meta, fragments, tbi, bed, barcodes, peaks -> [ meta, fragments, tbi, bed ] },
        ch_count_input.map { meta, fragments, tbi, bed, barcodes, peaks -> barcodes },
        ch_count_input.map { meta, fragments, tbi, bed, barcodes, peaks -> peaks }
    )
    ch_versions = ch_versions.mix(SCATAC_COUNT_ALLELES.out.versions.first())

    // Create AnnData H5AD output for scverse ecosystem (unless skipped)
    ch_anndata = Channel.empty()
    ch_cell_qc = Channel.empty()
    ch_zarr = Channel.empty()

    if (!params.skip_anndata) {
        SCATAC_CREATE_ANNDATA (
            SCATAC_COUNT_ALLELES.out.counts,
            params.create_zarr ?: false
        )
        ch_versions = ch_versions.mix(SCATAC_CREATE_ANNDATA.out.versions.first())
        ch_anndata = SCATAC_CREATE_ANNDATA.out.anndata
        ch_cell_qc = SCATAC_CREATE_ANNDATA.out.cell_qc
        ch_zarr = SCATAC_CREATE_ANNDATA.out.zarr
    }

    // Pseudo-bulk aggregation and statistical analysis (unless skipped)
    ch_pseudobulk = Channel.empty()
    ch_imbalance = Channel.empty()

    if (!params.skip_pseudobulk) {
        // Aggregate per-cell counts to pseudo-bulk
        SCATAC_PSEUDOBULK ( SCATAC_COUNT_ALLELES.out.counts )
        ch_versions = ch_versions.mix(SCATAC_PSEUDOBULK.out.versions.first())
        ch_pseudobulk = SCATAC_PSEUDOBULK.out.counts

        // Analyze allelic imbalance on pseudo-bulk counts
        WASP2_ANALYZE_IMBALANCE ( SCATAC_PSEUDOBULK.out.counts )
        ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions.first())
        ch_imbalance = WASP2_ANALYZE_IMBALANCE.out.results
    }

    emit:
    cell_counts = SCATAC_COUNT_ALLELES.out.counts   // channel: [ val(meta), path(counts.tsv) ]
    count_stats = SCATAC_COUNT_ALLELES.out.stats    // channel: [ val(meta), path(stats.tsv) ]
    anndata     = ch_anndata                        // channel: [ val(meta), path(*.h5ad) ]
    zarr        = ch_zarr                           // channel: [ val(meta), path(*.zarr) ]
    cell_qc     = ch_cell_qc                        // channel: [ val(meta), path(cell_qc.tsv) ]
    pseudobulk  = ch_pseudobulk                     // channel: [ val(meta), path(pseudobulk.tsv) ]
    imbalance   = ch_imbalance                      // channel: [ val(meta), path(results.tsv) ]
    versions    = ch_versions                       // channel: [ path(versions.yml) ]
}
