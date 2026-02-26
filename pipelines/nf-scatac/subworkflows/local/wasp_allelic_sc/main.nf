/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP_ALLELIC_SC SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP2 single-cell allelic imbalance analysis for scATAC-seq data.

    Supports two counting modes:
    1. BAM-based (allele-specific): True ref/alt counting with hap1/hap2 layers
    2. Fragment-based (overlap counting): Total coverage at SNP positions

    Includes per-cell counting, AnnData output, and optional pseudo-bulk analysis.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASP2_VCF_TO_BED              } from '../../../../nf-modules/modules/wasp2/vcf_to_bed/main'
include { WASP2_COUNT_SC                } from '../../../../nf-modules/modules/wasp2/count_sc/main'
include { SCATAC_COUNT_ALLELES          } from '../../../modules/local/scatac_count_alleles/main'
include { SCATAC_CREATE_ANNDATA         } from '../../../modules/local/scatac_create_anndata/main'
include { SCATAC_ADD_HAPLOTYPE_LAYERS   } from '../../../modules/local/scatac_add_haplotype_layers/main'
include { SCATAC_PSEUDOBULK             } from '../../../modules/local/scatac_pseudobulk/main'
include { WASP2_ANALYZE_IMBALANCE       } from '../../../../nf-modules/modules/wasp2/analyze_imbalance/main'

workflow WASP_ALLELIC_SC {
    take:
    ch_input   // channel: [ val(meta), path(fragments), path(fragments_tbi), path(barcodes), path(peaks), path(bam), path(bai) ]
    ch_vcf     // channel: [ val(meta), path(vcf), path(tbi) ]

    main:
    ch_versions = Channel.empty()
    ch_anndata = Channel.empty()
    ch_cell_qc = Channel.empty()
    ch_zarr = Channel.empty()
    ch_cell_counts = Channel.empty()
    ch_count_stats = Channel.empty()
    ch_pseudobulk = Channel.empty()
    ch_imbalance = Channel.empty()

    // Convert VCF to BED for heterozygous SNP positions (used by fragment-based counting)
    WASP2_VCF_TO_BED ( ch_vcf, '' )
    ch_versions = ch_versions.mix(WASP2_VCF_TO_BED.out.versions)

    // Branch samples by input type: BAM-based vs fragment-based
    // Note: .branch uses first-match semantics, so fragment_based is a fallback
    ch_branched = ch_input
        .branch {
            bam_based: it[5] != null && !it[5].name.startsWith('NO_FILE')  // has BAM file
            fragment_based: true                 // fallback: no BAM, use fragments
        }

    //==========================================================================
    // PATH 1: BAM-based allele-specific counting (true ref/alt/hap layers)
    //==========================================================================
    // Prepare BAM input channels for WASP2_COUNT_SC
    ch_bam_inputs = ch_branched.bam_based
        .multiMap { meta, fragments, fragments_tbi, barcodes, peaks, bam, bai ->
            main:     [ meta, bam, bai ]
            barcodes: barcodes
            peaks:    peaks
        }

    // Run allele-specific counting from BAM
    WASP2_COUNT_SC (
        ch_bam_inputs.main,
        ch_vcf,
        ch_bam_inputs.barcodes,
        ch_bam_inputs.peaks
    )
    ch_versions = ch_versions.mix(WASP2_COUNT_SC.out.versions.first().ifEmpty([]))

    // Add haplotype layers using phased VCF
    if (!params.skip_anndata) {
        SCATAC_ADD_HAPLOTYPE_LAYERS (
            WASP2_COUNT_SC.out.counts,
            ch_vcf,
            params.create_zarr ?: false
        )
        ch_versions = ch_versions.mix(SCATAC_ADD_HAPLOTYPE_LAYERS.out.versions.first().ifEmpty([]))
        ch_anndata = ch_anndata.mix(SCATAC_ADD_HAPLOTYPE_LAYERS.out.anndata)
        ch_cell_qc = ch_cell_qc.mix(SCATAC_ADD_HAPLOTYPE_LAYERS.out.cell_qc)
        ch_zarr = ch_zarr.mix(SCATAC_ADD_HAPLOTYPE_LAYERS.out.zarr)
    }

    // For BAM-based, use stats from WASP2_COUNT_SC
    ch_count_stats = ch_count_stats.mix(WASP2_COUNT_SC.out.stats.ifEmpty([]))

    //==========================================================================
    // PATH 2: Fragment-based overlap counting (total counts only)
    //==========================================================================
    // Combine fragments with SNP BED and prepare inputs for counting
    ch_frag_count_input = ch_branched.fragment_based
        .combine(WASP2_VCF_TO_BED.out.bed)
        .multiMap { meta, fragments, fragments_tbi, barcodes, peaks, bam, bai, var_meta, bed ->
            main:     [ meta, fragments, fragments_tbi, bed ]
            barcodes: barcodes
            peaks:    peaks
        }

    // Count fragment overlaps per cell at SNP positions
    SCATAC_COUNT_ALLELES (
        ch_frag_count_input.main,
        ch_frag_count_input.barcodes,
        ch_frag_count_input.peaks
    )
    ch_versions = ch_versions.mix(SCATAC_COUNT_ALLELES.out.versions.first().ifEmpty([]))
    ch_cell_counts = ch_cell_counts.mix(SCATAC_COUNT_ALLELES.out.counts)
    ch_count_stats = ch_count_stats.mix(SCATAC_COUNT_ALLELES.out.stats)

    // Create AnnData H5AD output for fragment-based counting
    if (!params.skip_anndata) {
        SCATAC_CREATE_ANNDATA (
            SCATAC_COUNT_ALLELES.out.counts,
            params.create_zarr ?: false
        )
        ch_versions = ch_versions.mix(SCATAC_CREATE_ANNDATA.out.versions.first().ifEmpty([]))
        ch_anndata = ch_anndata.mix(SCATAC_CREATE_ANNDATA.out.anndata)
        ch_cell_qc = ch_cell_qc.mix(SCATAC_CREATE_ANNDATA.out.cell_qc)
        ch_zarr = ch_zarr.mix(SCATAC_CREATE_ANNDATA.out.zarr)
    }

    //==========================================================================
    // SHARED: Pseudo-bulk aggregation and statistical analysis
    //==========================================================================
    if (!params.skip_pseudobulk) {
        // For fragment-based: use TSV counts
        SCATAC_PSEUDOBULK ( SCATAC_COUNT_ALLELES.out.counts )
        ch_versions = ch_versions.mix(SCATAC_PSEUDOBULK.out.versions.first().ifEmpty([]))
        ch_pseudobulk = ch_pseudobulk.mix(SCATAC_PSEUDOBULK.out.counts)

        WASP2_ANALYZE_IMBALANCE ( SCATAC_PSEUDOBULK.out.counts )
        ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions.first().ifEmpty([]))
        ch_imbalance = ch_imbalance.mix(WASP2_ANALYZE_IMBALANCE.out.results)

        // TODO: Add pseudo-bulk analysis for BAM-based path (aggregate from H5AD)
    }

    emit:
    cell_counts = ch_cell_counts                       // channel: [ val(meta), path(counts.tsv) ] (fragment-based only)
    count_stats = ch_count_stats                       // channel: [ val(meta), path(stats.tsv) ]
    anndata     = ch_anndata                           // channel: [ val(meta), path(*.h5ad) ]
    zarr        = ch_zarr                              // channel: [ val(meta), path(*.zarr) ]
    cell_qc     = ch_cell_qc                           // channel: [ val(meta), path(cell_qc.tsv) ]
    pseudobulk  = ch_pseudobulk                        // channel: [ val(meta), path(pseudobulk.tsv) ]
    imbalance   = ch_imbalance                         // channel: [ val(meta), path(results.tsv) ]
    versions    = ch_versions                          // channel: [ path(versions.yml) ]
}
