#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-scatac: Single-Cell ATAC-seq Allelic Imbalance Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Allelic imbalance analysis for single-cell ATAC-seq data using WASP2.

    Input modes:
    - BAM input: True allele-specific counting with ref/alt/hap1/hap2 layers
    - Fragments input: Overlap counting (total counts only)

    Features:
    - 10x Genomics CellRanger ATAC output support
    - Cell barcode filtering
    - Peak region filtering
    - Per-cell allele counts with AnnData/H5AD output
    - AnnData layers: X (total), ref, alt, hap1, hap2 (when BAM provided)
    - Zarr output for GenVarLoader integration
    - Pseudo-bulk aggregation for statistical power
    - Integration with ArchR/Signac via scverse ecosystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { SCATAC                   } from './workflows/scatac'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfscatac_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfscatac_pipeline'

workflow NFSCATAC {
    take:
    samplesheet  // channel: [ val(meta), path(fragments), path(fragments_tbi), path(barcodes), path(peaks), path(bam), path(bai) ]

    main:
    SCATAC ( samplesheet )

    emit:
    cell_counts = SCATAC.out.cell_counts  // channel: [ val(meta), path(counts.tsv) ]
    count_stats = SCATAC.out.count_stats  // channel: [ val(meta), path(stats.tsv) ]
    anndata     = SCATAC.out.anndata      // channel: [ val(meta), path(*.h5ad) ]
    zarr        = SCATAC.out.zarr         // channel: [ val(meta), path(*.zarr) ]
    cell_qc     = SCATAC.out.cell_qc      // channel: [ val(meta), path(cell_qc.tsv) ]
    pseudobulk  = SCATAC.out.pseudobulk   // channel: [ val(meta), path(pseudobulk.tsv) ]
    imbalance   = SCATAC.out.imbalance    // channel: [ val(meta), path(results.tsv) ]
    versions    = SCATAC.out.versions     // channel: [ path(versions.yml) ]
}

workflow {
    main:
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.input
    )

    NFSCATAC ( PIPELINE_INITIALISATION.out.samplesheet )

    PIPELINE_COMPLETION (
        params.outdir,
        Channel.empty()
    )
}
