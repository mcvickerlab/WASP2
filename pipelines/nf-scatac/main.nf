#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-scatac: Single-Cell ATAC-seq Allelic Imbalance Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { SCATAC                   } from './workflows/scatac'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfscatac_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfscatac_pipeline'

workflow NFSCATAC {
    take:
    samplesheet  // channel: [ val(meta), path(fragments), path(fragments_tbi) ]

    main:
    SCATAC ( samplesheet )

    emit:
    allele_counts = SCATAC.out.allele_counts  // channel: [ val(meta), path(counts.tsv) ]
    imbalance     = SCATAC.out.imbalance      // channel: [ val(meta), path(results.tsv) ]
    versions      = SCATAC.out.versions       // channel: [ path(versions.yml) ]
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
