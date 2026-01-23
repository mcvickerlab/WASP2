#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-scatac: Single-Cell ATAC-seq Allelic Imbalance Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Jaureguy760/WASP2-final
    Issue  : #32
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { SCATAC                   } from './workflows/scatac'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfscatac_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfscatac_pipeline'

workflow NFSCATAC {
    take:
    samplesheet

    main:
    SCATAC ( samplesheet )

    emit:
    allele_counts = SCATAC.out.allele_counts
    pseudobulk    = SCATAC.out.pseudobulk
    versions      = SCATAC.out.versions
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
        Channel.empty()  // MultiQC not implemented
    )
}
