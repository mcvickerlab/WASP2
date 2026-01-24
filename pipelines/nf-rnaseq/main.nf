#!/usr/bin/env nextflow
/*
========================================================================================
    WASP2 RNA-seq ASE Pipeline
========================================================================================
    Github : https://github.com/Jaureguy760/WASP2-exp
    Docs   : https://wasp2.readthedocs.io

    Pipeline for RNA-seq Allele-Specific Expression (ASE) analysis using WASP2
    for mapping bias correction.

    Workflow:
    FASTQ -> STAR align -> WASP2 make-reads -> STAR remap -> WASP2 filter -> count -> analyze
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// Print parameter summary
log.info """
==============================================
 WASP2 RNA-seq ASE Pipeline v${workflow.manifest.version}
==============================================
 input         : ${params.input}
 vcf           : ${params.vcf}
 star_index    : ${params.star_index}
 gtf           : ${params.gtf ?: 'not provided'}
 outdir        : ${params.outdir}
 skip_analysis : ${params.skip_analysis}
----------------------------------------------
"""

// Validate required parameters
if (!params.input)      { exit 1, "ERROR: --input samplesheet not specified" }
if (!params.vcf)        { exit 1, "ERROR: --vcf variants file not specified" }
if (!params.star_index) { exit 1, "ERROR: --star_index STAR index not specified" }

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { STAR_ALIGN as STAR_ALIGN_INITIAL } from '../nf-modules/modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_REMAP   } from '../nf-modules/modules/star/align/main'
include { WASP2_UNIFIED_MAKE_READS         } from '../nf-modules/modules/wasp2/unified_make_reads/main'
include { WASP2_FILTER_REMAPPED            } from '../nf-modules/modules/wasp2/filter_remapped/main'
include { WASP2_COUNT_ALLELES              } from '../nf-modules/modules/wasp2/count_alleles/main'
include { WASP2_ANALYZE_IMBALANCE          } from '../nf-modules/modules/wasp2/analyze_imbalance/main'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow RNASEQ_ASE {

    // Channel for version tracking
    ch_versions = Channel.empty()

    // Initialize output channel for conditional analysis step
    ch_ai_results = Channel.empty()

    //
    // Parse input samplesheet with validation
    //
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample)  { exit 1, "ERROR: Samplesheet missing 'sample' column" }
            if (!row.fastq_1) { exit 1, "ERROR: Samplesheet missing 'fastq_1' for sample: ${row.sample}" }

            def meta = [id: row.sample, single_end: !row.fastq_2, sample: row.sample]
            def fastqs = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
            tuple(meta, fastqs)
        }
        .set { ch_reads }

    //
    // Prepare VCF channel with index (singleton reference file)
    //
    ch_vcf = Channel.fromPath(params.vcf)
        .map { vcf ->
            def vcf_index = file("${vcf}.tbi")
            if (!vcf_index.exists()) {
                vcf_index = file("${vcf}.csi")
            }
            if (!vcf_index.exists()) {
                exit 1, "ERROR: VCF index not found. Expected ${vcf}.tbi or ${vcf}.csi"
            }
            tuple([id: 'reference'], vcf, vcf_index)
        }

    //
    // Load reference files
    //
    ch_star_index = file(params.star_index)
    ch_gtf = params.gtf ? file(params.gtf) : []

    //
    // STEP 1: Initial STAR alignment
    //
    STAR_ALIGN_INITIAL(
        ch_reads,
        ch_star_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_INITIAL.out.versions)

    //
    // STEP 2: Generate swapped allele FASTQs (Rust unified pipeline)
    // The unified pipeline handles VCF processing internally
    //
    WASP2_UNIFIED_MAKE_READS(
        STAR_ALIGN_INITIAL.out.bam,
        ch_vcf.first()
    )
    ch_versions = ch_versions.mix(WASP2_UNIFIED_MAKE_READS.out.versions)

    //
    // STEP 3: Remap swapped allele reads
    //
    ch_remap_reads = WASP2_UNIFIED_MAKE_READS.out.remap_fastq
        .map { meta, r1, r2 ->
            tuple(meta, [r1, r2])
        }

    STAR_ALIGN_REMAP(
        ch_remap_reads,
        ch_star_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_REMAP.out.versions)

    //
    // STEP 4: Filter remapped reads using WASP algorithm
    // Join channels by meta.id to ensure sample synchronization
    //
    ch_filter_input = STAR_ALIGN_REMAP.out.bam
        .join(WASP2_UNIFIED_MAKE_READS.out.to_remap_bam, by: [0])
        .join(WASP2_UNIFIED_MAKE_READS.out.keep_bam, by: [0])
        .join(WASP2_UNIFIED_MAKE_READS.out.wasp_json, by: [0])
        .multiMap { meta, remap_bam, remap_bai, to_remap, keep, json ->
            remapped: tuple(meta, remap_bam, remap_bai)
            to_remap: tuple(meta, to_remap)
            keep:     tuple(meta, keep)
            json:     tuple(meta, json)
        }

    WASP2_FILTER_REMAPPED(
        ch_filter_input.remapped,
        ch_filter_input.to_remap,
        ch_filter_input.keep,
        ch_filter_input.json
    )
    ch_versions = ch_versions.mix(WASP2_FILTER_REMAPPED.out.versions)

    //
    // STEP 5: Count alleles at heterozygous SNPs
    //
    WASP2_COUNT_ALLELES(
        WASP2_FILTER_REMAPPED.out.bam,
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
    // Collect and deduplicate versions
    //
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    wasp_bam = WASP2_FILTER_REMAPPED.out.bam
    counts   = WASP2_COUNT_ALLELES.out.counts
    results  = ch_ai_results
    versions = ch_versions
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    RNASEQ_ASE()
}

/*
========================================================================================
    COMPLETION SUMMARY
========================================================================================
*/

workflow.onComplete {
    log.info """
    ==============================================
     WASP2 RNA-seq ASE Pipeline Complete!
    ==============================================
     Started      : ${workflow.start}
     Completed    : ${workflow.complete}
     Duration     : ${workflow.duration}
     Success      : ${workflow.success}
     Work Dir     : ${workflow.workDir}
     Output Dir   : ${params.outdir}
    ==============================================
    """
}
