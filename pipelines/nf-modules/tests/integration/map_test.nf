#!/usr/bin/env nextflow
/*
 * Integration test for WASP2 MAP modules
 * Tests the full WASP mapping bias correction workflow
 */
nextflow.enable.dsl = 2

include { WASP2_MAP_MAKE_READS } from '../../modules/wasp2/map/main'
include { WASP2_MAP_FILTER } from '../../modules/wasp2/map/main'

params.outdir = 'results'

workflow {
    // Define test input channels
    bam_ch = Channel.of([
        [ id:'test', sample:'sample1' ],
        file("${projectDir}/../data/minimal.bam"),
        file("${projectDir}/../data/minimal.bam.bai")
    ])

    vcf_ch = Channel.of([
        [ id:'test_vcf' ],
        file("${projectDir}/../data/sample.vcf.gz"),
        file("${projectDir}/../data/sample.vcf.gz.tbi")
    ])

    // Step 1: Generate swapped allele reads
    WASP2_MAP_MAKE_READS(bam_ch, vcf_ch)

    // For this test, we simulate remapping by using the original BAM
    // In a real workflow, users would run their aligner here
    remapped_ch = Channel.of([
        [ id:'test' ],
        file("${projectDir}/../data/minimal_remap.bam"),
        file("${projectDir}/../data/minimal_remap.bam.bai")
    ])

    // Step 2: Filter remapped reads
    WASP2_MAP_FILTER(
        remapped_ch,
        WASP2_MAP_MAKE_READS.out.to_remap_bam,
        WASP2_MAP_MAKE_READS.out.keep_bam,
        WASP2_MAP_MAKE_READS.out.wasp_json
    )

    // Publish outputs
    WASP2_MAP_FILTER.out.bam
        .map { meta, bam, bai -> bam }
        .collectFile(name: 'test.wasp_filt.bam', storeDir: params.outdir)
}
