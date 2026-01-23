#!/usr/bin/env nextflow
/*
 * Integration test for WASP2_COUNT module
 * Runs actual wasp2-count on test data
 */
nextflow.enable.dsl = 2

include { WASP2_COUNT } from '../../modules/wasp2/count/main'

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

    // Run the count module
    WASP2_COUNT(bam_ch, vcf_ch, [])

    // Publish outputs
    WASP2_COUNT.out.counts
        .map { meta, counts -> counts }
        .collectFile(name: 'test.counts.tsv', storeDir: params.outdir)

    WASP2_COUNT.out.versions
        .collectFile(name: 'versions.yml', storeDir: params.outdir)
}
