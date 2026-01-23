#!/usr/bin/env nextflow
/*
 * Integration test for WASP2_ANALYZE module
 * Runs actual wasp2-analyze on test data
 */
nextflow.enable.dsl = 2

include { WASP2_ANALYZE } from '../../modules/wasp2/analyze/main'

params.outdir = 'results'

workflow {
    // Define test input channels with counts file
    counts_ch = Channel.of([
        [ id:'test', phased:false ],
        file("${projectDir}/../data/sample.counts.tsv")
    ])

    // Run the analyze module
    WASP2_ANALYZE(counts_ch)

    // Publish outputs
    WASP2_ANALYZE.out.stats
        .map { meta, stats -> stats }
        .collectFile(name: 'test.stats.tsv', storeDir: params.outdir)

    WASP2_ANALYZE.out.versions
        .collectFile(name: 'versions.yml', storeDir: params.outdir)
}
