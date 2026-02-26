#!/usr/bin/env nextflow
/*
 * Benchmark integration test for WASP2 Nextflow modules
 * Uses HG00731 chr21 data from 1000 Genomes (same as test_benchmarks suite)
 * Runs: count-variants → find-imbalance pipeline via Docker
 */
nextflow.enable.dsl = 2

include { WASP2_COUNT } from '../../modules/wasp2/count/main'
include { WASP2_ANALYZE } from '../../modules/wasp2/analyze/main'

params.outdir     = 'results'
params.data_dir   = "${projectDir}/../../../../test_benchmarks/data"
params.sample     = 'HG00731'

workflow {
    // Step 1: Count variants using benchmark BAM + VCF
    bam_ch = Channel.of([
        [ id:'chr21_benchmark', sample: params.sample ],
        file("${params.data_dir}/chr21.bam"),
        file("${params.data_dir}/chr21.bam.bai")
    ])

    vcf_ch = Channel.of([
        [ id:'chr21_vcf' ],
        file("${params.data_dir}/chr21.vcf.gz"),
        file("${params.data_dir}/chr21.vcf.gz.tbi")
    ])

    WASP2_COUNT(bam_ch, vcf_ch, [])

    // Step 2: Analyze imbalance on the counts
    counts_ch = WASP2_COUNT.out.counts.map { meta, counts ->
        [ [ id: meta.id, phased: false ], counts ]
    }

    WASP2_ANALYZE(counts_ch)

    // Publish outputs for comparison
    WASP2_COUNT.out.counts
        .map { meta, counts -> counts }
        .collectFile(name: 'nf_counts.tsv', storeDir: params.outdir)

    WASP2_ANALYZE.out.stats
        .map { meta, stats -> stats }
        .collectFile(name: 'nf_ai_results.tsv', storeDir: params.outdir)
}
