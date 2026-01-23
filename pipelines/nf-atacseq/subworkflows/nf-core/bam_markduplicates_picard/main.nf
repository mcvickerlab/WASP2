//
// Mark duplicates with Picard and run BAM stats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam    // channel: [ val(meta), path(bam) ]
    ch_fasta  // channel: path(fasta)
    ch_fai    // channel: path(fasta_fai)

    main:
    ch_versions = Channel.empty()

    //
    // Mark duplicates with Picard
    //
    PICARD_MARKDUPLICATES (
        ch_bam,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //
    // Join BAM and BAI for stats
    //
    ch_bam_bai = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai, by: [0], failOnMismatch: true)

    //
    // Run BAM stats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = PICARD_MARKDUPLICATES.out.bai      // channel: [ val(meta), path(bai) ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics  // channel: [ val(meta), path(metrics) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat    // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats    // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                        // channel: path(versions.yml)
}
