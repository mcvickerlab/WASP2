//
// Sort, index BAM file and run samtools stats, flagstat, and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam    // channel: [ val(meta), path(bam) ]
    ch_fasta  // channel: path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Sort BAM file
    //
    SAMTOOLS_SORT ( ch_bam, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // Index sorted BAM file
    //
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Join BAM and BAI for stats - failOnMismatch ensures every BAM has a corresponding index
    //
    ch_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], failOnMismatch: true)

    //
    // Run samtools stats, flagstat, and idxstats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = SAMTOOLS_SORT.out.bam              // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), path(bai) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat    // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats    // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                        // channel: path(versions.yml)
}
