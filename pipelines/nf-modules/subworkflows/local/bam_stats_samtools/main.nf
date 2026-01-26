//
// BAM_STATS_SAMTOOLS: Generate BAM statistics using samtools
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    fasta    // channel: [ val(meta), path(fasta) ] (optional, for CRAM)

    main:
    ch_versions = Channel.empty()

    // Run samtools stats
    SAMTOOLS_STATS ( bam_bai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    // Run samtools flagstat
    SAMTOOLS_FLAGSTAT ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // Run samtools idxstats
    SAMTOOLS_IDXSTATS ( bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
