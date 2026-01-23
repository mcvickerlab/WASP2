//
// Alignment with BWA-MEM
//

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS      } from '../bam_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {
    take:
    ch_reads   // channel: [ val(meta), path(reads) ]
    ch_index   // channel: path(index)
    ch_fasta   // channel: path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Align reads with BWA-MEM (outputs sorted BAM)
    //
    BWA_MEM (
        ch_reads,
        ch_index,
        ch_fasta,
        true  // sort_bam
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // Index BAM file
    //
    SAMTOOLS_INDEX ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Join BAM and BAI
    //
    ch_bam_bai = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], failOnMismatch: true)

    //
    // Run BAM stats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BWA_MEM.out.bam                    // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), path(bai) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat    // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats    // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                        // channel: path(versions.yml)
}
