//
// Alignment with Bowtie2
//

include { BOWTIE2_ALIGN           } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS      } from '../bam_stats_samtools/main'

workflow FASTQ_ALIGN_BOWTIE2 {
    take:
    ch_reads   // channel: [ val(meta), path(reads) ]
    ch_index   // channel: path(index)
    ch_fasta   // channel: path(fasta)

    main:
    ch_versions = Channel.empty()

    //
    // Align reads with Bowtie2 (outputs sorted BAM)
    //
    BOWTIE2_ALIGN (
        ch_reads,
        ch_index,
        ch_fasta,
        false,  // save_unaligned
        true    // sort_bam
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Index BAM file
    //
    SAMTOOLS_INDEX ( BOWTIE2_ALIGN.out.aligned )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Join BAM and BAI
    //
    ch_bam_bai = BOWTIE2_ALIGN.out.aligned
        .join(SAMTOOLS_INDEX.out.bai, by: [0], failOnMismatch: true)

    //
    // Run BAM stats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BOWTIE2_ALIGN.out.aligned          // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), path(bai) ]
    log_out  = BOWTIE2_ALIGN.out.log              // channel: [ val(meta), path(log) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat    // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats    // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                        // channel: path(versions.yml)
}
