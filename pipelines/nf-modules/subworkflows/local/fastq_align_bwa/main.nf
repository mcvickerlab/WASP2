//
// FASTQ_ALIGN_BWA: Align FASTQ reads using BWA MEM with sorting and statistics
//
// Standardized nf-core subworkflow pattern for BWA alignment
//

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {
    take:
    reads     // channel: [ val(meta), path(reads) ]
    index     // channel: [ val(meta), path(index) ]
    fasta     // channel: [ val(meta), path(fasta) ]
    sort_bam  // value: boolean - whether to sort output BAM

    main:
    ch_versions = Channel.empty()

    // Align reads with BWA MEM
    BWA_MEM (
        reads,
        index,
        fasta,
        sort_bam
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Sort BAM and generate statistics
    BAM_SORT_STATS_SAMTOOLS (
        BWA_MEM.out.bam,
        fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam       // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai       // channel: [ val(meta), path(bai) ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi       // channel: [ val(meta), path(csi) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats     // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats  // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                            // channel: [ versions.yml ]
}
