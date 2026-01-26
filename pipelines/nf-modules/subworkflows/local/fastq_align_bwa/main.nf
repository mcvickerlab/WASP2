//
// Alignment with BWA-MEM and BAM statistics
//
// Shared FASTQ_ALIGN_BWA subworkflow following nf-core patterns.
//

include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {

    take:
    ch_reads   // channel: [ val(meta), path(reads) ]
    ch_index   // channel: [ val(meta), path(index) ]
    ch_fasta   // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // Extract paths from index/fasta channels for BWA_MEM
    // Using .first() ensures proper value channel semantics for multi-sample runs
    //
    ch_index_path = ch_index.map { meta, idx -> idx }.first()
    ch_fasta_path = ch_fasta.map { meta, fa -> fa }.first()

    //
    // Align reads with BWA-MEM
    //
    BWA_MEM(
        ch_reads,
        ch_index_path,
        ch_fasta_path,
        false  // sort_bam - let BAM_SORT_STATS_SAMTOOLS handle sorting
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // Sort BAM, index, and collect statistics
    //
    BAM_SORT_STATS_SAMTOOLS(
        BWA_MEM.out.bam,
        ch_fasta_path
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam       // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai       // channel: [ val(meta), path(bai) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats     // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats  // channel: [ val(meta), path(idxstats) ]
    versions = ch_versions                           // channel: [ path(versions.yml) ]
}
