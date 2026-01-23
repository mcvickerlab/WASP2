/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP_MAPPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Performs WASP2 mapping bias correction:
    1. Generate swapped-allele reads for remapping
    2. Remap reads with aligner
    3. Filter reads that don't map to same position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASP2_MAKE_READS      } from '../../modules/local/wasp2_make_reads'
include { WASP2_FILTER_REMAPPED } from '../../modules/local/wasp2_filter_remapped'
include { BWA_MEM               } from '../../modules/nf-core/bwa/mem/main'
include { BOWTIE2_ALIGN         } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/samtools/index/main'

workflow WASP_MAPPING {
    take:
    ch_bam           // channel: [ val(meta), path(bam), path(bai) ]
    ch_vcf           // channel: path(vcf)
    ch_aligner_index // channel: path(index)
    ch_fasta         // channel: path(fasta)
    aligner          // string: 'bwa' or 'bowtie2'

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Generate reads with swapped alleles for remapping
    //
    WASP2_MAKE_READS(
        ch_bam,
        ch_vcf
    )
    ch_versions = ch_versions.mix(WASP2_MAKE_READS.out.versions.first())

    //
    // Prepare FASTQ channel for remapping
    // Transform from [meta, fq1, fq2] to [meta_remap, [fq1, fq2]]
    //
    ch_remap_reads = WASP2_MAKE_READS.out.fastq
        .map { meta, fq1, fq2 ->
            def meta_remap = meta.clone()
            meta_remap.id = "${meta.id}_remap"
            meta_remap.single_end = false
            [ meta_remap, [ fq1, fq2 ] ]
        }

    //
    // MODULE: Remap swapped-allele reads
    //
    if (aligner == 'bwa') {
        BWA_MEM(
            ch_remap_reads,
            ch_aligner_index,
            ch_fasta,
            true  // sort_bam
        )
        ch_remapped_raw = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    } else {
        BOWTIE2_ALIGN(
            ch_remap_reads,
            ch_aligner_index,
            ch_fasta,
            false,  // save_unaligned
            true    // sort_bam
        )
        ch_remapped_raw = BOWTIE2_ALIGN.out.aligned
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    }

    //
    // MODULE: Index remapped BAM (aligners already sort when sort_bam=true)
    //
    SAMTOOLS_INDEX(ch_remapped_raw)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Combine BAM with index
    ch_remapped = ch_remapped_raw
        .join(SAMTOOLS_INDEX.out.bai, by: [0], failOnMismatch: true)

    //
    // Join remapped BAM with WASP intermediate files for filtering
    //
    ch_wasp_intermediates = WASP2_MAKE_READS.out.to_remap_bam
        .join(WASP2_MAKE_READS.out.keep_bam, by: [0], failOnMismatch: true)
        .join(WASP2_MAKE_READS.out.json, by: [0], failOnMismatch: true)
        .map { meta, to_remap, keep, json -> [ meta.id, meta, to_remap, keep, json ] }

    ch_filter_input = ch_remapped
        .map { meta, bam, bai -> [ meta.id.replace('_remap', ''), bam, bai ] }
        .join(ch_wasp_intermediates, by: [0], failOnMismatch: true)
        .map { _id, bam, bai, meta, to_remap, keep, json -> [ meta, bam, bai, to_remap, keep, json ] }

    WASP2_FILTER_REMAPPED(
        ch_filter_input  // Already structured as [meta, bam, bai, to_remap, keep, json]
    )
    ch_versions = ch_versions.mix(WASP2_FILTER_REMAPPED.out.versions.first())

    emit:
    bam      = WASP2_FILTER_REMAPPED.out.bam      // channel: [ val(meta), path(bam), path(bai) ]
    stats    = WASP2_FILTER_REMAPPED.out.stats    // channel: [ val(meta), path(stats) ]
    versions = ch_versions                         // channel: path(versions.yml)
}
