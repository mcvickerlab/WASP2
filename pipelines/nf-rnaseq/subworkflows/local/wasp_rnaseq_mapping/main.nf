/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP_RNASEQ_MAPPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Performs WASP2 mapping bias correction for RNA-seq data using STAR aligner:
    1. Generate swapped-allele reads for remapping
    2. Remap reads with STAR
    3. Filter reads that don't map to same position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_ALIGN as STAR_ALIGN_REMAP } from '../../../modules/local/star_align/main'
include { WASP2_UNIFIED_MAKE_READS       } from '../../../modules/local/wasp2_unified_make_reads/main'
include { WASP2_FILTER_REMAPPED          } from '../../../modules/local/wasp2_filter_remapped/main'

workflow WASP_RNASEQ_MAPPING {
    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ]
    ch_vcf        // channel: [ val(meta), path(vcf), path(vcf_index) ]
    ch_star_index // path: STAR index directory
    ch_gtf        // path: GTF annotation file (optional)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Generate reads with swapped alleles for remapping
    // The unified pipeline handles VCF processing internally
    //
    WASP2_UNIFIED_MAKE_READS(
        ch_bam,
        ch_vcf.first()
    )
    ch_versions = ch_versions.mix(WASP2_UNIFIED_MAKE_READS.out.versions)

    //
    // Prepare FASTQ channel for remapping
    //
    ch_remap_reads = WASP2_UNIFIED_MAKE_READS.out.remap_fastq
        .map { meta, r1, r2 ->
            def meta_remap = meta.clone()
            meta_remap.id = "${meta.id}_remap"
            meta_remap.single_end = false
            tuple(meta_remap, [r1, r2])
        }

    //
    // MODULE: Remap swapped-allele reads with STAR
    //
    STAR_ALIGN_REMAP(
        ch_remap_reads,
        ch_star_index,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_REMAP.out.versions)

    //
    // Join remapped BAM with WASP intermediate files for filtering
    //
    ch_filter_input = STAR_ALIGN_REMAP.out.bam
        .map { meta, bam, bai -> [ meta.id.replace('_remap', ''), meta, bam, bai ] }
        .join(
            WASP2_UNIFIED_MAKE_READS.out.to_remap_bam
                .map { meta, bam -> [ meta.id, bam ] },
            by: [0]
        )
        .join(
            WASP2_UNIFIED_MAKE_READS.out.keep_bam
                .map { meta, bam -> [ meta.id, bam ] },
            by: [0]
        )
        .join(
            WASP2_UNIFIED_MAKE_READS.out.wasp_json
                .map { meta, json -> [ meta.id, json ] },
            by: [0]
        )
        .map { _id, meta, bam, bai, to_remap, keep, json ->
            // Restore original meta (without _remap suffix)
            def meta_orig = meta.clone()
            meta_orig.id = _id
            [ meta_orig, bam, bai, to_remap, keep, json ]
        }
        .multiMap { meta, remap_bam, remap_bai, to_remap, keep, json ->
            remapped: tuple(meta, remap_bam, remap_bai)
            to_remap: tuple(meta, to_remap)
            keep:     tuple(meta, keep)
            json:     tuple(meta, json)
        }

    //
    // MODULE: Filter remapped reads using WASP algorithm
    //
    WASP2_FILTER_REMAPPED(
        ch_filter_input.remapped,
        ch_filter_input.to_remap,
        ch_filter_input.keep,
        ch_filter_input.json
    )
    ch_versions = ch_versions.mix(WASP2_FILTER_REMAPPED.out.versions)

    emit:
    bam      = WASP2_FILTER_REMAPPED.out.bam      // channel: [ val(meta), path(bam), path(bai) ]
    stats    = WASP2_FILTER_REMAPPED.out.stats    // channel: [ val(meta), path(stats) ]
    versions = ch_versions                         // channel: path(versions.yml)
}
