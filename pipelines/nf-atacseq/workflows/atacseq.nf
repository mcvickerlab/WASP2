/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MACS2_CALLPEAK         } from '../modules/nf-core/macs2/callpeak/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

// Local WASP2 modules
include { WASP2_COUNT_VARIANTS   } from '../modules/local/wasp2_count_variants/main'
include { WASP2_FIND_IMBALANCE   } from '../modules/local/wasp2_find_imbalance/main'

// nf-core subworkflows (standardized alignment interfaces)
include { FASTQ_ALIGN_BWA             } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { FASTQ_ALIGN_BOWTIE2         } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { BAM_MARKDUPLICATES_PICARD   } from '../subworkflows/nf-core/bam_markduplicates_picard/main'

// Local subworkflows
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome/main'
include { WASP_MAPPING           } from '../subworkflows/local/wasp_mapping/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ATACSEQ {

    take:
    ch_samplesheet // channel: [ val(meta), [ fastq_1, fastq_2 ] ]

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Validate aligner parameter
    //
    def valid_aligners = ['bwa', 'bowtie2']
    if (!valid_aligners.contains(params.aligner)) {
        error "Invalid aligner '${params.aligner}'. Must be one of: ${valid_aligners.join(', ')}"
    }

    //
    // SUBWORKFLOW: Prepare genome reference and indices
    //
    PREPARE_GENOME ()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_fasta = PREPARE_GENOME.out.fasta
    ch_vcf   = params.vcf ? Channel.fromPath(params.vcf, checkIfExists: true).collect() : Channel.empty()

    //
    // MODULE: FastQC - Raw read QC
    //
    if (!params.skip_fastqc) {
        FASTQC ( ch_samplesheet )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    }

    //
    // MODULE: Fastp - Adapter trimming and QC
    //
    if (!params.skip_trimming) {
        FASTP (
            ch_samplesheet,
            [],     // adapter_fasta
            false,  // save_trimmed_fail
            false   // save_merged
        )
        ch_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
    } else {
        ch_reads = ch_samplesheet
    }

    //
    // SUBWORKFLOW: Alignment (BWA-MEM or Bowtie2)
    // Uses standardized nf-core subworkflow interface
    //
    ch_aligned_bam    = Channel.empty()
    ch_aligned_bai    = Channel.empty()
    ch_align_stats    = Channel.empty()
    ch_align_flagstat = Channel.empty()
    ch_align_idxstats = Channel.empty()

    if (params.aligner == 'bwa') {
        FASTQ_ALIGN_BWA (
            ch_reads,
            PREPARE_GENOME.out.bwa_index,
            ch_fasta
        )
        ch_aligned_bam    = FASTQ_ALIGN_BWA.out.bam
        ch_aligned_bai    = FASTQ_ALIGN_BWA.out.bai
        ch_align_stats    = FASTQ_ALIGN_BWA.out.stats
        ch_align_flagstat = FASTQ_ALIGN_BWA.out.flagstat
        ch_align_idxstats = FASTQ_ALIGN_BWA.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions)
    } else if (params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            ch_reads,
            PREPARE_GENOME.out.bowtie2_index,
            ch_fasta
        )
        ch_aligned_bam    = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_aligned_bai    = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_align_stats    = FASTQ_ALIGN_BOWTIE2.out.stats
        ch_align_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
        ch_align_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)
    }

    // Add alignment stats to MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(ch_align_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_align_flagstat.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_align_idxstats.collect { it[1] })

    // Combine BAM with index
    ch_bam_indexed = ch_aligned_bam
        .join(ch_aligned_bai, by: [0], failOnMismatch: true)

    //
    // SUBWORKFLOW: Mark duplicates with Picard and run BAM stats (optional)
    //
    ch_fasta_fai = PREPARE_GENOME.out.fasta_fai

    if (!params.skip_dedup) {
        BAM_MARKDUPLICATES_PICARD (
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam] },
            ch_fasta,
            ch_fasta_fai
        )
        ch_bam_dedup = BAM_MARKDUPLICATES_PICARD.out.bam
            .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], failOnMismatch: true)
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        // Add deduplication stats to MultiQC
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect { it[1] })
    } else {
        ch_bam_dedup = ch_bam_indexed
    }

    //
    // MODULE: Peak calling (MACS2)
    //
    if (!params.skip_peak_calling) {
        MACS2_CALLPEAK (
            ch_bam_dedup.map { meta, bam, bai -> [meta, bam] },
            params.macs_gsize
        )
        ch_peaks = MACS2_CALLPEAK.out.peak
        ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions.first())
    } else {
        // Use provided peaks file - validate it exists
        if (!params.peaks) {
            error "ERROR: --peaks is required when --skip_peak_calling is enabled"
        }
        ch_peaks = Channel.fromPath(params.peaks, checkIfExists: true)
            .map { peaks -> [[id: 'provided_peaks'], peaks] }
    }

    //
    // SUBWORKFLOW: WASP2 Mapping Bias Correction
    //
    // Determine which index to pass to WASP_MAPPING
    ch_aligner_index = params.aligner == 'bwa'
        ? PREPARE_GENOME.out.bwa_index
        : PREPARE_GENOME.out.bowtie2_index

    if (!params.skip_wasp && params.vcf) {
        WASP_MAPPING (
            ch_bam_dedup,
            ch_vcf,
            ch_aligner_index,
            ch_fasta,
            params.aligner
        )
        ch_wasp_bam = WASP_MAPPING.out.bam
        ch_versions = ch_versions.mix(WASP_MAPPING.out.versions)
    } else {
        ch_wasp_bam = ch_bam_dedup
    }

    //
    // MODULE: WASP2 Allele Counting at Peaks
    //
    ch_counts     = Channel.empty()
    ch_ai_results = Channel.empty()

    if (params.vcf) {
        if (params.skip_peak_calling && params.peaks) {
            // Single consensus peaks file for all samples
            ch_peaks_for_count = Channel.fromPath(params.peaks, checkIfExists: true).first()
            WASP2_COUNT_VARIANTS (
                ch_wasp_bam,
                ch_vcf,
                ch_peaks_for_count
            )
        } else {
            // Sample-specific peaks from MACS2
            ch_bam_with_peaks = ch_wasp_bam
                .join(ch_peaks, by: [0], failOnMismatch: true)
                .map { meta, bam, bai, peaks -> [ meta, bam, bai, peaks ] }

            WASP2_COUNT_VARIANTS (
                ch_bam_with_peaks.map { meta, bam, bai, peaks -> [meta, bam, bai] },
                ch_vcf,
                ch_bam_with_peaks.map { meta, bam, bai, peaks -> peaks }
            )
        }
        ch_counts = WASP2_COUNT_VARIANTS.out.counts
        ch_versions = ch_versions.mix(WASP2_COUNT_VARIANTS.out.versions.first())

        //
        // MODULE: WASP2 Allelic Imbalance Analysis
        //
        WASP2_FIND_IMBALANCE (
            ch_counts,
            params.wasp_min_count,
            params.wasp_pseudocount
        )
        ch_ai_results = WASP2_FIND_IMBALANCE.out.results
        ch_versions = ch_versions.mix(WASP2_FIND_IMBALANCE.out.versions.first())
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: false).ifEmpty([])

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            [],  // extra_multiqc_config
            []   // multiqc_logo
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    bam            = ch_wasp_bam                              // channel: [ val(meta), path(bam), path(bai) ]
    peaks          = ch_peaks                                  // channel: [ val(meta), path(peaks) ]
    counts         = ch_counts                                 // channel: [ val(meta), path(counts) ]
    ai_results     = ch_ai_results                             // channel: [ val(meta), path(results) ]
    multiqc_report = ch_multiqc_report                         // channel: path(report)
    versions       = ch_versions                               // channel: path(versions.yml)
}
