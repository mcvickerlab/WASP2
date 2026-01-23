/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { MACS2_CALLPEAK         } from '../modules/nf-core/macs2/callpeak/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

// Local WASP2 modules
include { WASP2_MAKE_READS       } from '../modules/local/wasp2_make_reads'
include { WASP2_FILTER_REMAPPED  } from '../modules/local/wasp2_filter_remapped'
include { WASP2_COUNT_VARIANTS   } from '../modules/local/wasp2_count_variants'
include { WASP2_FIND_IMBALANCE   } from '../modules/local/wasp2_find_imbalance'

// Subworkflows
include { WASP_MAPPING           } from '../subworkflows/local/wasp_mapping'

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
    ch_counts        = Channel.empty()
    ch_ai_results    = Channel.empty()

    //
    // Prepare reference files
    //
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta, checkIfExists: true).collect() : Channel.empty()
    ch_vcf   = params.vcf   ? Channel.fromPath(params.vcf, checkIfExists: true).collect()   : Channel.empty()

    // Validate aligner parameter
    def valid_aligners = ['bwa', 'bowtie2']
    if (!valid_aligners.contains(params.aligner)) {
        error "Invalid aligner '${params.aligner}'. Must be one of: ${valid_aligners.join(', ')}"
    }

    // Prepare aligner index
    if (params.aligner == 'bwa') {
        if (params.bwa_index) {
            ch_aligner_index = Channel.fromPath(params.bwa_index, checkIfExists: true).collect()
        } else {
            BWA_INDEX(ch_fasta.map { fasta -> [[id: 'genome'], fasta] })
            ch_aligner_index = BWA_INDEX.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    } else if (params.aligner == 'bowtie2') {
        if (params.bowtie2_index) {
            ch_aligner_index = Channel.fromPath(params.bowtie2_index, checkIfExists: true).collect()
        } else {
            BOWTIE2_BUILD(ch_fasta.map { fasta -> [[id: 'genome'], fasta] })
            ch_aligner_index = BOWTIE2_BUILD.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    //
    // MODULE: FastQC - Raw read QC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_samplesheet)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    }

    //
    // MODULE: Fastp - Adapter trimming and QC
    //
    if (!params.skip_trimming) {
        FASTP(
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
    // MODULE: Alignment (BWA-MEM or Bowtie2)
    //
    if (params.aligner == 'bwa') {
        BWA_MEM(
            ch_reads,
            ch_aligner_index,
            ch_fasta,
            true  // sort_bam
        )
        ch_aligned = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    } else {
        BOWTIE2_ALIGN(
            ch_reads,
            ch_aligner_index,
            ch_fasta,
            false,  // save_unaligned
            true    // sort_bam
        )
        ch_aligned = BOWTIE2_ALIGN.out.aligned
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    }

    //
    // MODULE: Index BAM (aligners already sort when sort_bam=true)
    //
    SAMTOOLS_INDEX(ch_aligned)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Combine BAM with index
    ch_bam_indexed = ch_aligned
        .join(SAMTOOLS_INDEX.out.bai, by: [0])

    //
    // MODULE: Mark duplicates (optional)
    //
    if (!params.skip_dedup) {
        PICARD_MARKDUPLICATES(
            ch_bam_indexed.map { meta, bam, bai -> [meta, bam] },
            ch_fasta,
            []  // fasta_fai
        )
        ch_bam_dedup = PICARD_MARKDUPLICATES.out.bam
            .join(PICARD_MARKDUPLICATES.out.bai, by: [0])
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect { it[1] })
    } else {
        ch_bam_dedup = ch_bam_indexed
    }

    //
    // MODULE: BAM statistics
    //
    SAMTOOLS_STATS(ch_bam_dedup, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { it[1] })

    SAMTOOLS_FLAGSTAT(ch_bam_dedup)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect { it[1] })

    //
    // MODULE: Peak calling (MACS2)
    //
    if (!params.skip_peak_calling) {
        MACS2_CALLPEAK(
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
    if (!params.skip_wasp && params.vcf) {
        WASP_MAPPING(
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
    if (params.vcf) {
        // Determine peaks to use for counting
        // Option 1: User-provided consensus peaks file (single file for all samples)
        // Option 2: Sample-specific peaks from MACS2 (joined by sample ID)
        if (params.skip_peak_calling && params.peaks) {
            // Single consensus peaks file for all samples
            ch_peaks_for_count = Channel.fromPath(params.peaks, checkIfExists: true).first()
            WASP2_COUNT_VARIANTS(
                ch_wasp_bam,
                ch_vcf,
                ch_peaks_for_count
            )
        } else {
            // Sample-specific peaks from MACS2 - join BAM with peaks by sample meta
            ch_bam_with_peaks = ch_wasp_bam
                .join(ch_peaks, by: [0])
                .map { meta, bam, bai, peaks ->
                    [ meta, bam, bai, peaks ]
                }

            WASP2_COUNT_VARIANTS(
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
        WASP2_FIND_IMBALANCE(
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
    if (!params.skip_multiqc) {
        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: false).ifEmpty([])

        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            [],  // extra_multiqc_config
            []   // multiqc_logo
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    } else {
        ch_multiqc_report = Channel.empty()
    }

    emit:
    bam            = ch_wasp_bam                              // channel: [ val(meta), path(bam), path(bai) ]
    peaks          = ch_peaks                                  // channel: [ val(meta), path(peaks) ]
    counts         = params.vcf ? ch_counts : Channel.empty()  // channel: [ val(meta), path(counts) ]
    ai_results     = params.vcf ? ch_ai_results : Channel.empty() // channel: [ val(meta), path(results) ]
    multiqc_report = ch_multiqc_report                         // channel: path(report)
    versions       = ch_versions                               // channel: path(versions.yml)
}
