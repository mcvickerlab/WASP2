//
// Prepare genome reference files and aligner indices
//

include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD  } from '../../../modules/nf-core/bowtie2/index/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    main:
    ch_versions = Channel.empty()

    //
    // Validate required parameters
    //
    if (!params.fasta) {
        error "ERROR: --fasta is required. Please provide a reference FASTA file."
    }

    //
    // Load FASTA reference with metadata
    //
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id: fasta.baseName], fasta] }

    //
    // Generate or load FASTA index
    //
    ch_fasta_fai = Channel.empty()
    if (params.fasta_fai) {
        ch_fasta_fai = Channel.fromPath(params.fasta_fai, checkIfExists: true)
            .map { fai -> [[id: file(params.fasta).baseName], fai] }
    } else {
        SAMTOOLS_FAIDX ( ch_fasta )
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    //
    // Prepare BWA index (only when aligner is 'bwa')
    //
    ch_bwa_index = Channel.empty()
    if (params.aligner == 'bwa') {
        if (params.bwa_index) {
            ch_bwa_index = Channel.fromPath(params.bwa_index, checkIfExists: true).collect()
        } else {
            BWA_INDEX ( ch_fasta )
            ch_bwa_index = BWA_INDEX.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    //
    // Prepare Bowtie2 index (only when aligner is 'bowtie2')
    //
    ch_bowtie2_index = Channel.empty()
    if (params.aligner == 'bowtie2') {
        if (params.bowtie2_index) {
            ch_bowtie2_index = Channel.fromPath(params.bowtie2_index, checkIfExists: true).collect()
        } else {
            BOWTIE2_BUILD ( ch_fasta )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta.map { meta, fasta -> fasta }.collect()  // channel: path(fasta)
    fasta_fai     = ch_fasta_fai.map { meta, fai -> fai }.collect()  // channel: path(fasta.fai)
    bwa_index     = ch_bwa_index                                      // channel: path(bwa_index)
    bowtie2_index = ch_bowtie2_index                                  // channel: path(bowtie2_index)

    versions      = ch_versions                                       // channel: path(versions.yml)
}
