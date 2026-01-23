//
// Prepare genome reference files and aligner indices
//

include { BWA_INDEX     } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/index/main'

workflow PREPARE_GENOME {
    take:
    aligner  // val: 'bwa' or 'bowtie2'

    main:
    ch_versions = Channel.empty()

    //
    // Load FASTA reference
    //
    ch_fasta = params.fasta
        ? Channel.fromPath(params.fasta, checkIfExists: true).collect()
        : Channel.empty()

    //
    // Prepare BWA index
    //
    ch_bwa_index = Channel.empty()
    if (aligner == 'bwa') {
        if (params.bwa_index) {
            ch_bwa_index = Channel.fromPath(params.bwa_index, checkIfExists: true).collect()
        } else {
            BWA_INDEX ( ch_fasta.map { fasta -> [[id: 'genome'], fasta] } )
            ch_bwa_index = BWA_INDEX.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    //
    // Prepare Bowtie2 index
    //
    ch_bowtie2_index = Channel.empty()
    if (aligner == 'bowtie2') {
        if (params.bowtie2_index) {
            ch_bowtie2_index = Channel.fromPath(params.bowtie2_index, checkIfExists: true).collect()
        } else {
            BOWTIE2_BUILD ( ch_fasta.map { fasta -> [[id: 'genome'], fasta] } )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index.map { meta, index -> index }
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta          // channel: path(fasta)
    bwa_index     = ch_bwa_index      // channel: path(bwa_index)
    bowtie2_index = ch_bowtie2_index  // channel: path(bowtie2_index)

    versions      = ch_versions       // channel: [ path(versions.yml) ]
}
