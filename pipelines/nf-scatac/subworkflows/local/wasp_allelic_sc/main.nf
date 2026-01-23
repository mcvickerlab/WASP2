/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP_ALLELIC_SC SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WASP2 single-cell allelic imbalance analysis for scATAC-seq data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASP2_VCF_TO_BED        } from '../../../../nf-modules/modules/wasp2/vcf_to_bed/main'
include { SCATAC_COUNT_ALLELES    } from '../../../modules/local/scatac_count_alleles/main'
include { WASP2_ANALYZE_IMBALANCE } from '../../../../nf-modules/modules/wasp2/analyze_imbalance/main'

workflow WASP_ALLELIC_SC {
    take:
    ch_fragments   // channel: [ val(meta), path(fragments.tsv.gz), path(fragments.tbi) ]
    ch_vcf         // channel: [ val(meta), path(vcf), path(tbi) ]

    main:
    ch_versions = Channel.empty()

    // Convert VCF to BED for heterozygous SNP positions
    WASP2_VCF_TO_BED ( ch_vcf, '' )
    ch_versions = ch_versions.mix(WASP2_VCF_TO_BED.out.versions)

    // Combine fragments with SNP BED (keep sample meta, discard variant meta)
    ch_count_input = ch_fragments
        .combine(WASP2_VCF_TO_BED.out.bed)
        .map { meta, fragments, fragments_tbi, var_meta, bed ->
            [ meta, fragments, fragments_tbi, bed ]
        }

    // Count fragment overlaps per cell at SNP positions
    SCATAC_COUNT_ALLELES ( ch_count_input )
    ch_versions = ch_versions.mix(SCATAC_COUNT_ALLELES.out.versions.first())

    // Analyze allelic imbalance (pseudo-bulk aggregation for statistical power)
    WASP2_ANALYZE_IMBALANCE ( SCATAC_COUNT_ALLELES.out.counts )
    ch_versions = ch_versions.mix(WASP2_ANALYZE_IMBALANCE.out.versions.first())

    emit:
    cell_counts = SCATAC_COUNT_ALLELES.out.counts       // channel: [ val(meta), path(counts.tsv) ]
    imbalance   = WASP2_ANALYZE_IMBALANCE.out.results   // channel: [ val(meta), path(results.tsv) ]
    versions    = ch_versions                           // channel: [ path(versions.yml) ]
}
