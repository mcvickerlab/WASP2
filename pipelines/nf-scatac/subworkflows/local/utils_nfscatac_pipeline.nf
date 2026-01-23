/*
    Utility subworkflows for nf-scatac pipeline
*/

workflow PIPELINE_INITIALISATION {
    take:
    version
    help
    validate_params
    input

    main:
    if (version) {
        log.info "nf-scatac v${workflow.manifest.version}"
        System.exit(0)
    }

    if (help) {
        log.info helpMessage()
        System.exit(0)
    }

    ch_samplesheet = Channel.fromPath(input)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample]
            def fragments = row.fragments ? file(row.fragments) : file("${row.cellranger_dir}/outs/fragments.tsv.gz")
            def fragments_tbi = file("${fragments}.tbi")
            [ meta, fragments, fragments_tbi ]
        }

    emit:
    samplesheet = ch_samplesheet
}

workflow PIPELINE_COMPLETION {
    take:
    outdir
    multiqc_report

    main:
    log.info """
    =============================================================
     nf-scatac COMPLETE
    =============================================================
    Output: ${outdir}
    """.stripIndent()
}

def helpMessage() {
    return """
    nf-scatac - Single-Cell ATAC-seq Allelic Imbalance

    Usage:
        nextflow run nf-scatac -profile docker --input samplesheet.csv --vcf variants.vcf.gz

    Required:
        --input   Samplesheet CSV
        --vcf     VCF/BCF variant file

    See: https://github.com/Jaureguy760/WASP2-final/issues/32
    """.stripIndent()
}
