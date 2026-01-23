/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility subworkflows for nf-scatac pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    // Parse samplesheet with scATAC-specific meta map fields
    ch_samplesheet = Channel.fromPath(input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) {
                error "Samplesheet error: 'sample' column is required. Found columns: ${row.keySet()}"
            }
            if (!row.fragments && !row.cellranger_dir) {
                error "Samplesheet error for '${row.sample}': provide 'fragments' or 'cellranger_dir'"
            }

            def meta = [
                id: row.sample,
                single_end: false,
                cell_barcode_tag: row.barcode_tag ?: 'CB',
                chemistry: row.chemistry ?: '10x-atac-v2'
            ]

            def fragments = row.fragments
                ? file(row.fragments, checkIfExists: true)
                : file("${row.cellranger_dir}/outs/fragments.tsv.gz", checkIfExists: true)
            def fragments_tbi = file("${fragments}.tbi", checkIfExists: true)

            [ meta, fragments, fragments_tbi ]
        }

    emit:
    samplesheet = ch_samplesheet  // channel: [ val(meta), path(fragments), path(fragments_tbi) ]
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

    Results:
      - allele_counts/: Per-cell allele counts at het SNPs
      - imbalance/: Allelic imbalance analysis results
    """.stripIndent()
}

def helpMessage() {
    return """
    nf-scatac - Single-Cell ATAC-seq Allelic Imbalance

    Usage:
        nextflow run nf-scatac -profile docker --input samplesheet.csv --vcf variants.vcf.gz

    Required:
        --input       Samplesheet CSV (columns: sample, fragments or cellranger_dir)
        --vcf         Indexed VCF/BCF with heterozygous SNPs

    Optional:
        --outdir                  Output directory [default: ./results]
        --min_fragments_per_cell  Minimum fragments per cell [default: 1000]

    Samplesheet:
        sample,fragments,cellranger_dir,barcode_tag,chemistry
        sample1,/path/to/fragments.tsv.gz,,,CB,10x-atac-v2
        sample2,,/path/to/cellranger/output,,CB,10x-atac-v2

    Profiles:
        -profile docker        Run with Docker
        -profile singularity   Run with Singularity
        -profile conda         Run with Conda
        -profile test_stub     Run stub test
    """.stripIndent()
}
