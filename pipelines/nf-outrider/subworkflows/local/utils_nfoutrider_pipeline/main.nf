/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE UTILITY SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Common utility functions for nf-outrider pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version        // boolean: show version
    help           // boolean: show help
    input          // string: path to samplesheet

    main:
    //
    // Print help message if requested
    //
    if (help) {
        log.info helpMessage()
        System.exit(0)
    }

    //
    // Print version if requested
    //
    if (version) {
        log.info "nf-outrider version ${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Validate inputs
    //
    if (!input) {
        error "ERROR: --input samplesheet is required"
    }

    //
    // Parse samplesheet
    //
    ch_samplesheet = Channel
        .fromPath(input, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            validateSamplesheetRow(row)
        }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    outdir          // string: output directory
    multiqc_report  // channel: multiqc report

    main:
    // Completion message
    workflow.onComplete {
        if (workflow.success) {
            log.info "Pipeline completed successfully!"
            log.info "Results are available in: ${outdir}"
        } else {
            log.error "Pipeline completed with errors"
        }
    }

    // Error handling
    workflow.onError {
        log.error "Pipeline execution stopped with an error"
        log.error "Error message: ${workflow.errorMessage}"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateSamplesheetRow(row) {
    /**
     * Validate a samplesheet row and return a channel entry
     *
     * Expected columns:
     *   sample  - Sample ID (required)
     *   bam     - Path to BAM file (required)
     *   bai     - Path to BAM index (optional, auto-detected if missing)
     */

    // Check required columns
    if (!row.sample) {
        error "ERROR: 'sample' column is required in samplesheet"
    }
    if (!row.bam) {
        error "ERROR: 'bam' column is required in samplesheet for sample: ${row.sample}"
    }

    // Create meta map
    def meta = [:]
    meta.id = row.sample
    meta.sample = row.sample

    // Resolve paths: handle ${projectDir} in samplesheets and relative paths
    def resolve_path = { p ->
        p = p.replace('${projectDir}', projectDir.toString())
        def f = file(p)
        if (!f.exists() && !p.startsWith('/')) { f = file("${projectDir}/${p}") }
        return f
    }

    // Validate BAM file exists
    def bam = resolve_path(row.bam)
    if (!bam.exists()) {
        error "ERROR: BAM file not found for ${row.sample}: ${row.bam}"
    }

    // Find BAM index
    def bai = null
    if (row.bai && row.bai.trim() != '') {
        bai = resolve_path(row.bai)
        if (!bai.exists()) {
            error "ERROR: BAI file not found for ${row.sample}: ${row.bai}"
        }
    } else {
        // Auto-detect BAI
        def bai_path = "${bam}.bai"
        def alt_bai_path = bam.toString().replaceAll(/\.bam$/, '.bai')
        if (file(bai_path).exists()) {
            bai = file(bai_path)
        } else if (file(alt_bai_path).exists()) {
            bai = file(alt_bai_path)
        } else {
            error "ERROR: BAM index not found for ${row.bam}. " +
                  "Please provide .bai file or specify 'bai' column in samplesheet."
        }
    }

    return [meta, bam, bai]
}

def helpMessage() {
    """
    =========================================
     nf-outrider v${workflow.manifest.version}
    =========================================
    WASP2 + OUTRIDER for Aberrant Expression Detection

    Usage:
      nextflow run nf-outrider -profile <profile> --input <samplesheet.csv> --vcf <variants.vcf.gz> --gtf <annotation.gtf>

    Required:
      --input         Path to samplesheet CSV with columns: sample, bam, [bai]
      --vcf           Path to VCF file with heterozygous variants
      --gtf           Path to gene annotation GTF

    Optional:
      --outdir        Output directory [default: ./results]
      --skip_mae      Skip mono-allelic expression analysis

    OUTRIDER options:
      --outrider_padj         Adjusted p-value cutoff [default: 0.05]
      --outrider_zScore       Z-score cutoff [default: 2]
      --outrider_q            Encoding dimension [default: auto]

    Profiles:
      -profile docker         Run with Docker
      -profile singularity    Run with Singularity
      -profile conda          Run with Conda
      -profile test           Run with minimal test data

    For more information, see: https://github.com/your-org/WASP2
    """.stripIndent()
}
