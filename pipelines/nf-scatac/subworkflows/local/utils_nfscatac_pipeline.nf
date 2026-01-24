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
    // Supports optional barcodes and peaks columns for filtering
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

            // Resolve fragments file path
            def fragments = row.fragments
                ? file(row.fragments, checkIfExists: true)
                : file("${row.cellranger_dir}/outs/fragments.tsv.gz", checkIfExists: true)
            def fragments_tbi = file("${fragments}.tbi", checkIfExists: true)

            // Optional: cell barcode whitelist file
            def barcodes = row.barcodes && row.barcodes.trim()
                ? file(row.barcodes, checkIfExists: true)
                : file('NO_FILE')

            // Optional: peak BED file for restricting analysis to peak regions
            def peaks = row.peaks && row.peaks.trim()
                ? file(row.peaks, checkIfExists: true)
                : file('NO_FILE')

            [ meta, fragments, fragments_tbi, barcodes, peaks ]
        }

    emit:
    samplesheet = ch_samplesheet  // channel: [ val(meta), path(fragments), path(fragments_tbi), path(barcodes), path(peaks) ]
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
      - anndata/: AnnData H5AD files for scverse ecosystem
      - cell_qc/: Cell QC metrics
      - imbalance/: Allelic imbalance analysis results (pseudo-bulk)
    """.stripIndent()
}

def helpMessage() {
    return """
    nf-scatac - Single-Cell ATAC-seq Allelic Imbalance Pipeline

    Usage:
        nextflow run nf-scatac -profile docker --input samplesheet.csv --vcf variants.vcf.gz

    Required:
        --input       Samplesheet CSV (see format below)
        --vcf         Indexed VCF/BCF with heterozygous SNPs

    Optional:
        --outdir                  Output directory [default: ./results]
        --min_fragments_per_cell  Minimum fragments per cell to include [default: 1000]
        --min_cells_per_snp       Minimum cells per SNP for pseudo-bulk [default: 3]
        --create_zarr             Also output Zarr format for GenVarLoader [default: false]
        --skip_pseudobulk         Skip pseudo-bulk aggregation [default: false]
        --skip_anndata            Skip AnnData H5AD creation [default: false]

    Samplesheet format:
        sample,fragments,cellranger_dir,barcode_tag,chemistry,barcodes,peaks
        sample1,/path/to/fragments.tsv.gz,,CB,10x-atac-v2,/path/to/barcodes.txt,/path/to/peaks.bed
        sample2,,/path/to/cellranger/output,CB,10x-atac-v2,,

    Column descriptions:
        sample         - Sample identifier (required)
        fragments      - Path to fragments.tsv.gz (required if no cellranger_dir)
        cellranger_dir - Path to CellRanger ATAC output (required if no fragments)
        barcode_tag    - BAM tag for cell barcodes [default: CB]
        chemistry      - Library chemistry [default: 10x-atac-v2]
        barcodes       - Optional: file with valid cell barcodes (one per line)
        peaks          - Optional: BED file with peak regions to restrict analysis

    Profiles:
        -profile docker        Run with Docker
        -profile singularity   Run with Singularity
        -profile conda         Run with Conda
        -profile test_stub     Run stub tests
        -profile test_real     Run integration tests with real data

    Output:
        results/
        ├── allele_counts/      # Per-cell allele counts TSV
        ├── anndata/            # AnnData H5AD files
        ├── zarr/               # Zarr files (if --create_zarr)
        ├── cell_qc/            # Cell QC metrics
        ├── pseudobulk/         # Pseudo-bulk aggregated counts
        ├── imbalance/          # Allelic imbalance results
        └── pipeline_info/      # Execution reports
    """.stripIndent()
}
