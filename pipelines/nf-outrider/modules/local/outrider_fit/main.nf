process OUTRIDER_FIT {
    tag "outrider"
    label 'process_high'
    label 'process_high_memory'

    conda "bioconda::bioconductor-outrider=1.26.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-outrider:1.26.3--r44he5774e6_0' :
        'quay.io/biocontainers/bioconductor-outrider:1.26.3--r44he5774e6_0' }"

    input:
    path count_matrix
    val padj_cutoff
    val zscore_cutoff
    val encoding_dim
    val max_iterations
    val convergence
    val min_count

    output:
    path "outrider_model.rds"   , emit: model
    path "outrider_results.tsv" , emit: results
    path "outrider_summary.html", emit: summary, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def q_arg = encoding_dim ? "q_val <- ${encoding_dim}" : "ods <- estimateBestQ(ods); q_val <- getBestQ(ods)"
    """
    set -euo pipefail

    Rscript - <<'REOF'
    suppressPackageStartupMessages({
        library(OUTRIDER)
        library(BiocParallel)
        library(data.table)
    })

    message("=== OUTRIDER Aberrant Expression Detection ===")

    # Parameters from Nextflow
    counts_file    <- "${count_matrix}"
    padj_cut       <- ${padj_cutoff}
    zscore_cut     <- ${zscore_cutoff}
    max_iter       <- ${max_iterations}
    conv_thresh    <- ${convergence}
    n_threads      <- ${task.cpus}

    # Set parallel processing
    if (n_threads > 1) {
        register(MulticoreParam(n_threads))
    } else {
        register(SerialParam())
    }

    # Load count matrix
    message(sprintf("Loading: %s", counts_file))
    counts <- fread(counts_file, data.table = FALSE)
    rownames(counts) <- counts[[1]]
    counts <- counts[, -1, drop = FALSE]
    message(sprintf("Loaded: %d genes x %d samples", nrow(counts), ncol(counts)))

    # Filter low-expressed genes
    min_count_thresh <- ${min_count}
    min_samples <- max(2, floor(ncol(counts) * 0.5))
    row_sums <- rowSums(counts >= min_count_thresh)
    keep_genes <- row_sums >= min_samples
    counts_filtered <- counts[keep_genes, , drop = FALSE]
    message(sprintf("After filtering: %d genes", nrow(counts_filtered)))

    if (nrow(counts_filtered) < 10) {
        message(sprintf("ERROR: Too few genes (%d) for OUTRIDER. Minimum 10 required.", nrow(counts_filtered)))
        quit(status = 1)
    }

    # Create OutriderDataSet
    ods <- OutriderDataSet(countData = as.matrix(counts_filtered))
    ods <- filterExpression(ods, minCounts = TRUE, filterGenes = FALSE)
    ods <- estimateSizeFactors(ods)

    # Encoding dimension
    ${q_arg}
    message(sprintf("Encoding dimension q: %d", q_val))

    # Fit OUTRIDER model
    message("Fitting OUTRIDER autoencoder...")
    ods <- OUTRIDER(ods, q = q_val, controlData = TRUE,
                    iterations = max_iter, convergence = conv_thresh)

    saveRDS(ods, "outrider_model.rds")

    # Extract results
    res <- results(ods, padjCutoff = padj_cut, zScoreCutoff = zscore_cut, all = TRUE)
    res_dt <- as.data.table(res)
    res_dt[, aberrant := (padjust < padj_cut) & (abs(zScore) > zscore_cut)]
    fwrite(res_dt, "outrider_results.tsv", sep = "\\t")

    n_outliers <- sum(res_dt[["aberrant"]], na.rm = TRUE)
    message(sprintf("Aberrant expression events: %d / %d tested", n_outliers, nrow(res_dt)))
    message("=== OUTRIDER Complete ===")
REOF

    # Validate output was created
    if [ ! -f "outrider_results.tsv" ]; then
        echo "ERROR: OUTRIDER results file was not created" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //' | cut -d' ' -f1)
        OUTRIDER: \$(Rscript -e "cat(as.character(packageVersion('OUTRIDER')))")
    END_VERSIONS
    """

    stub:
    """
    touch outrider_model.rds
    cat <<-END_HEADER > outrider_results.tsv
    geneID	sampleID	pValue	padjust	zScore	l2fc	rawcounts	normcounts	meanCorrected	theta	aberrant
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.4.3
        OUTRIDER: 1.26.3
    END_VERSIONS
    """
}
