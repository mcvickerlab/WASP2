process MAE_DETECT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)
    val min_count
    val padj_cutoff
    val alt_ratio_threshold

    output:
    tuple val(meta), path("*.mae_results.tsv"), emit: mae_results
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    set -euo pipefail

    # Detect mono-allelic expression using binomial test
    python3 << 'EOF'
import pandas as pd
import sys
from scipy.stats import binomtest, false_discovery_control

# Output columns for consistent schema
OUTPUT_COLUMNS = [
    'chrom', 'pos', 'gene_id', 'ref_count', 'alt_count',
    'total_count', 'alt_ratio', 'pvalue', 'mae_status',
    'padj', 'significant', 'sample_id'
]

def binomial_pvalue(ref_count, alt_count, expected_ratio=0.5):
    \"\"\"Calculate p-value for allelic imbalance using binomial test.\"\"\"
    total = ref_count + alt_count
    if total == 0:
        return 1.0
    result = binomtest(int(alt_count), int(total), expected_ratio)
    return result.pvalue

def detect_mae(counts_df, min_count, alt_ratio_threshold):
    \"\"\"Detect mono-allelic expression sites.\"\"\"
    results = []

    for _, row in counts_df.iterrows():
        ref = row['ref_count']
        alt = row['alt_count']
        total = ref + alt

        # Skip zero total (division by zero guard)
        if total == 0:
            continue

        if total < min_count:
            continue

        alt_ratio = alt / total
        ref_ratio = ref / total

        # Calculate p-value for imbalance
        pval = binomial_pvalue(ref, alt)

        # Determine MAE status
        is_mae_alt = alt_ratio >= alt_ratio_threshold
        is_mae_ref = ref_ratio >= alt_ratio_threshold
        mae_status = 'MAE_ALT' if is_mae_alt else ('MAE_REF' if is_mae_ref else 'BIALLELIC')

        results.append({
            'chrom': row.get('chrom', 'NA'),
            'pos': row.get('pos', 0),
            'gene_id': row.get('gene_id', 'NA'),
            'ref_count': ref,
            'alt_count': alt,
            'total_count': total,
            'alt_ratio': alt_ratio,
            'pvalue': pval,
            'mae_status': mae_status
        })

    return pd.DataFrame(results)

def main():
    # Load counts with error handling
    try:
        counts_df = pd.read_csv('${counts}', sep='\\t')
    except pd.errors.EmptyDataError:
        print(f"ERROR: Input file '${counts}' is empty", file=sys.stderr)
        sys.exit(1)
    except pd.errors.ParserError as e:
        print(f"ERROR: Failed to parse '${counts}': {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"ERROR: Input file '${counts}' not found", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required_cols = ['ref_count', 'alt_count']
    missing_cols = [c for c in required_cols if c not in counts_df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
        print(f"Found columns: {list(counts_df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Detect MAE
    mae_df = detect_mae(counts_df, ${min_count}, ${alt_ratio_threshold})

    # Handle empty results with consistent schema
    if len(mae_df) == 0:
        print("WARNING: No variants passed filtering. Check min_count threshold.", file=sys.stderr)
        mae_df = pd.DataFrame(columns=OUTPUT_COLUMNS[:-1])  # Exclude sample_id, added below
    else:
        # Apply multiple testing correction
        mae_df['padj'] = false_discovery_control(mae_df['pvalue'].values, method='bh')
        mae_df['significant'] = (mae_df['padj'] < ${padj_cutoff}) & (mae_df['mae_status'] != 'BIALLELIC')

    # Add sample info
    mae_df['sample_id'] = '${meta.id}'

    # Save results
    mae_df.to_csv('${prefix}.mae_results.tsv', sep='\\t', index=False)

    # Summary statistics
    n_mae = int(mae_df['significant'].sum()) if len(mae_df) > 0 else 0
    print(f"Found {n_mae} significant MAE sites out of {len(mae_df)} tested")

if __name__ == '__main__':
    main()
EOF

    # Validate output was created
    if [ ! -f "${prefix}.mae_results.tsv" ]; then
        echo "ERROR: Output file ${prefix}.mae_results.tsv was not created" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    cat <<-END_HEADER > ${prefix}.mae_results.tsv
    chrom	pos	gene_id	ref_count	alt_count	total_count	alt_ratio	pvalue	mae_status	padj	significant	sample_id
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        scipy: 1.12.0
    END_VERSIONS
    """
}
