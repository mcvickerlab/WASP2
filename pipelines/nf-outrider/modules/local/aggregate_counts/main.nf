process AGGREGATE_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jaureguy760/wasp2:latest' :
        'ghcr.io/jaureguy760/wasp2:latest' }"

    input:
    tuple val(meta), path(counts)
    path gtf
    val method  // 'sum', 'mean', or 'max'

    output:
    tuple val(meta), path("*.gene_counts.tsv"), emit: gene_counts
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    set -euo pipefail

    # Aggregate variant-level allele counts to gene level
    python3 << 'EOF'
import pandas as pd
import sys

def parse_gtf_gene_map(gtf_path):
    \"\"\"Parse GTF to map variant positions to genes.\"\"\"
    gene_map = {}  # chrom:pos -> gene_id
    chromosomes_in_gtf = set()
    try:
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if len(fields) < 9:
                    continue
                if fields[2] != 'gene':
                    continue
                chrom = fields[0]
                chromosomes_in_gtf.add(chrom)
                start = int(fields[3])
                end = int(fields[4])
                # Parse attributes
                attrs = {}
                for attr in fields[8].rstrip(';').split(';'):
                    attr = attr.strip()
                    if ' ' in attr:
                        key, val = attr.split(' ', 1)
                        attrs[key] = val.strip('"')
                gene_id = attrs.get('gene_id', 'unknown')
                # Store gene boundaries
                for pos in range(start, end + 1, 1000):  # Sample positions
                    key = f"{chrom}:{pos // 1000}"
                    if key not in gene_map:
                        gene_map[key] = []
                    gene_map[key].append((start, end, gene_id))
    except FileNotFoundError:
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"ERROR: Permission denied reading GTF: {gtf_path}", file=sys.stderr)
        sys.exit(1)
    except UnicodeDecodeError as e:
        print(f"ERROR: GTF file encoding error (expected UTF-8): {e}", file=sys.stderr)
        sys.exit(1)
    except (ValueError, IndexError) as e:
        print(f"ERROR: Malformed GTF line: {e}", file=sys.stderr)
        sys.exit(1)
    return gene_map, chromosomes_in_gtf

def find_gene(chrom, pos, gene_map):
    \"\"\"Find gene containing a variant position.\"\"\"
    key = f"{chrom}:{pos // 1000}"
    if key in gene_map:
        for start, end, gene_id in gene_map[key]:
            if start <= pos <= end:
                return gene_id
    return None

def main():
    # Load counts with error handling
    try:
        counts_df = pd.read_csv('${counts}', sep='\\t')
    except pd.errors.EmptyDataError:
        print(f"ERROR: Input counts file is empty: ${counts}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.ParserError as e:
        print(f"ERROR: Failed to parse counts file: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"ERROR: Input counts file not found: ${counts}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required_cols = ['chrom', 'pos', 'ref_count', 'alt_count']
    missing_cols = [c for c in required_cols if c not in counts_df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
        print(f"Found columns: {list(counts_df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Parse GTF and map variants to genes
    gene_map, gtf_chroms = parse_gtf_gene_map('${gtf}')

    if len(gene_map) == 0:
        print("ERROR: No genes found in GTF file", file=sys.stderr)
        sys.exit(1)

    # Check for chromosome naming convention mismatch
    variant_chroms = set(counts_df['chrom'].unique())
    common_chroms = variant_chroms & gtf_chroms
    if len(common_chroms) == 0 and len(variant_chroms) > 0:
        print(f"ERROR: Chromosome naming mismatch - no overlap detected!", file=sys.stderr)
        print(f"  Variants use: {sorted(list(variant_chroms))[:5]}...", file=sys.stderr)
        print(f"  GTF uses: {sorted(list(gtf_chroms))[:5]}...", file=sys.stderr)
        print("  Solution: Ensure consistent chromosome naming convention.", file=sys.stderr)
        print("  Common fix: Add or remove 'chr' prefix from VCF or GTF.", file=sys.stderr)
        sys.exit(1)

    # Add gene column
    counts_df['gene_id'] = counts_df.apply(
        lambda row: find_gene(row['chrom'], row['pos'], gene_map),
        axis=1
    )

    # Filter variants without gene assignment
    n_total = len(counts_df)
    counts_df = counts_df[counts_df['gene_id'].notna()]
    n_mapped = len(counts_df)

    if n_mapped == 0:
        print(f"ERROR: No variants mapped to genes (0/{n_total})", file=sys.stderr)
        print("Check chromosome naming conventions between VCF and GTF", file=sys.stderr)
        sys.exit(1)

    mapping_rate = n_mapped / n_total * 100 if n_total > 0 else 0
    if mapping_rate < 50:
        print(f"WARNING: Low gene mapping rate: {n_mapped}/{n_total} ({mapping_rate:.1f}%)", file=sys.stderr)
        print("  Possible causes:", file=sys.stderr)
        print("  - Variants outside annotated gene regions (intronic/intergenic)", file=sys.stderr)
        print("  - GTF annotation version mismatch with genome build", file=sys.stderr)
        print("  - VCF contains non-genic variants (e.g., regulatory regions)", file=sys.stderr)

    # Aggregate by gene using specified method
    agg_method = '${method}'
    agg_cols = {'ref_count': agg_method, 'alt_count': agg_method}
    if 'other_count' in counts_df.columns:
        agg_cols['other_count'] = agg_method
    gene_counts = counts_df.groupby('gene_id').agg(agg_cols).reset_index()

    # Calculate total count and allelic ratio (with division by zero protection)
    gene_counts['total_count'] = gene_counts['ref_count'] + gene_counts['alt_count']
    gene_counts['alt_ratio'] = (gene_counts['alt_count'] / gene_counts['total_count']).fillna(0.0)

    # Add sample info
    gene_counts['sample_id'] = '${meta.id}'

    # Save
    gene_counts.to_csv('${prefix}.gene_counts.tsv', sep='\\t', index=False)
    print(f"Aggregated {n_mapped} variants to {len(gene_counts)} genes ({mapping_rate:.1f}% mapped)")

if __name__ == '__main__':
    main()
EOF

    # Validate output was created
    if [ ! -f "${prefix}.gene_counts.tsv" ]; then
        echo "ERROR: Output file ${prefix}.gene_counts.tsv was not created" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = (task.ext.prefix ?: "${meta.id}").replaceAll(/[^a-zA-Z0-9._-]/, '_')
    """
    cat <<-END_HEADER > ${prefix}.gene_counts.tsv
    gene_id	ref_count	alt_count	total_count	alt_ratio	sample_id
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
    END_VERSIONS
    """
}
