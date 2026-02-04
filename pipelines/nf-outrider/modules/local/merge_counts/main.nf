process MERGE_COUNTS {
    tag "merge_counts"
    label 'process_medium'

    conda "${moduleDir}/../../../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wasp2:1.2.1--pyhdfd78af_0' :
        'biocontainers/wasp2:1.2.1--pyhdfd78af_0' }"

    input:
    path gene_counts  // Collection of gene count files

    output:
    path "count_matrix.tsv", emit: count_matrix
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Merge individual sample gene counts into a matrix for OUTRIDER
    python3 << 'EOF'
import pandas as pd
import glob
import sys

def main():
    # Find all gene count files
    count_files = glob.glob('*.gene_counts.tsv')
    print(f"Found {len(count_files)} count files")

    if len(count_files) == 0:
        print("ERROR: No gene count files found in working directory", file=sys.stderr)
        sys.exit(1)

    # Load and merge all samples with error handling
    all_counts = []
    for f in count_files:
        try:
            df = pd.read_csv(f, sep='\\t')
            if df.empty:
                print(f"WARNING: Empty count file: {f}", file=sys.stderr)
                continue
            # Validate required columns
            required_cols = ['gene_id', 'sample_id', 'total_count']
            missing_cols = [c for c in required_cols if c not in df.columns]
            if missing_cols:
                print(f"ERROR: File {f} missing required columns: {missing_cols}", file=sys.stderr)
                print(f"Found columns: {list(df.columns)}", file=sys.stderr)
                sys.exit(1)
            all_counts.append(df)
        except pd.errors.EmptyDataError:
            print(f"WARNING: Empty file skipped: {f}", file=sys.stderr)
            continue
        except pd.errors.ParserError as e:
            print(f"ERROR: Failed to parse {f}: {e}", file=sys.stderr)
            sys.exit(1)
        except FileNotFoundError:
            print(f"ERROR: Count file not found: {f}", file=sys.stderr)
            sys.exit(1)

    if len(all_counts) == 0:
        print("ERROR: No valid count files to merge", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_counts, ignore_index=True)

    # Pivot to create gene x sample matrix using total counts
    # OUTRIDER expects genes in rows, samples in columns
    try:
        count_matrix = combined.pivot_table(
            index='gene_id',
            columns='sample_id',
            values='total_count',
            aggfunc='sum'
        ).fillna(0).astype(int)
    except KeyError as e:
        print(f"ERROR: Missing column during pivot: {e}", file=sys.stderr)
        sys.exit(1)
    except (TypeError, ValueError) as e:
        print(f"ERROR: Non-numeric data in total_count column: {e}", file=sys.stderr)
        sys.exit(1)
    except OverflowError as e:
        print(f"ERROR: Count values overflow integer range: {e}", file=sys.stderr)
        sys.exit(1)

    if count_matrix.empty:
        print("ERROR: Count matrix is empty after pivot", file=sys.stderr)
        sys.exit(1)

    # Save main count matrix for OUTRIDER
    count_matrix.to_csv('count_matrix.tsv', sep='\\t')

    print(f"Created count matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")

if __name__ == '__main__':
    main()
EOF

    # Validate output was created
    if [ ! -f "count_matrix.tsv" ]; then
        echo "ERROR: Count matrix was not created" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_HEADER > count_matrix.tsv
    gene_id	sample1	sample2	sample3
    ENSG00000000001	100	200	150
    END_HEADER

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pandas: 2.0.0
    END_VERSIONS
    """
}
