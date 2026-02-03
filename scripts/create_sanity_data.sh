#!/usr/bin/env bash
#
# Create chr21 sanity test data from HG00731 RNA-seq data.
#
# This script extracts chr21 subset from the existing HG00731 benchmark data
# and generates expected outputs by running the WASP2 pipeline.
#
# Usage:
#   ./scripts/create_sanity_data.sh [OUTPUT_DIR]
#
# Requirements:
#   - samtools
#   - bcftools (with tabix)
#   - WASP2 environment activated
#
# Source data:
#   /iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/
#   WASP2_current/cvpc/WASP2-exp/paper/figure2/data/hg00731/

set -euo pipefail

# Configuration
SRC_DIR="/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/paper/figure2/data/hg00731"
OUT_DIR="${1:-/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/benchmarking/sanity_test}"
VERSION="v1"
CHROMOSOME="chr21"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."

    for cmd in samtools bcftools tabix bgzip; do
        if ! command -v "$cmd" &> /dev/null; then
            log_error "$cmd is required but not installed."
            exit 1
        fi
    done

    # Check source files exist
    if [[ ! -f "$SRC_DIR/original.bam" ]]; then
        log_error "Source BAM not found: $SRC_DIR/original.bam"
        exit 1
    fi

    if [[ ! -f "$SRC_DIR/HG00731_het_only_chr.vcf.v4.2.gz" ]]; then
        log_error "Source VCF not found: $SRC_DIR/HG00731_het_only_chr.vcf.v4.2.gz"
        exit 1
    fi

    log_info "All prerequisites satisfied."
}

# Create output directory
create_output_dir() {
    log_info "Creating output directory: $OUT_DIR"
    mkdir -p "$OUT_DIR"
}

# Extract chr21 BAM subset
extract_chr21_bam() {
    log_info "Extracting $CHROMOSOME reads from BAM..."
    local out_bam="$OUT_DIR/chr21.bam"

    samtools view -b "$SRC_DIR/original.bam" "$CHROMOSOME" > "$out_bam"
    samtools index "$out_bam"

    local read_count
    read_count=$(samtools view -c "$out_bam")
    log_info "Extracted $read_count reads to $out_bam"
}

# Extract chr21 VCF subset
extract_chr21_vcf() {
    log_info "Extracting $CHROMOSOME variants from VCF..."
    local out_vcf="$OUT_DIR/chr21.vcf.gz"

    tabix -h "$SRC_DIR/HG00731_het_only_chr.vcf.v4.2.gz" "$CHROMOSOME" | bgzip > "$out_vcf"
    tabix -p vcf "$out_vcf"

    local variant_count
    variant_count=$(bcftools view -H "$out_vcf" | wc -l)
    log_info "Extracted $variant_count variants to $out_vcf"
}

# Generate expected outputs using WASP2 pipeline
generate_expected_outputs() {
    log_info "Generating expected outputs using WASP2 pipeline..."

    local bam="$OUT_DIR/chr21.bam"
    local vcf="$OUT_DIR/chr21.vcf.gz"
    local temp_dir="$OUT_DIR/temp"
    mkdir -p "$temp_dir"

    # Generate expected counts
    log_info "Running counting pipeline..."
    WASP2 counting count-variants \
        "$bam" \
        "$vcf" \
        --out_file "$OUT_DIR/expected_counts.tsv" \
        --temp_loc "$temp_dir" \
        2>&1 | tee "$OUT_DIR/counting.log" || {
        log_warn "Counting pipeline may have warnings, check log"
    }

    # Generate expected remapping outputs
    log_info "Running mapping pipeline..."
    WASP2 mapping make-reads \
        "$bam" \
        "$vcf" \
        --out_dir "$OUT_DIR/wasp_output" \
        --temp_loc "$temp_dir" \
        2>&1 | tee "$OUT_DIR/mapping.log" || {
        log_warn "Mapping pipeline may have warnings, check log"
    }

    # Move FASTQ outputs to expected locations
    if [[ -d "$OUT_DIR/wasp_output" ]]; then
        # Find and copy R1/R2 FASTQ files
        find "$OUT_DIR/wasp_output" -name "*_R1*.fq*" -o -name "*_1.fq*" | head -1 | while read -r f; do
            if [[ -n "$f" ]]; then
                if [[ "$f" == *.gz ]]; then
                    cp "$f" "$OUT_DIR/expected_r1.fq.gz"
                else
                    gzip -c "$f" > "$OUT_DIR/expected_r1.fq.gz"
                fi
            fi
        done

        find "$OUT_DIR/wasp_output" -name "*_R2*.fq*" -o -name "*_2.fq*" | head -1 | while read -r f; do
            if [[ -n "$f" ]]; then
                if [[ "$f" == *.gz ]]; then
                    cp "$f" "$OUT_DIR/expected_r2.fq.gz"
                else
                    gzip -c "$f" > "$OUT_DIR/expected_r2.fq.gz"
                fi
            fi
        done
    fi

    # Cleanup temp
    rm -rf "$temp_dir"

    log_info "Expected outputs generated."
}

# Create metadata file
create_metadata() {
    log_info "Creating metadata.json..."

    local bam_size vcf_size counts_lines
    bam_size=$(stat -c%s "$OUT_DIR/chr21.bam" 2>/dev/null || stat -f%z "$OUT_DIR/chr21.bam")
    vcf_size=$(stat -c%s "$OUT_DIR/chr21.vcf.gz" 2>/dev/null || stat -f%z "$OUT_DIR/chr21.vcf.gz")
    counts_lines=$(wc -l < "$OUT_DIR/expected_counts.tsv" 2>/dev/null || echo "0")

    cat > "$OUT_DIR/metadata.json" << EOF
{
    "version": "$VERSION",
    "created": "$(date -Iseconds)",
    "source": {
        "sample": "HG00731",
        "data_type": "RNA-seq",
        "aligner": "STAR",
        "chromosome": "$CHROMOSOME",
        "source_dir": "$SRC_DIR"
    },
    "files": {
        "bam": {
            "name": "chr21.bam",
            "size_bytes": $bam_size
        },
        "vcf": {
            "name": "chr21.vcf.gz",
            "size_bytes": $vcf_size
        },
        "expected_counts": {
            "name": "expected_counts.tsv",
            "lines": $counts_lines
        }
    },
    "wasp2_version": "$(WASP2 --version 2>&1 | head -1 || echo 'unknown')"
}
EOF

    log_info "Metadata created."
}

# Create README
create_readme() {
    log_info "Creating README.md..."

    cat > "$OUT_DIR/README.md" << 'EOF'
# WASP2 Sanity Test Data (chr21)

## Overview

This directory contains chr21 subset data from HG00731 RNA-seq for WASP2 sanity testing.
The data is used to validate that the WASP2 pipeline produces consistent, reproducible results.

## Files

| File | Description | Size |
|------|-------------|------|
| `chr21.bam` | Chr21 aligned reads (STAR) | ~100MB |
| `chr21.bam.bai` | BAM index | ~100KB |
| `chr21.vcf.gz` | Het variants for chr21 | ~2MB |
| `chr21.vcf.gz.tbi` | VCF index | ~50KB |
| `expected_counts.tsv` | Expected allele counts | ~2MB |
| `expected_r1.fq.gz` | Expected R1 FASTQ for remapping | ~20MB |
| `expected_r2.fq.gz` | Expected R2 FASTQ for remapping | ~20MB |
| `metadata.json` | Data provenance metadata | ~1KB |

## Data Source

- **Sample**: HG00731
- **Data Type**: RNA-seq
- **Aligner**: STAR
- **Original Location**: WASP2 paper figure2 benchmark data

## Usage

### Download (CI/local)
```bash
make download-sanity-data
```

### Run Sanity Tests
```bash
pytest tests/sanity/ -v --tb=short
```

### Regenerate Expected Outputs
```bash
./scripts/create_sanity_data.sh
```

## Statistics

- Reads: ~855K
- Variants: ~37K het SNPs
- Processing time: ~30 seconds

## Version

See `metadata.json` for version and creation details.
EOF

    log_info "README created."
}

# Create tarball for release
create_tarball() {
    log_info "Creating release tarball..."

    local tarball="wasp2-sanity-chr21-$VERSION.tar.xz"
    local tarball_path="$OUT_DIR/../$tarball"

    # Create tarball with compression
    tar -cJf "$tarball_path" \
        -C "$(dirname "$OUT_DIR")" \
        "$(basename "$OUT_DIR")" \
        --transform "s|$(basename "$OUT_DIR")|wasp2-sanity-chr21-$VERSION|"

    local tarball_size
    tarball_size=$(stat -c%s "$tarball_path" 2>/dev/null || stat -f%z "$tarball_path")
    local tarball_mb=$((tarball_size / 1024 / 1024))

    log_info "Tarball created: $tarball_path ($tarball_mb MB)"

    # Generate checksum
    sha256sum "$tarball_path" > "$tarball_path.sha256"
    log_info "Checksum: $(cat "$tarball_path.sha256")"
}

# Main
main() {
    echo "=========================================="
    echo "WASP2 Sanity Data Generation Script"
    echo "=========================================="
    echo ""

    check_prerequisites
    create_output_dir
    extract_chr21_bam
    extract_chr21_vcf
    generate_expected_outputs
    create_metadata
    create_readme
    create_tarball

    echo ""
    log_info "Sanity data generation complete!"
    log_info "Output directory: $OUT_DIR"
    echo ""
    echo "Next steps:"
    echo "  1. Upload tarball to GitHub release"
    echo "  2. Update SANITY_DATA_URL in tests/sanity/conftest.py"
    echo "  3. Run: pytest tests/sanity/ -v"
}

main "$@"
