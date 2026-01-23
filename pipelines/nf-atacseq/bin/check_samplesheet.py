#!/usr/bin/env python3
"""
Validate nf-atacseq samplesheet format.
"""

import argparse
import csv
import sys
from pathlib import Path


def validate_samplesheet(samplesheet_path: str) -> bool:
    """
    Validate samplesheet CSV format and content.

    Expected format:
    sample,fastq_1,fastq_2,sample_name
    """
    required_columns = ['sample', 'fastq_1', 'fastq_2']
    optional_columns = ['sample_name', 'single_end']

    errors = []
    warnings = []

    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)

        # Check columns
        if not reader.fieldnames:
            print("ERROR: Empty samplesheet or invalid CSV format", file=sys.stderr)
            return False

        for col in required_columns:
            if col not in reader.fieldnames:
                errors.append(f"Missing required column: '{col}'")

        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return False

        # Validate rows
        sample_ids = set()
        for row_num, row in enumerate(reader, start=2):
            sample_id = row.get('sample', '').strip()
            fastq_1 = row.get('fastq_1', '').strip()
            fastq_2 = row.get('fastq_2', '').strip()

            # Check sample ID
            if not sample_id:
                errors.append(f"Row {row_num}: Missing sample ID")
            elif sample_id in sample_ids:
                errors.append(f"Row {row_num}: Duplicate sample ID '{sample_id}'")
            else:
                sample_ids.add(sample_id)

            # Validate sample ID characters
            if sample_id and not sample_id.replace('_', '').replace('-', '').isalnum():
                errors.append(f"Row {row_num}: Sample ID '{sample_id}' contains invalid characters (use only alphanumeric, underscore, hyphen)")

            # Check FASTQ files
            if not fastq_1:
                errors.append(f"Row {row_num}: Missing fastq_1 for sample '{sample_id}'")
            else:
                if not Path(fastq_1).exists():
                    warnings.append(f"Row {row_num}: fastq_1 file not found: {fastq_1}")

            if not fastq_2:
                # Single-end data
                pass
            else:
                if not Path(fastq_2).exists():
                    warnings.append(f"Row {row_num}: fastq_2 file not found: {fastq_2}")

    # Print results
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)

    if errors:
        return False

    print(f"Samplesheet validation passed: {len(sample_ids)} samples", file=sys.stderr)
    return True


def main():
    parser = argparse.ArgumentParser(description='Validate nf-atacseq samplesheet')
    parser.add_argument('samplesheet', help='Path to samplesheet CSV')
    args = parser.parse_args()

    if not Path(args.samplesheet).exists():
        print(f"ERROR: Samplesheet not found: {args.samplesheet}", file=sys.stderr)
        sys.exit(1)

    if validate_samplesheet(args.samplesheet):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
