#!/usr/bin/env python3
"""
Validate nf-scatac samplesheet format.
"""

import argparse
import csv
import sys
from pathlib import Path


def validate_samplesheet(samplesheet_path: str) -> bool:
    """
    Validate samplesheet CSV format and content.

    Expected format:
    sample,fragments,cellranger_dir,barcode_tag,chemistry

    Either fragments or cellranger_dir must be provided.
    """
    required_columns = ["sample"]
    source_columns = ["fragments", "cellranger_dir"]
    optional_columns = ["barcode_tag", "chemistry"]
    valid_chemistries = {"10x-atac-v1", "10x-atac-v2", "custom"}

    errors = []
    warnings = []

    with open(samplesheet_path) as f:
        reader = csv.DictReader(f)

        # Check columns
        if not reader.fieldnames:
            print("ERROR: Empty samplesheet or invalid CSV format", file=sys.stderr)
            return False

        for col in required_columns:
            if col not in reader.fieldnames:
                errors.append(f"Missing required column: '{col}'")

        # Check at least one source column exists
        has_source_col = any(col in reader.fieldnames for col in source_columns)
        if not has_source_col:
            errors.append(f"Missing data source column. Need at least one of: {source_columns}")

        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return False

        # Validate rows
        sample_ids = set()
        for row_num, row in enumerate(reader, start=2):
            sample_id = row.get("sample", "").strip()
            fragments = row.get("fragments", "").strip()
            cellranger_dir = row.get("cellranger_dir", "").strip()
            barcode_tag = row.get("barcode_tag", "CB").strip()
            chemistry = row.get("chemistry", "10x-atac-v2").strip()

            # Check sample ID
            if not sample_id:
                errors.append(f"Row {row_num}: Missing sample ID")
            elif sample_id in sample_ids:
                errors.append(f"Row {row_num}: Duplicate sample ID '{sample_id}'")
            else:
                sample_ids.add(sample_id)

            # Validate sample ID characters
            if sample_id and not sample_id.replace("_", "").replace("-", "").isalnum():
                errors.append(
                    f"Row {row_num}: Sample ID '{sample_id}' contains invalid characters (use only alphanumeric, underscore, hyphen)"
                )

            # Check that either fragments or cellranger_dir is provided
            if not fragments and not cellranger_dir:
                errors.append(
                    f"Row {row_num}: Must provide either 'fragments' or 'cellranger_dir' for sample '{sample_id}'"
                )

            # Check fragments file if provided
            if fragments:
                if not fragments.endswith((".tsv", ".tsv.gz")):
                    errors.append(
                        f"Row {row_num}: Fragments file must have extension '.tsv' or '.tsv.gz': {fragments}"
                    )
                elif not Path(fragments).exists():
                    warnings.append(f"Row {row_num}: Fragments file not found: {fragments}")

            # Check cellranger_dir if provided
            if cellranger_dir:
                if not Path(cellranger_dir).exists():
                    warnings.append(
                        f"Row {row_num}: CellRanger directory not found: {cellranger_dir}"
                    )
                elif not Path(cellranger_dir).is_dir():
                    errors.append(
                        f"Row {row_num}: CellRanger path is not a directory: {cellranger_dir}"
                    )

            # Validate barcode_tag format (2 uppercase letters)
            if barcode_tag:
                if len(barcode_tag) != 2 or not barcode_tag.isupper() or not barcode_tag.isalpha():
                    errors.append(
                        f"Row {row_num}: Barcode tag must be exactly 2 uppercase letters: '{barcode_tag}'"
                    )

            # Validate chemistry
            if chemistry and chemistry not in valid_chemistries:
                errors.append(
                    f"Row {row_num}: Invalid chemistry '{chemistry}'. Must be one of: {valid_chemistries}"
                )

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
    parser = argparse.ArgumentParser(description="Validate nf-scatac samplesheet")
    parser.add_argument("samplesheet", help="Path to samplesheet CSV")
    args = parser.parse_args()

    if not Path(args.samplesheet).exists():
        print(f"ERROR: Samplesheet not found: {args.samplesheet}", file=sys.stderr)
        sys.exit(1)

    if validate_samplesheet(args.samplesheet):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
