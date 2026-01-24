# Expected Test Output Baselines

This directory contains expected output baselines for regression testing.

## Purpose

Baseline files serve as reference outputs for detecting unexpected changes in pipeline behavior. When integration tests run, the actual outputs can be compared against these baselines to catch regressions.

## Generating Baselines

After the first successful integration test run:

1. Run the integration tests:
   ```bash
   cd pipelines/nf-rnaseq
   nf-test test tests/integration.nf.test --profile test_integration,conda
   ```

2. Copy the output counts file to this directory:
   ```bash
   cp .nf-test/tests/*/output/results/counts/SAMPLE1_counts.tsv \
      tests/data/expected/sample1_counts.tsv
   ```

3. Generate checksums for regression validation:
   ```bash
   cd tests/data/expected
   sha256sum *.tsv > checksums.sha256
   ```

## Expected Files

| File | Description |
|------|-------------|
| `sample1_counts.tsv` | Allele counts output from WASP2 count step |
| `checksums.sha256` | SHA256 checksums for all baseline files |

## Updating Baselines

If intentional changes are made to the pipeline that affect output format:

1. Re-run the integration tests
2. Review the changes to ensure they are expected
3. Copy new outputs to this directory
4. Update checksums
5. Commit the updated baselines with a clear commit message explaining the change

## Note

Baseline files should only be updated when:
- The output format intentionally changes
- Bug fixes change expected output
- New features add new columns/fields

Never update baselines to "fix" a failing test without understanding why the output changed.
