# Bulk Python/Rust compatibility validation

Validation date: 2026-10-14. Changes are proposed in PR #133, which includes
the compatibility changes from PR #135. Neither PR is merged or released.

## Reference and data

The reference is archived **Python**, commit
`41354fc087404ba1fa99542803ee9cb9b7e0e6bf`, with NumPy 1.26.4 and SciPy 1.13.1.
Tests execute its functions; Rust outputs are not used as expected results.

Public real ATAC inputs are GM12878 experiment ENCSR637XSC, replicate-2
BAM ENCFF962FMH, peaks ENCFF452ZBX, and matched HG001 GIAB NISTv3.3.1
GRCh38 genotypes. Source MD5 checks passed. Chromosome 21 contains 1,760,836
alignment records, 27,182 eligible SNPs and 556 peak-overlapping SNP rows.
Counting matches the unmodified archived function at thread settings 1 and 4.
The full counting CLI independently reproduced the 556 peak rows and memberships.
Source links, compressed count tables and checksums are in
[the public-data fixture](../tests/analysis/fixtures/public_atac/README.md).

## Results and reproduction

- The 32-run matrix passes locally on macOS arm64: 16 legacy runs across
  SNP/peak, single/linear, phased/unphased and pseudocount 0/1, plus 16 donor
  feature runs across global/per-donor dispersion. Reported raw mu is compared.
- Counts require exact agreement. Statistical tolerances remain
  `rtol=atol=5e-5`; no proportion flipping or region-specific exceptions occur.
- Every legacy matrix case also executes CLI inference and checks its output.
  Mixed-genotype tests compare the actual archive's successful and error paths.
- Heterogeneous donor tests independently fit dispersion with Python and compare
  donor effects, chi-square LRT p-values and donor-local BH correction. They cover
  overlapping features, both dispersion models, phase modes and pseudocounts.
- Local suite: 443 Python tests passed, 82 skipped, 1 deselected; 132 Rust tests
  passed, 7 ignored. Ruff, formatting and diff checks pass. Rust 1.88 Clippy
  passes with the pre-existing `uninlined_format_args` lint disabled.

In a full clone with the built extension and pinned numerical dependencies:

```bash
python -m pytest tests/analysis/test_public_atac_matrix.py \
  tests/analysis/test_bulk_edge_parity.py \
  tests/analysis/test_donor_archive_parity.py -q
```

Python 3.11 Linux CI pins the numerical dependencies and runs the complete real
matrix. Other Python versions run the remaining compatibility regressions.
Final-head Linux and container outcomes must be reviewed separately in GitHub
Actions; local results do not establish those outcomes.

The first full Linux run at `cb7b834` exposed six dynamic-reference failures:
three public unphased peak cases differed in raw mu, and three heterogeneous
per-donor linear cases differed in fitted statistics. Earlier frozen regressions
passed, so they did not establish full Linux compatibility. Both native amd64
and arm64 container builds, ATAC smoke tests and label checks passed on that head.
These strict failures remain integration gates while platform arithmetic is
investigated. `scripts/diagnose_bulk_numerics.py` compares fused/unfused kernels
and NumPy scalar/vector arithmetic without modifying production source.
The probe confirms that the pinned Linux x86_64 SciPy build uses separate
multiply/add arithmetic at the tested Cephes components, while the macOS arm64
build uses fused arithmetic. Rust now follows that measured compiler policy;
the full strict matrix must verify the resulting inference on each platform.

## Deliberate extensions and reference repairs

Global/per-donor dispersion is the intended extension. Shared dispersion does
not pool effect tests or multiple-testing correction across donors. Archived
linear array shapes require explicit reference adapters; phased linear tests
explicitly forward the phase flag that the archived linear wrapper ignored.
The shape adapter preserves the original sequential probability recurrence;
exact scalar/vector recurrence tests guard against replacing it with a
mathematically equivalent but differently rounded product reduction.

Interior rho values are preserved, including tested values from `1e-12` to
`1-1e-12`. Exact endpoints retain a finite regularization because archived
beta-binomial shapes are undefined there. Long probability products use
log-space fallback after losing normal range; this extends the archive's
underflowing arithmetic rather than reproducing its failure.

Archived phasing checks use distinct original GT tokens and the presence of
`0|1`/`1|0`, then parse every retained first allele. Numeric homozygous tokens
follow that behavior. Malformed mixed tokens raise an error instead of silently
changing a region's likelihood. One public peak SNP has `0/1`: the unfiltered
phased command is rejected by both implementations. Phased comparisons explicitly
select truly phased heterozygous SNPs; unphased comparisons retain all SNPs.

## Integration and release gates

Review compatibility PR #135 first, then refresh #133 against the resulting dev
branch and rerun required checks. No merge or tag is included in this validation.
All active version metadata currently uses 1.5.0. CI/release image builds pass
the package version and source SHA to OCI labels; container smoke tests verify
those labels and Python/Rust package versions and exercise ATAC grouping/phasing.

Review the exact final head's tests and containers before integration. A release
dry run remains necessary after the reviewed code reaches the permitted release
branch, including the separate Singularity preflight. Existing published images
and packages will remain old until a release is actually published.

These comparisons cover one real donor/chromosome and controlled heterogeneous
donors. They do not establish biological calibration, a new WASP remapping pass,
or universal bitwise agreement across numerical environments. A symmetric
unphased likelihood has no identifiable biological direction even when its
reported mu agrees numerically with the pinned reference.
