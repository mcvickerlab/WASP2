# Public ATAC bulk compatibility data

These compressed tables contain measured allele counts from GRCh38 chromosome 21
of GM12878 ATAC experiment ENCSR637XSC, replicate-2 BAM
[ENCFF962FMH](https://www.encodeproject.org/files/ENCFF962FMH/), intersected with
[ENCFF452ZBX peaks](https://www.encodeproject.org/files/ENCFF452ZBX/) and matched
[HG001 GIAB NISTv3.3.1 phased genotypes](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.1/GRCh38/).

Source downloads passed MD5 checks. Counting was independently compared with the
unmodified archived Python function at thread settings 1 and 4. There are 27,182
SNP rows and 556 peak-overlapping rows; the complete counting CLI reproduced all
556 peak rows exactly. `metadata.json` records source and table checksums.

Tests execute archived Python commit `41354fc` on these inputs rather than storing
Rust outputs as expected results. The original 32-run matrix supplies wide-integer
DataFrames to the archive, bypassing its UInt16 CSV reader. It therefore establishes
wide-input compatibility, not exact reproduction of legacy file-input rounding.
The bulk tests additionally execute untouched CSV input and a separately labeled
CSV adapter changing only the three count-column dtypes to UInt64. Per-row strict
legacy failures and widened-reference comparisons are retained in test properties;
a legacy failure is not relabeled as a strict parity pass. Phased comparisons select truly phased SNPs.
Linear-reference shape repairs and forwarding of its ignored phased flag are
explicit in the test. Strict statistical tolerances remain `rtol=atol=5e-5`;
counts must match exactly. Sixteen legacy runs plus sixteen donor feature runs
form the 32-run matrix, with NumPy 1.26.4 and SciPy 1.13.1 pinned in Python 3.11 CI.

This is one donor and one chromosome, using published aligned reads. It does not
validate a new WASP remapping pass, multi-donor calibration or biological claims.
