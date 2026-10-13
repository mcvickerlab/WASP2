`archive_count_alleles.py.txt` is a verbatim copy of
`src/counting/count_alleles.py` from WASP2 commit
`41354fc087404ba1fa99542803ee9cb9b7e0e6bf` (`archive/master-pre-v1.3.0`).

The counting parity tests compile only `count_snp_alleles` from this source.
They execute its original SNP-ordered BAM fetches and chromosome-wide read-name
set without translating the algorithm or changing expected counts.

The Rust bulk biallelic SNP path now uses the same per-SNP fetch order. This can
require more index seeks and repeated BAM decoding than a chromosome-span scan;
large-cohort performance has not been measured here. Non-SNP counting keeps its
existing implementation. Count storage remains UInt32 to avoid legacy UInt16
overflow, and an explicitly requested base-quality threshold remains available;
archive comparisons use the default `min_qual=0`.
