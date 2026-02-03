# Sanity Test Data Hosting Strategy

## Overview

The WASP2 sanity test uses real chr21 HG00731 data (~35MB compressed) to validate
pipeline reproducibility in CI. This document describes the data hosting strategy.

## Hosting Approach: GitHub Releases + Zenodo

| Source | Purpose | URL Pattern |
|--------|---------|-------------|
| **GitHub Releases** (Primary) | CI testing, fast download | `releases/download/v1.3.0/wasp2-sanity-chr21-v1.tar.xz` |
| **Zenodo** (Archival) | DOI citation, long-term preservation | `zenodo.org/records/XXXXXXX` |

### Why This Approach?

1. **GitHub Releases** (Primary for CI)
   - Integrated with GitHub Actions caching
   - Fast CDN-backed downloads
   - No external dependencies
   - Free for public repositories
   - 2GB file limit (sufficient for 35MB tarball)

2. **Zenodo** (Archival backup)
   - Provides DOI for academic citation
   - CERN-backed long-term preservation
   - GitHub integration for automatic versioning
   - Free, 50GB file limit

### Alternative Options Considered

| Option | Verdict | Reason |
|--------|---------|--------|
| Git LFS | Not used | Bandwidth limits, adds complexity |
| AWS S3 | Not needed | Overkill for 35MB, requires cost management |
| Figshare | Alternative | Similar to Zenodo, less GitHub integration |
| In-repo | Not suitable | Bloats repo, slow clones |

## File Inventory

```
wasp2-sanity-chr21-v1.tar.xz (35MB)
├── chr21.bam           (32MB)  - HG00731 RNA-seq chr21 subset
├── chr21.bam.bai       (46KB)  - BAM index
├── chr21.vcf.gz        (530KB) - Het variants
├── chr21.vcf.gz.tbi    (18KB)  - VCF index
├── expected_counts.tsv (807KB) - Expected allele counts
├── expected_r1.fq.gz   (786KB) - Expected R1 FASTQ
├── expected_r2.fq.gz   (813KB) - Expected R2 FASTQ
├── expected_analysis.tsv (24KB) - Expected analysis output
├── metadata.json       (1.4KB) - Dataset metadata
└── README.md           (1.7KB) - Dataset documentation
```

## Updating the Data

To regenerate sanity data (e.g., after pipeline changes):

```bash
# 1. Generate new expected outputs
cd /path/to/sanity_test
wasp2-count --bam chr21.bam --vcf chr21.vcf.gz --output expected_counts.tsv
# ... (see implementation plan for full commands)

# 2. Create new tarball with incremented version
tar -cJf wasp2-sanity-chr21-v2.tar.xz sanity_test/

# 3. Upload to GitHub release
gh release upload v1.4.0 wasp2-sanity-chr21-v2.tar.xz

# 4. Update SANITY_DATA_VERSION in conftest.py
# 5. Optionally upload to Zenodo for archival DOI
```

## References

- [GitHub Releases documentation](https://docs.github.com/en/repositories/releasing-projects-on-github)
- [Zenodo GitHub integration](https://help.zenodo.org/docs/github/)
- [PHA4GE Pipeline Best Practices](https://github.com/pha4ge/public-health-pipeline-best-practices)
