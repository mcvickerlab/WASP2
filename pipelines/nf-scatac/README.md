# wasp2-nf-scatac

Single-Cell ATAC-seq Allelic Imbalance Pipeline

## Features
- 10x Genomics scATAC fragment support
- Allelic imbalance analysis at heterozygous SNPs
- Pseudo-bulk aggregation per sample

## Quick Start
```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --vcf variants.vcf.gz
```

## Testing

### Stub Tests (CI/CD)
Run fast stub tests that validate workflow structure without real computation:
```bash
# Using nf-test
nf-test test --profile test_stub

# Or direct Nextflow stub run
nextflow run . -profile test_stub -stub-run
```

### Integration Tests (Real Data)
Run full pipeline with GM12878 scATAC-seq data:
```bash
nextflow run . -profile test_real,singularity
```

Test data locations:
- **BAM**: `/iblm/netapp/data3/aho/project_data/wasp2/10x_cellranger_atac/gm12878_el4/`
- **VCF**: `/iblm/netapp/data1/aho/variants/NA12878.vcf.gz`

### Test Structure
```
tests/
├── main.nf.test       # nf-test suite
├── tags.yml           # Test tags for filtering
├── stub/              # Minimal stub test data
│   ├── samplesheet.csv
│   ├── fragments.tsv.gz
│   └── variants.vcf.gz
└── real/              # Integration test samplesheet
    └── samplesheet.csv
```

## Output
- Per-cell allele counts (TSV)
- Pseudo-bulk imbalance results (TSV)

## References
- Issue: [#48](https://github.com/Jaureguy760/WASP2-final/issues/48) - Validation & Test Suite
- Issue: [#32](https://github.com/Jaureguy760/WASP2-final/issues/32) - scATAC Pipeline
- nf-test docs: https://code.askimed.com/nf-test/
