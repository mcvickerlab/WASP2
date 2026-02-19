# WASP2 Shared Core Test Data

Unified test dataset for all WASP2 pipelines, Galaxy tools, CLI smoke tests, and container validation.

## Regeneration

```bash
conda activate WASP2_dev2
cd tests/shared_data
bash generate_core_data.sh
```

## Contents

| File | Description | Size |
|------|-------------|------|
| `chr_test.fa` + `.fai` | 20kb synthetic reference genome (single contig `chr_test`) | ~20K |
| `variants.vcf` + `.gz` + `.tbi` | 10 het SNPs across 2 samples (SAMPLE1, SAMPLE2) | ~2K |
| `annotation.gtf` | 2 genes, 6 exons (INTGENE001 + strand, INTGENE002 - strand) | ~1.5K |
| `regions.bed` | Exonic regions from GTF | ~300B |
| `sample{1,2,3}.bam` + `.bai` | BWA-aligned wgsim reads (500 pairs each, seeds 42/43/44) | ~50K each |
| `sample{1,2,3}_R{1,2}.fq.gz` | Compressed FASTQ reads for pipeline input | ~7K each |
| `bwa_index/` | BWA index for chr_test.fa | ~60K |

Total: ~700K

## Genome Layout

```
chr_test (19,800 bp)
├── Gene 1 (INTGENE001, + strand, 500-5500)
│   ├── Exon 1: 500-1500   [SNPs at 750, 1200]
│   ├── Exon 2: 2500-3500  [SNPs at 2800, 3200]
│   └── Exon 3: 4500-5500  [SNP at 5000]
└── Gene 2 (INTGENE002, - strand, 10500-15500)
    ├── Exon 1: 10500-11500 [SNPs at 10800, 11200]
    ├── Exon 2: 12500-13500 [SNPs at 12800, 13200]
    └── Exon 3: 14500-15500 [SNP at 15000]
```

## Downstream Usage

Per-pipeline generators symlink or copy from this directory:

- `pipelines/nf-atacseq/tests/data/generate_test_data.sh`
- `pipelines/nf-scatac/tests/data/generate_test_data.sh`
- `pipelines/nf-outrider/tests/data/generate_test_data.sh`
- `galaxy/tools/wasp2/test-data/generate_test_data.sh`

## Dependencies

samtools, bgzip, tabix, wgsim, bwa, bcftools (all available in `WASP2_dev2` conda env)
