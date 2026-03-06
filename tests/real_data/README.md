# WASP2 Real Data Test Setup

Real genomic data from the 1000 Genomes Project for validating WASP2 pipelines beyond the synthetic test data in `tests/shared_data/`.

## Quick Start

```bash
# Preview what will be downloaded (no actual download)
DRY_RUN=1 bash tests/real_data/download_chr21.sh

# Download with default sample (NA12878)
bash tests/real_data/download_chr21.sh

# Download with HG00731 instead
SAMPLE=HG00731 bash tests/real_data/download_chr21.sh

# Skip STAR index (saves 1.5 GB, only needed for nf-rnaseq)
SKIP_STAR_INDEX=1 bash tests/real_data/download_chr21.sh
```

## Data Sources

| Resource | Source | URL | Size |
|----------|--------|-----|------|
| Reference FASTA | UCSC hg38 chr21 | `hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr21.fa.gz` | 12 MB compressed, ~40 MB uncompressed |
| VCF (phased) | 1000G 2022 SHAPEIT5 panel | `ftp.1000genomes.ebi.ac.uk/.../20220422_3202_phased_SNV_INDEL_SV/` | 407 MB full (subset to ~5 MB) |
| VCF (alternate) | 1000G 2020 SHAPEIT2 panel | `ftp.1000genomes.ebi.ac.uk/.../20201028_3202_phased/` | 473 MB full |
| NA12878 CRAM | NYGC 30x, ENA ERR3239334 | `ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram` | ~30 GB full (stream chr21 only) |
| HG00731 CRAM | HGSV SV discovery | `ftp.1000genomes.ebi.ac.uk/.../hgsv_sv_discovery/data/PUR/HG00731/` | ~30 GB full (stream chr21 only) |
| GRCh38 reference | 1000G technical reference | `ftp.1000genomes.ebi.ac.uk/.../GRCh38_full_analysis_set_plus_decoy_hla.fa` | 3 GB (temp, for CRAM decode) |

## Samples

| Sample | Population | Description | Why Use It |
|--------|------------|-------------|------------|
| **NA12878** | CEU (Utah/European) | GIAB reference, most benchmarked human genome | Gold standard, widely validated truth sets available |
| **HG00731** | PUR (Puerto Rican) | 1000G trio child, used in WASP2 benchmarks | Consistent with existing `test_benchmarks/` data |

## Disk Space

| Component | Size | Required For |
|-----------|------|-------------|
| chr21 reference FASTA | ~40 MB | All pipelines |
| BWA index | ~50 MB | nf-atacseq |
| chr21 BAM (30x) | ~1-2 GB | nf-outrider (direct), others (via FASTQ extract) |
| chr21 FASTQ (paired) | ~500 MB | nf-atacseq, nf-rnaseq |
| chr21 VCF (het SNPs) | ~5 MB | All pipelines |
| STAR index | ~1.5 GB | nf-rnaseq only |
| **Total minimum** | **~2 GB** | Without STAR index |
| **Total with STAR** | **~3.5 GB** | With STAR index |

During download, peak disk usage is ~5-6 GB due to the temporary GRCh38 full reference needed for CRAM-to-BAM conversion. The full reference is deleted after extraction.

## Directory Structure (after download)

```
tests/real_data/
├── download_chr21.sh          # This download script
├── README.md                  # This file
├── data/                      # Downloaded data (gitignored)
│   ├── chr21.fa               # Reference FASTA
│   ├── chr21.fa.fai           # FASTA index
│   ├── chr21.dict             # Sequence dictionary
│   ├── bwa_index/             # BWA index for chr21
│   ├── star_index/            # STAR index for chr21 (optional)
│   ├── NA12878_chr21.bam      # Chr21 reads
│   ├── NA12878_chr21.bam.bai  # BAM index
│   ├── NA12878_chr21_R1.fastq.gz  # Paired FASTQ (extracted from BAM)
│   ├── NA12878_chr21_R2.fastq.gz
│   ├── NA12878_chr21_het.vcf.gz   # Het SNPs only
│   ├── NA12878_chr21_het.vcf.gz.tbi
│   ├── NA12878_chr21_all.vcf.gz   # All variant genotypes
│   └── NA12878_chr21_all.vcf.gz.tbi
├── samplesheets/              # Per-pipeline samplesheets
│   ├── atacseq_samplesheet.csv
│   ├── rnaseq_samplesheet.csv
│   ├── outrider_samplesheet.csv
│   └── scatac_samplesheet.csv     # Template only (needs 10x data)
└── configs/                   # Nextflow configs for each pipeline
    ├── test_real_atacseq.config
    ├── test_real_rnaseq.config
    ├── test_real_outrider.config
    └── test_real_scatac.config    # Template only
```

## Running Pipelines with Real Data

After downloading, run each pipeline:

```bash
# ATAC-seq (bulk)
nextflow run pipelines/nf-atacseq \
    -c tests/real_data/configs/test_real_atacseq.config \
    -profile docker

# RNA-seq ASE
nextflow run pipelines/nf-rnaseq \
    -c tests/real_data/configs/test_real_rnaseq.config \
    -profile docker

# OUTRIDER (expression outlier detection)
nextflow run pipelines/nf-outrider \
    -c tests/real_data/configs/test_real_outrider.config \
    -profile docker
```

## Pipeline-Specific Notes

### nf-atacseq
Works directly with the downloaded data. Uses FASTQ input with BWA alignment.

### nf-rnaseq
The downloaded data is WGS (whole-genome sequencing), not RNA-seq. This is useful for testing pipeline mechanics (alignment, WASP2 filtering, counting) but will not produce biologically meaningful results since WGS reads do not reflect transcription patterns. For biological validation, use GTEx or ENCODE chr21-subsetted RNA-seq BAMs.

### nf-outrider
Requires multiple samples for outlier detection statistics. Download data for several samples by running the script multiple times with different `SAMPLE` values and combining samplesheet rows. Minimum recommended: 5 samples.

### nf-scatac
Cannot use bulk WGS data. Requires 10x Chromium scATAC-seq data processed by CellRanger-ATAC (fragments.tsv.gz). The generated samplesheet is a template that must be updated with real scATAC-seq data paths. The existing `test_real.config` in `pipelines/nf-scatac/conf/` references GM12878 (NA12878 cell line) scATAC-seq data on the IBLM filesystem.

## VCF Details

Two VCF subsets are created per sample:

1. **`{SAMPLE}_chr21_het.vcf.gz`** -- Heterozygous SNPs only. This is what WASP2 needs for allelic bias correction (het sites are where reference bias manifests).

2. **`{SAMPLE}_chr21_all.vcf.gz`** -- All variant genotypes (het + hom-alt, SNPs + indels). Useful for broader analyses or pipeline steps that use all known variants.

The source VCF is the 2022 SHAPEIT5-phased panel from the 1000 Genomes Project (3,202 samples, GRCh38 coordinates). This panel has higher phasing accuracy than the earlier 2020 SHAPEIT2 panel.

## Comparison with Existing Test Data

| Aspect | `tests/shared_data/` (synthetic) | `tests/real_data/` (this) |
|--------|----------------------------------|---------------------------|
| Purpose | CI/CD, fast smoke tests | Full pipeline validation |
| Genome | 20 KB synthetic contig | chr21 (~46 MB, real) |
| Variants | 10 synthetic het SNPs | ~30K-50K real het SNPs |
| Reads | 500 simulated pairs | Millions of real reads |
| Coverage | ~5x synthetic | ~30x real |
| Runtime | Seconds | 10-60 minutes |
| Disk | ~700 KB | ~2-3.5 GB |
| CI suitable | Yes | No (too large, too slow) |

## Prerequisites

All tools available in the `WASP2_dev2` conda environment:

```bash
conda activate WASP2_dev2
# Or install individually:
#   samtools >= 1.10, bcftools, tabix, bgzip, bwa, curl
#   Optional: STAR (for nf-rnaseq index)
```

## Troubleshooting

**CRAM download hangs or times out:**
The CRAM files are served from EBI/ENA FTP servers. If `samtools view` hangs when streaming the remote CRAM, try downloading the full CRAM first:
```bash
curl -L -o NA12878.final.cram "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram"
curl -L -o NA12878.final.cram.crai "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram.crai"
# Then subset locally:
samtools view -b -T GRCh38_ref.fa NA12878.final.cram chr21 -o NA12878_chr21.bam
```

**"Failed to open reference" during CRAM conversion:**
CRAM files require the exact reference genome used during alignment. The script downloads the 1000 Genomes GRCh38 full reference for this purpose. If the download fails, get it manually:
```bash
curl -L -o GRCh38_full_analysis_set_plus_decoy_hla.fa \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
```

**VCF sample not found:**
The 2022 SHAPEIT5 panel contains 3,202 samples from the 1000 Genomes Project. Verify your sample is included:
```bash
bcftools query -l 1kGP_chr21_full.vcf.gz | grep YOUR_SAMPLE
```
