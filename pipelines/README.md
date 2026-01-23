# WASP2 Nextflow Pipelines

Modular Nextflow DSL2 pipelines for allele-specific analysis.

## Pipelines

| Pipeline | Description | Status |
|----------|-------------|--------|
| **nf-rnaseq** | RNA-seq allele-specific expression (ASE) | 🚧 Planned |
| **nf-atacseq** | ATAC-seq allelic imbalance (AI) | 🚧 Planned |
| **nf-scatac** | Single-cell ATAC-seq AI | 🚧 Planned |
| **nf-modules** | Shared DSL2 modules | 🚧 Planned |

## Architecture

```
pipelines/
├── nf-modules/          # Shared modules (WASP2 counting, filtering)
│   └── modules/
│       ├── wasp2_count/
│       ├── wasp2_filter/
│       └── vcf_processing/
├── nf-rnaseq/           # RNA-seq ASE pipeline
│   ├── main.nf
│   ├── nextflow.config
│   └── conf/
├── nf-atacseq/          # ATAC-seq AI pipeline
└── nf-scatac/           # Single-cell ATAC pipeline
```

## Usage

```bash
# RNA-seq ASE
nextflow run pipelines/nf-rnaseq -profile docker --input samplesheet.csv

# ATAC-seq AI
nextflow run pipelines/nf-atacseq -profile singularity --input samplesheet.csv
```

## nf-core Compatibility

These pipelines follow nf-core standards where practical:
- DSL2 modules with meta maps
- MultiQC integration
- Conda/Docker/Singularity support
- Tower compatibility
