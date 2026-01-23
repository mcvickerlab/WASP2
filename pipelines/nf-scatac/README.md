# wasp2-nf-scatac

Single-Cell ATAC-seq Allelic Imbalance Pipeline

## Features
- 10x Genomics scATAC support
- Cell-type-resolved allelic analysis
- Pseudo-bulk aggregation
- AnnData/Zarr output for ML

## Quick Start
```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --vcf variants.vcf.gz \
  --cellranger_output /path/to/cellranger/
```

## Integration
- Compatible with ArchR, Signac, SnapATAC2
- Exports to GenVarLoader format
