# WASP2 workflow integrations

The repository includes Nextflow DSL2 pipelines and Galaxy wrappers. Pipeline presence and stub tests do not establish biological validation or nf-core certification.

| Component | Purpose | Current boundary |
|-----------|---------|------------------|
| [nf-rnaseq](nf-rnaseq/) | STAR alignment, WASP correction, allele counting and bulk inference | Optional GTF exon/gene grouping is preserved |
| [nf-atacseq](nf-atacseq/) | ATAC alignment, WASP correction, peak-associated counting and bulk inference | Original region membership is passed to analysis |
| [nf-scatac](nf-scatac/) | Per-cell BAM allele counts or fragment overlap coverage; AnnData output | Fragment coverage has no observed alleles; downstream imbalance inference is not wired for either input path |
| [nf-outrider](nf-outrider/) | OUTRIDER aberrant-expression workflow | Separate expression-outlier model, not an allelic imbalance backend |
| [nf-modules](nf-modules/) | Shared WASP2 and alignment modules | Bulk analysis and single-cell counting are separate modules |
| [Galaxy wrappers](../galaxy/tools/wasp2/) | Mapping, bulk counting and bulk inference | No single-cell Galaxy wrappers are included |

The Python command line intentionally calls native Rust for supported counting and bulk inference. Single-cell inference is implemented in Python/SciPy; the bulk `analyze_imbalance` workflow module does not provide single-cell inference.

## Running workflows

Run from the repository root using the profile appropriate for the installed runtime:

```bash
nextflow run pipelines/nf-rnaseq -profile docker --input samplesheet.csv --vcf variants.vcf.gz --star_index star_index --gtf genes.gtf
nextflow run pipelines/nf-atacseq -profile docker --input samplesheet.csv --vcf variants.vcf.gz --fasta reference.fa
```

Read each pipeline's usage documentation for its required reference files, samplesheet format and supported profiles. Container tags in source describe the intended release; development changes require an image built from the matching source.

## Shared module inventory

- `wasp2/count`, `wasp2/count_alleles`: bulk allele counting.
- `wasp2/count_sc`: per-cell allele counting from BAM.
- `wasp2/analyze`, `wasp2/analyze_imbalance`: bulk inference.
- `wasp2/unified_make_reads`, `wasp2/filter_remapped`: mapping-bias correction.
- `wasp2/vcf_to_bed`: variant-coordinate conversion.
- `wasp2/ml_output`: format conversion; it does not change the statistical model.
- Alignment and BAM-processing modules integrate tools such as STAR, BWA and samtools.

Modules use DSL2 channels and metadata maps. Meaningful validation includes actual process execution and schema/parameter assertions in addition to stub tests. Detailed output formats are documented beside each pipeline.
