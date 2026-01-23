# wasp2-nf-modules

Shared Nextflow DSL2 modules for WASP2 pipelines.

## Modules

| Module | Description |
|--------|-------------|
| `wasp2_count` | Allele counting with Rust acceleration |
| `wasp2_filter` | WASP mapping bias filter |
| `vcf_het` | Extract heterozygous variants |
| `beta_binomial` | Statistical testing |

## Usage

```groovy
include { WASP2_COUNT } from '../nf-modules/modules/wasp2_count'

workflow {
    WASP2_COUNT(bam_ch, vcf_ch, regions_ch)
}
```
