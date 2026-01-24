# WASP2 Galaxy Tools

Galaxy tool wrappers for WASP2 allele-specific analysis.

## Tools Included

| Tool | Description | CLI Command |
|------|-------------|-------------|
| `wasp2_count_variants` | Count allele-specific reads | `wasp2-count count-variants` |
| `wasp2_make_reads` | Generate reads for WASP remapping | `wasp2-map make-reads` |
| `wasp2_filter_remapped` | Filter remapped reads | `wasp2-map filter-remapped` |
| `wasp2_find_imbalance` | Statistical analysis | `wasp2-analyze find-imbalance` |

## Installation

### From Galaxy Tool Shed

```
Search for "wasp2" in the Galaxy Tool Shed
```

### Manual Installation

1. Copy the `wasp2` directory to your Galaxy `tools/` directory
2. Add to `tool_conf.xml`:

```xml
<section id="wasp2" name="WASP2 Allele-Specific">
    <tool file="wasp2/wasp2_count_variants.xml"/>
    <tool file="wasp2/wasp2_make_reads.xml"/>
    <tool file="wasp2/wasp2_filter_remapped.xml"/>
    <tool file="wasp2/wasp2_find_imbalance.xml"/>
</section>
```

3. Install the wasp2 conda package:
```bash
conda install -c bioconda wasp2
```

## Testing with Planemo

```bash
# Install planemo
pip install planemo

# Lint tools
planemo lint wasp2/

# Run tests
planemo test wasp2/

# Serve locally for testing
planemo serve wasp2/
```

## Workflow: WASP Bias Correction

```
                    ┌─────────────────┐
                    │   Input BAM     │
                    │   + VCF         │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │ WASP2 Make      │
                    │ Reads           │
                    └────────┬────────┘
                             │
            ┌────────────────┼────────────────┐
            ▼                ▼                ▼
     ┌───────────┐    ┌───────────┐    ┌───────────┐
     │ to_remap  │    │ to_remap  │    │  keep     │
     │ R1.fq.gz  │    │ R2.fq.gz  │    │  .bam     │
     └─────┬─────┘    └─────┬─────┘    └─────┬─────┘
           │                │                │
           └───────┬────────┘                │
                   ▼                         │
           ┌───────────────┐                 │
           │  Your Aligner │                 │
           │  (BWA, STAR)  │                 │
           └───────┬───────┘                 │
                   ▼                         │
           ┌───────────────┐                 │
           │ remapped.bam  │                 │
           └───────┬───────┘                 │
                   │                         │
                   └─────────┬───────────────┘
                             ▼
                    ┌────────────────┐
                    │ WASP2 Filter   │
                    │ Remapped       │
                    └────────┬───────┘
                             ▼
                    ┌────────────────┐
                    │ filtered.bam   │
                    │ (unbiased)     │
                    └────────┬───────┘
                             │
                    ┌────────▼────────┐
                    │ WASP2 Count     │
                    │ Variants        │
                    └────────┬────────┘
                             ▼
                    ┌────────────────┐
                    │ WASP2 Find     │
                    │ Imbalance      │
                    └────────────────┘
```

## Test Data

Test data files should be placed in `test-data/`:
- `test.bam` + `test.bam.bai`: Small aligned BAM
- `test.vcf`: Matching VCF with heterozygous variants
- `counts.tsv`: Example count output
- `wasp_data.json`: Example WASP metadata

## Support

- Issues: https://github.com/Jaureguy760/WASP2-exp/issues
- Documentation: https://Jaureguy760.github.io/WASP2-exp/
