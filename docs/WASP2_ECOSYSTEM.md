# WASP2 Nextflow Pipeline Ecosystem

> Tracking document for [EPIC #25](https://github.com/Jaureguy760/WASP2-final/issues/25)

## Status Matrix

### Core Pipelines

| Component | Issue | Status | Infrastructure |
|-----------|-------|--------|----------------|
| wasp2-nf-modules | [#29](../../issues/29) | ✅ Complete | 9 modules, nf-test |
| wasp2-nf-rnaseq | [#30](../../issues/30) | ✅ Complete | docs, tests, assets |
| wasp2-nf-atacseq | [#31](../../issues/31) | ✅ Complete | docs, tests, assets, bin |
| wasp2-nf-scatac | [#32](../../issues/32) | ✅ Complete | docs, tests, assets, bin |
| wasp2-nf-outrider | [#35](../../issues/35) | ✅ Complete | docs, tests, assets, bin |

### Integrations

| Component | Issue | Status |
|-----------|-------|--------|
| ML Output Formats | [#36](../../issues/36) | ✅ Complete |
| GenVarLoader | [#37](../../issues/37) | ✅ Complete |
| nf-core Compliance | [#38](../../issues/38) | ✅ Complete |
| Seqera AI | [#39](../../issues/39) | ✅ Complete |

## Module Inventory

| Module | Function | Performance |
|--------|----------|-------------|
| WASP2_COUNT | Allelic read counting | Rust: 61× faster |
| WASP2_MAP | Read remapping/filtering | Rust: 5× faster |
| WASP2_ANALYZE | Statistical analysis | Rust-backed |
| WASP2_COUNT_ALLELES | Single-cell counting | Rust |
| WASP2_ANALYZE_IMBALANCE | SC imbalance | Rust |
| WASP2_ML_OUTPUT | ML format conversion | Zarr, Parquet, AnnData |
| VCF_TO_BED | VCF conversion | Rust: 7-25× faster |
| STAR_ALIGN | STAR 2-pass | Native |

## Pipeline Directory Structure

All pipelines follow a consistent nf-core-inspired structure:

```
pipelines/
├── nf-modules/              # Shared DSL2 modules
│   └── modules/wasp2/       # WASP2-specific modules
├── nf-rnaseq/               # RNA-seq allelic imbalance
├── nf-atacseq/              # ATAC-seq allelic imbalance
├── nf-scatac/               # Single-cell ATAC-seq AI
│   ├── main.nf
│   ├── nextflow.config
│   ├── workflows/
│   ├── subworkflows/
│   ├── modules/local/
│   ├── conf/
│   ├── assets/              # samplesheet schema, multiqc config
│   ├── bin/                 # helper scripts
│   ├── docs/                # usage.md, output.md
│   └── tests/               # nf-test, stub data
└── nf-outrider/             # OUTRIDER aberrant expression
    ├── main.nf
    ├── nextflow.config
    ├── workflows/
    ├── subworkflows/
    ├── modules/local/
    ├── conf/
    ├── assets/
    ├── bin/
    ├── docs/
    └── tests/
```

## Dependency Graph

```
                wasp2-nf-modules (#29) ✅
                        │
        ┌───────────────┼───────────────┐
        ▼               ▼               ▼
   nf-rnaseq ✅    nf-atacseq ✅   ML Formats ✅
        │               │               │
        ▼               ▼               ▼
   nf-outrider ✅   nf-scatac ✅   GenVarLoader ✅
                        │
                        ▼
                nf-core Compliance ✅
                        │
                        ▼
                   Seqera AI ✅
```

## Implementation Roadmap

### Phase 1: Foundation ✅
- [x] Core DSL2 modules (9 modules)
- [x] nf-rnaseq pipeline
- [x] nf-atacseq pipeline
- [x] Docker builds, nf-test infrastructure

### Phase 2: Expansion ✅
- [x] nf-scatac (#32) - Single-cell ATAC-seq allelic imbalance
- [x] nf-outrider (#35) - OUTRIDER aberrant expression + MAE
- [x] ML output formats (#36) - Zarr, Parquet, AnnData

### Phase 3: Integration ✅
- [x] GenVarLoader integration (#37) - Via Zarr output format
- [x] nf-core compliance (#38) - Pipeline structure compliance
- [x] Seqera AI compatibility (#39) - [Integration guide](./source/seqera_ai_integration.md)

## ML Output Formats

All pipelines support optional ML-ready output formats via the `--output_format` parameter:

```bash
# Single format
nextflow run . --output_format zarr

# Multiple formats (comma-separated)
nextflow run . --output_format zarr,parquet,anndata
```

### Available Formats

| Format | Description | Ecosystem |
|--------|-------------|-----------|
| **Zarr** | Chunked cloud-native arrays | GenVarLoader, xarray |
| **Parquet** | Columnar analytics format | Polars, DuckDB, pandas |
| **AnnData** | H5AD with layers | Scanpy, ArchR, scverse |

### GenVarLoader Compatibility

Zarr outputs are directly compatible with [GenVarLoader](https://genvarloader.readthedocs.io/) for ML training:

```python
import genvarloader as gvl
loader = gvl.VariantLoader(zarr_path="sample.zarr")
```

## Testing

All pipelines support:
- **Stub tests**: Fast CI/CD validation with `-profile test_stub -stub-run`
- **Integration tests**: Real data with `-profile test_real` or `-profile test`
- **nf-test framework**: Modular testing at workflow, subworkflow, and module levels

Run stub tests:
```bash
cd pipelines/nf-scatac
nextflow run . -profile test_stub -stub-run

cd pipelines/nf-outrider
nextflow run . -profile test_stub -stub-run
```

## References

- [nf-core/drop](https://nf-co.re/drop/dev/) - Reference OUTRIDER implementation
- [GenVarLoader](https://genvarloader.readthedocs.io/) - ML variant loading
- [Seqera AI](https://seqera.io/blog/seqera-ai-new-features-june-2025/) - Pipeline AI assistant

---
*Milestone: v1.3.0 - Pipeline Ecosystem*
*Last updated: 2026-02-03*
