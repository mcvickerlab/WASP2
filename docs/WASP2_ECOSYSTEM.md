# WASP2 Nextflow Pipeline Ecosystem

> Tracking document for [EPIC #25](https://github.com/Jaureguy760/WASP2-final/issues/25)

## Status Matrix

### Core Pipelines

| Component | Issue | Status |
|-----------|-------|--------|
| wasp2-nf-modules | [#29](../../issues/29) | ✅ Complete |
| wasp2-nf-rnaseq | [#30](../../issues/30) | ✅ Complete |
| wasp2-nf-atacseq | [#31](../../issues/31) | ✅ Complete |
| wasp2-nf-scatac | [#32](../../issues/32) | 🔄 In Progress |
| wasp2-nf-outrider | [#35](../../issues/35) | 🔄 Open |

### Integrations

| Component | Issue | Status |
|-----------|-------|--------|
| ML Output Formats | [#36](../../issues/36) | 🔄 Open |
| GenVarLoader | [#37](../../issues/37) | 🔄 Open |
| nf-core Compliance | [#38](../../issues/38) | 🔄 Open |
| Seqera AI | [#39](../../issues/39) | 🔄 Open |

## Module Inventory

| Module | Function | Performance |
|--------|----------|-------------|
| WASP2_COUNT | Allelic read counting | Rust: 5-10x |
| WASP2_MAP | Read remapping/filtering | Rust: 5x |
| WASP2_ANALYZE | Statistical analysis | Rust-backed |
| WASP2_COUNT_ALLELES | Single-cell counting | Rust |
| WASP2_ANALYZE_IMBALANCE | SC imbalance | Rust |
| VCF_TO_BED | VCF conversion | Rust: 7-25x |
| STAR_ALIGN | STAR 2-pass | Native |

## Dependency Graph

```
                wasp2-nf-modules (#29) ✅
                        │
        ┌───────────────┼───────────────┐
        ▼               ▼               ▼
   nf-rnaseq ✅    nf-atacseq ✅   ML Formats
        │               │               │
        ▼               ▼               ▼
   nf-outrider     nf-scatac       GenVarLoader
                        │
                        ▼
                nf-core Compliance
                        │
                        ▼
                   Seqera AI
```

## Implementation Roadmap

### Phase 1: Foundation ✅
- [x] Core DSL2 modules
- [x] nf-rnaseq, nf-atacseq pipelines
- [x] Docker builds, nf-test

### Phase 2: Expansion (Current)
- [ ] nf-scatac (#32)
- [ ] nf-outrider (#35)
- [ ] ML output formats (#36)

### Phase 3: Integration
- [ ] GenVarLoader (#37)
- [ ] nf-core compliance (#38)
- [ ] Seqera AI (#39)

## References

- [nf-core/drop](https://nf-co.re/drop/dev/)
- [GenVarLoader](https://genvarloader.readthedocs.io/)
- [Seqera AI](https://seqera.io/blog/seqera-ai-new-features-june-2025/)

---
*Milestone: v1.3.0 - Pipeline Ecosystem*
