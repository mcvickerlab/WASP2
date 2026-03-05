# Seqera AI Development Integration

> WASP2 Pipeline Development with AI-Assisted Tooling

This guide documents the integration of Seqera AI tools into the WASP2 pipeline development workflow. It complements Claude Code for complex logic while leveraging Seqera AI's specialized Nextflow DSL2 capabilities.

## Overview

WASP2 pipeline development benefits from a multi-tool AI strategy:

| Tool | Strength | Use Case |
|------|----------|----------|
| **Claude Code** | Architecture, complex logic, code review | Design decisions, debugging, refactoring |
| **Seqera AI** | Nextflow DSL2 syntax, nf-test generation | Process definitions, pipeline scaffolding |
| **Nextflow Tooling** | Environment validation, nf-test, nf-core lint | Pre-flight checks, compliance verification |

## Seqera AI Capabilities

### 1. VS Code Integration

The Seqera AI VS Code extension provides:

- **@Seqera chat**: Nextflow code generation via chat interface
- **Pipeline Mode**: Contextual debugging with error explanations
- **nf-test generation**: Automatic test template creation

Installation:
```bash
# VS Code Extension Marketplace
# Search: "Seqera AI"
# Or visit: https://marketplace.visualstudio.com/items?itemName=seqera.seqera-ai
```

### 2. Chat Features

Use `@Seqera` in VS Code chat for:

```
@Seqera Create a process that runs STAR alignment with WASP filtering
@Seqera Generate nf-test for this module
@Seqera Debug this error: "WASP2_COUNT failed with exit code 1"
```

### 3. Pipeline Mode

Enable Pipeline Mode when working on WASP2 pipelines for:

- Contextual understanding of pipeline structure
- Error message interpretation
- Fix suggestions with DSL2 syntax

## Development Workflow

### Recommended 4-Phase Approach

```
┌─────────────────────────────────────────────────────────────┐
│  Phase 1: DESIGN        │  Phase 2: GENERATE               │
│  ───────────────────    │  ──────────────────              │
│  Tool: Claude Code      │  Tool: Seqera AI                 │
│                         │                                   │
│  • Architecture design  │  • DSL2 process definitions      │
│  • Module structure     │  • Workflow scaffolding          │
│  • Error handling       │  • nf-test templates             │
│  • Complex algorithms   │  • Config generation             │
├─────────────────────────┼───────────────────────────────────│
│  Phase 3: VALIDATE      │  Phase 4: REVIEW                 │
│  ───────────────────    │  ──────────────────              │
│  Tool: Nextflow/nf-test │  Tool: Claude Code               │
│                         │                                   │
│  • nextflow -preview    │  • Security review               │
│  • nf-test execution    │  • Code review                   │
│  • nf-core lint         │  • Performance optimization      │
└─────────────────────────┴───────────────────────────────────┘
```

### Phase 1: Design with Claude Code

Start by designing the module structure:

```bash
# Example: Planning a new WASP2 subworkflow
claude "Design a subworkflow for allelic imbalance analysis that:
- Takes BAM + VCF input
- Runs WASP2_COUNT for allele counting
- Runs WASP2_ANALYZE for statistical testing
- Outputs TSV with beta-binomial p-values"
```

Output: Architecture document with file structure, input/output specs, error handling strategy.

### Phase 2: Generate DSL2 with Seqera AI

Use Seqera AI for Nextflow-specific code:

```
@Seqera Create a DSL2 process called WASP2_COUNT that:
- Takes meta map, bam, bai, vcf as input
- Runs wasp2-count count-variants
- Outputs meta map and counts TSV
- Uses container 'ghcr.io/mcvickerlab/wasp2:latest'
```

Example generated output:
```nextflow
process WASP2_COUNT {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/mcvickerlab/wasp2:latest'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts
    path "versions.yml"                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wasp2-count count-variants \\
        ${bam} \\
        ${vcf} \\
        --output ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: \$(wasp2-count --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp2: 1.0.0
    END_VERSIONS
    """
}
```

Generate nf-test:
```
@Seqera Generate nf-test for the WASP2_COUNT process
```

### Phase 3: Validate Environment

Before running the pipeline, validate the environment using standard Nextflow checks:

```bash
# Verify Nextflow installation and version
nextflow -v

# Check Docker/Singularity availability
docker --version || singularity --version

# Validate pipeline configuration (dry-run)
nextflow run . -profile test -preview

# Verify container accessibility
docker pull ghcr.io/mcvickerlab/wasp2:latest
```

> **Note**: The [Anthropic life-sciences](https://github.com/anthropics/life-sciences) plugin provides additional validation scripts (`check_environment.py`, `generate_samplesheet.py`, `manage_genomes.py`) when installed. See their documentation for setup instructions.

### Phase 4: Review with Claude Code

Final review and integration:

```bash
# Code review
claude "Review this Nextflow module for nf-core compliance:
$(cat modules/local/wasp2_count/main.nf)"

# Integration testing
claude "Help debug this nf-test failure:
$(nf-test test modules/local/wasp2_count/tests/main.nf.test 2>&1)"
```

## Tool Comparison Matrix

| Capability | Claude Code | Seqera AI | Nextflow Tooling |
|------------|:-----------:|:---------:|:----------------:|
| Architecture design | ★★★ | ★ | - |
| DSL2 syntax | ★★ | ★★★ | - |
| nf-test generation | ★★ | ★★★ | ★★ |
| Error debugging | ★★★ | ★★★ | ★ |
| Environment validation | ★ | ★★ | ★★★ |
| Samplesheet validation | ★★ | ★ | ★★★ |
| Code review | ★★★ | ★ | - |
| Complex algorithms | ★★★ | ★ | - |
| nf-core compliance | ★★★ | ★★ | ★★★ |

## Example Workflows

### Creating a New WASP2 Module

1. **Design** (Claude Code):
   ```bash
   claude "Design a WASP2 module for differential allelic imbalance
   between conditions. Include input/output specs and algorithm outline."
   ```

2. **Generate** (Seqera AI):
   ```
   @Seqera Create DSL2 process for differential_allelic_imbalance
   that compares two conditions using beta-binomial regression
   ```

3. **Test** (Seqera AI):
   ```
   @Seqera Generate nf-test with test data for this module
   ```

4. **Validate** (nf-test):
   ```bash
   nextflow run . -profile test -preview
   nf-test test modules/local/differential_ai/tests/
   ```

5. **Review** (Claude Code):
   ```bash
   claude "Review for nf-core compliance and optimize performance"
   ```

### Debugging Pipeline Failures

1. **Enable Pipeline Mode** in VS Code Seqera AI settings

2. **Ask Seqera AI**:
   ```
   @Seqera Pipeline Mode: Debug this error from WASP2_COUNT:
   "Error: VCF index not found for sample001.vcf.gz"
   ```

3. **Complex fixes** with Claude Code:
   ```bash
   claude "The VCF index issue suggests a race condition in our
   subworkflow. Review the channel logic in subworkflows/allelic_analysis.nf"
   ```

## Configuration

### Seqera Platform Integration

For running WASP2 pipelines on Seqera Platform:

```yaml
# seqera.yml
manifest:
  name: 'wasp2/nf-rnaseq'
  version: '1.0.0'
  description: 'WASP2 RNA-seq allelic imbalance pipeline'

compute:
  aws:
    region: 'us-west-2'
    queue: 'wasp2-production'

params:
  outdir: 's3://wasp2-results/rnaseq'
  genome: 'GRCh38'
```

### VS Code Settings

```json
{
  "seqera.ai.pipelineMode": true,
  "seqera.ai.workspace": "wasp2-pipelines",
  "seqera.nextflow.path": "~/.nextflow/bin/nextflow"
}
```

## Best Practices

### 1. Use the Right Tool for the Task

- **Architecture and design** → Claude Code
- **Nextflow syntax and boilerplate** → Seqera AI
- **Complex debugging** → Both tools together
- **Environment validation** → Standard Nextflow tooling

### 2. Leverage Pipeline Mode

Enable Seqera AI Pipeline Mode when:
- Debugging pipeline failures
- Writing new processes
- Understanding error messages

### 3. Maintain nf-core Compliance

All WASP2 modules should follow nf-core standards:
- Meta map propagation
- Proper container definitions
- versions.yml output
- Stub run support
- nf-test coverage

### 4. Document AI Assistance

When using AI-generated code, document the source:

```nextflow
/*
 * Module: WASP2_DIFFERENTIAL_AI
 *
 * Generated: Seqera AI (VS Code)
 * Reviewed: Claude Code
 * Author: WASP2 Team
 */
```

### 5. Security Considerations

When using AI-assisted development:

- **Review all generated code** before committing - AI may introduce insecure patterns
- **Never commit credentials** - Use Nextflow secrets or environment variables for sensitive data
- **Pin container versions** - Avoid `latest` tags in production (`ghcr.io/mcvickerlab/wasp2:1.3.0`)
- **Validate inputs** - AI-generated processes may not include proper input validation
- **Use signed containers** - Enable container signature verification when available
- **Audit dependencies** - Review any new dependencies suggested by AI tools

```nextflow
// SECURE: Pin container version, use secrets
process SECURE_EXAMPLE {
    container 'ghcr.io/mcvickerlab/wasp2:1.3.0'
    secret 'API_KEY'

    // ...
}

// INSECURE: Avoid these patterns
// container 'ghcr.io/mcvickerlab/wasp2:latest'  // Unpinned
// script: "curl ${params.api_key}"  // Exposed credential
```

## References

- [Seqera AI Features](https://seqera.io/blog/seqera-ai-new-features-june-2025/)
- [Seqera AI VS Code Extension](https://seqera.io/blog/seqera-ai--nextflow-vs-code/)
- [Anthropic life-sciences Plugin](https://github.com/anthropics/life-sciences)
- [nf-core Developer Docs](https://nf-co.re/docs/contributing/modules)
- [WASP2 Pipeline Documentation](./WASP2_ECOSYSTEM.md)

## Related Issues

- **Parent**: [EPIC #25 - Nextflow Pipeline Ecosystem](https://github.com/mcvickerlab/WASP2/issues/25)
- **Supports**: [#58 - nf-core subworkflow compliance](https://github.com/mcvickerlab/WASP2/issues/58)
- **Complements**: [#95 - Anthropic life-sciences integration](https://github.com/mcvickerlab/WASP2/issues/95)

---

*Milestone: v1.3.0 - Pipeline Ecosystem*
*Last updated: 2026-02-03*
