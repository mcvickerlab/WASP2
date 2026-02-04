# Illumination: Performance Transformation
# Episode: 003 - The Rust Metamorphosis

## Performance Comparison

```mermaid
xychart-beta
    title "WASP2 Performance Gains (seconds)"
    x-axis ["BAM-BED Intersect", "Statistical Analysis", "Full Pipeline"]
    y-axis "Time (seconds)" 0 --> 550
    bar [152, 2.7, 500]
    bar [3, 0.5, 50]
```

## The Speedup Table

```mermaid
graph LR
    subgraph "Before: Python"
        P1["BAM-BED: 152s"]
        P2["Analysis: 2.7s"]
        P3["Pipeline: ~500s"]
    end

    subgraph "After: Rust"
        R1["BAM-BED: 2-3s"]
        R2["Analysis: 0.5s"]
        R3["Pipeline: ~50s"]
    end

    subgraph "Speedup"
        S1["50-75x"]
        S2["5x"]
        S3["10x"]
    end

    P1 -.-> S1
    P2 -.-> S2
    P3 -.-> S3
    S1 -.-> R1
    S2 -.-> R2
    S3 -.-> R3

    style P1 fill:#FFB6C1
    style P2 fill:#FFB6C1
    style P3 fill:#FFB6C1
    style R1 fill:#90EE90
    style R2 fill:#90EE90
    style R3 fill:#90EE90
    style S1 fill:#FFD700
    style S2 fill:#FFD700
    style S3 fill:#FFD700
```

## Rust Module Architecture

```mermaid
graph TB
    subgraph "Python Layer"
        CLI["CLI<br/>(Typer)"]
        ORCH["Orchestration"]
        IO["I/O dispatch"]
    end

    subgraph "Rust Layer (via PyO3)"
        BAM_INT["bam_intersect.rs<br/>COITree intervals"]
        BAM_CNT["bam_counter.rs<br/>Parallel counting"]
        BAM_RMP["bam_remapper.rs<br/>CIGAR manipulation"]
        ANAL["analysis.rs<br/>Beta-binomial"]
    end

    CLI --> ORCH
    ORCH --> IO
    IO --> BAM_INT
    IO --> BAM_CNT
    IO --> BAM_RMP
    ORCH --> ANAL

    style CLI fill:#87CEEB
    style BAM_INT fill:#FF8C00
    style BAM_CNT fill:#FF8C00
    style BAM_RMP fill:#FF8C00
    style ANAL fill:#FF8C00
```

## The COITree Secret Weapon

```mermaid
graph TD
    subgraph "Old: pybedtools"
        OLD1["BAM file"] --> OLD2["Write temp BED"]
        OLD2 --> OLD3["bedtools intersect<br/>(subprocess)"]
        OLD3 --> OLD4["Parse output"]
        OLD4 --> OLD5["152 seconds"]
    end

    subgraph "New: COITree"
        NEW1["BAM file"] --> NEW2["Build interval tree<br/>(O(n log n))"]
        NEW2 --> NEW3["Query per read<br/>(O(log n + k))"]
        NEW3 --> NEW4["2-3 seconds"]
    end

    style OLD5 fill:#FFB6C1
    style NEW4 fill:#90EE90
```

## Parallel Processing Architecture

```mermaid
flowchart TB
    subgraph "Input"
        BAM["BAM File"]
    end

    subgraph "Chunking"
        BAM --> C1["Chunk 1<br/>chr1:1-10M"]
        BAM --> C2["Chunk 2<br/>chr1:10M-20M"]
        BAM --> C3["Chunk 3<br/>chr1:20M-30M"]
        BAM --> C4["..."]
    end

    subgraph "Parallel Workers (Rayon)"
        C1 --> W1["Worker 1"]
        C2 --> W2["Worker 2"]
        C3 --> W3["Worker 3"]
        C4 --> W4["Worker N"]
    end

    subgraph "Aggregation"
        W1 --> AGG["Lock-free<br/>aggregation"]
        W2 --> AGG
        W3 --> AGG
        W4 --> AGG
        AGG --> OUT["Final counts"]
    end

    style W1 fill:#FF8C00
    style W2 fill:#FF8C00
    style W3 fill:#FF8C00
    style W4 fill:#FF8C00
    style OUT fill:#90EE90
```

## The Python/Rust Boundary

```mermaid
graph TB
    subgraph "Stays in Python"
        P1["CLI argument parsing"]
        P2["Configuration handling"]
        P3["High-level workflow"]
        P4["User messages"]
        P5["I/O format detection"]
    end

    subgraph "Moves to Rust"
        R1["Inner loops over reads"]
        R2["Interval tree operations"]
        R3["Log-likelihood calculations"]
        R4["CIGAR string parsing"]
        R5["Allele swapping"]
    end

    P3 --> R1
    P3 --> R2
    P3 --> R3

    style P1 fill:#87CEEB
    style P2 fill:#87CEEB
    style P3 fill:#87CEEB
    style P4 fill:#87CEEB
    style P5 fill:#87CEEB
    style R1 fill:#FF8C00
    style R2 fill:#FF8C00
    style R3 fill:#FF8C00
    style R4 fill:#FF8C00
    style R5 fill:#FF8C00
```

## New Capabilities Enabled

```mermaid
mindmap
    root((Rust<br/>Metamorphosis))
        INDEL Support
            Full insertions
            Full deletions
            Not just SNPs
        Multi-Format
            VCF native
            BCF native
            PGEN native
            Auto-detection
        Scale
            Millions of cells
            Streaming processing
            Constant memory
        Statistics
            Beta-binomial
            More accurate
            Proper overdispersion
```

## Format Speedups

```mermaid
graph LR
    subgraph "VCF Parsing"
        V1["Standard Python<br/>1x baseline"]
        V2["cyvcf2 (C-backed)<br/>6.9x faster"]
    end

    subgraph "Genotype Format"
        G1["VCF/BCF<br/>1x baseline"]
        G2["PGEN format<br/>25x faster"]
    end

    V1 -.->|"6.9x"| V2
    G1 -.->|"25x"| G2

    style V2 fill:#90EE90
    style G2 fill:#90EE90
```

## The 80/20 Principle Applied

```mermaid
pie title "Code Distribution vs Runtime Impact"
    "Python (90% of code)" : 5
    "Rust (10% of code)" : 95
```

*10% of the code was responsible for 95% of the runtime. Rewrite those 10%.*

---

## Deployment Ecosystem

```mermaid
graph TB
    subgraph "Source"
        CODE["WASP2<br/>Python + Rust"]
    end

    subgraph "Build Systems"
        MATURIN["maturin<br/>(Rust→Python)"]
        DOCKER["Docker<br/>Multi-stage"]
    end

    subgraph "Distribution"
        PYPI["PyPI<br/>pip install wasp2"]
        BIOCONDA["Bioconda<br/>conda install"]
        DOCKERHUB["Docker Hub<br/>jaureguy760/wasp2"]
        SINGULARITY["Singularity<br/>HPC clusters"]
    end

    subgraph "Workflows"
        NF_RNA["nf-rnaseq"]
        NF_ATAC["nf-atacseq"]
        NF_SC["nf-scatac"]
        NF_OUT["nf-outrider"]
    end

    CODE --> MATURIN
    CODE --> DOCKER
    MATURIN --> PYPI
    MATURIN --> BIOCONDA
    DOCKER --> DOCKERHUB
    DOCKERHUB --> SINGULARITY

    PYPI --> NF_RNA
    DOCKERHUB --> NF_RNA
    PYPI --> NF_ATAC
    DOCKERHUB --> NF_ATAC
    PYPI --> NF_SC
    DOCKERHUB --> NF_SC
    PYPI --> NF_OUT
    DOCKERHUB --> NF_OUT

    style CODE fill:#87CEEB
    style PYPI fill:#90EE90
    style BIOCONDA fill:#90EE90
    style DOCKERHUB fill:#FF8C00
    style SINGULARITY fill:#FF8C00
```

## Nextflow Pipeline Architecture

```mermaid
flowchart LR
    subgraph "Input"
        BAM["BAM files"]
        VCF["VCF/BCF/PGEN"]
        META["Sample sheet"]
    end

    subgraph "Nextflow Pipeline"
        NF["nextflow run<br/>wasp2/nf-rnaseq"]

        subgraph "Processes"
            P1["WASP2_COUNT"]
            P2["WASP2_MAP"]
            P3["WASP2_ANALYZE"]
        end
    end

    subgraph "Execution"
        LOCAL["Local"]
        SLURM["SLURM"]
        AWS["AWS Batch"]
        DOCKER2["Docker"]
        SING["Singularity"]
    end

    subgraph "Output"
        COUNTS["Allele counts"]
        FILTERED["Filtered BAMs"]
        RESULTS["QTL results"]
        REPORT["MultiQC report"]
    end

    BAM --> NF
    VCF --> NF
    META --> NF
    NF --> P1
    P1 --> P2
    P2 --> P3

    NF --> LOCAL
    NF --> SLURM
    NF --> AWS
    NF --> DOCKER2
    NF --> SING

    P1 --> COUNTS
    P2 --> FILTERED
    P3 --> RESULTS
    P3 --> REPORT

    style NF fill:#87CEEB
    style SLURM fill:#FFD700
    style SING fill:#FF8C00
```

---

## Episode Reference
- **Episode**: 003 - The Rust Metamorphosis
- **Topic**: Rust acceleration and deployment ecosystem (2024-2026)
- **Version**: 1.3.0
- **Rust Lines**: 10,551+
- **Pipelines**: nf-rnaseq, nf-atacseq, nf-scatac, nf-outrider
