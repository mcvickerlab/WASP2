# Illumination: WASP2 Architecture
# Episode: 002 - Building the New Hive

## The Modernization: Old vs New

```mermaid
graph LR
    subgraph "Original WASP (2015)"
        O_VCF["VCF"] --> O_CONV["snp2h5<br/>(conversion)"]
        O_CONV --> O_H5["HDF5"]
        O_H5 --> O_SCRIPTS["Multiple<br/>Python scripts"]
        O_SCRIPTS --> O_OUT["Output"]
    end

    subgraph "WASP2 (2021+)"
        N_VCF["VCF/BCF"] --> N_CLI["Unified CLI<br/>(wasp2-*)"]
        N_CLI --> N_OUT["Parquet/AnnData"]
    end

    style O_CONV fill:#FFB6C1
    style O_H5 fill:#FFA07A
    style N_CLI fill:#90EE90
```

## Module Organization

```mermaid
graph TB
    subgraph "src/wasp2/"
        CLI["cli/<br/>Typer-based commands"]
        COUNT["counting/<br/>Allele counting"]
        MAP["mapping/<br/>Read filtering"]
        ANAL["analysis/<br/>Statistical tests"]
        IO["io/<br/>Format handlers"]
    end

    CLI --> COUNT
    CLI --> MAP
    CLI --> ANAL
    COUNT --> IO
    MAP --> IO
    ANAL --> IO

    style CLI fill:#87CEEB
    style COUNT fill:#98FB98
    style MAP fill:#DDA0DD
    style ANAL fill:#FFD700
    style IO fill:#F0E68C
```

## The Unified CLI

```mermaid
flowchart LR
    subgraph "Command Structure"
        WASP2["wasp2"]
        WASP2 --> COUNT["wasp2-count<br/>Allele counting"]
        WASP2 --> MAP["wasp2-map<br/>Bias correction"]
        WASP2 --> ANALYZE["wasp2-analyze<br/>QTL discovery"]
    end

    subgraph "Features"
        COUNT --> F1["• VCF/BCF native<br/>• No conversion<br/>• Parquet output"]
        MAP --> F2["• WASP filtering<br/>• Multi-sample<br/>• Remapping"]
        ANALYZE --> F3["• CHT<br/>• Beta-binomial<br/>• Single-cell"]
    end

    style COUNT fill:#98FB98
    style MAP fill:#DDA0DD
    style ANALYZE fill:#FFD700
```

## Data Flow

```mermaid
flowchart TB
    subgraph "Inputs"
        BAM["BAM/CRAM<br/>Alignments"]
        VCF["VCF/BCF<br/>Variants"]
        META["Sample<br/>Metadata"]
    end

    subgraph "WASP2 Processing"
        BAM --> WASP["WASP2"]
        VCF --> WASP
        META --> WASP

        WASP --> FILT["Filtered reads<br/>(bias removed)"]
        WASP --> COUNTS["Allele counts<br/>(per variant)"]
        WASP --> STATS["Statistical tests<br/>(QTL calls)"]
    end

    subgraph "Outputs"
        FILT --> O_BAM["Filtered BAM"]
        COUNTS --> O_PQ["Parquet tables"]
        STATS --> O_RES["Results TSV"]
        COUNTS --> O_AD["AnnData<br/>(single-cell)"]
    end

    style WASP fill:#87CEEB
    style O_PQ fill:#90EE90
    style O_AD fill:#FFD700
```

## Technology Stack Comparison

```mermaid
graph TB
    subgraph "Original WASP"
        O1["Python 3.x"]
        O2["C extensions"]
        O3["HDF5/PyTables"]
        O4["NumPy/SciPy"]
        O5["pysam"]
        O1 --> O2
        O2 --> O3
        O3 --> O4
        O4 --> O5
    end

    subgraph "WASP2"
        N1["Python 3.8+"]
        N2["Typer CLI"]
        N3["Rich terminal"]
        N4["Parquet/Arrow"]
        N5["AnnData"]
        N6["cyvcf2"]
        N1 --> N2
        N2 --> N3
        N1 --> N4
        N4 --> N5
        N1 --> N6
    end

    style O3 fill:#FFB6C1
    style N2 fill:#90EE90
    style N5 fill:#90EE90
```

## Single-Cell Integration

```mermaid
flowchart LR
    subgraph "WASP2 Output"
        COUNTS["Allele counts<br/>per cell × variant"]
    end

    subgraph "AnnData Structure"
        X["X: count matrix"]
        VAR["var: variant info"]
        OBS["obs: cell metadata"]
        LAYERS["layers: ref/alt counts"]
    end

    subgraph "Downstream"
        SCANPY["scanpy"]
        DIFF["Differential AI<br/>analysis"]
    end

    COUNTS --> X
    COUNTS --> VAR
    COUNTS --> LAYERS
    X --> SCANPY
    LAYERS --> DIFF

    style COUNTS fill:#87CEEB
    style SCANPY fill:#90EE90
    style DIFF fill:#FFD700
```

## Design Principles

```mermaid
mindmap
    root((WASP2))
        No Conversion
            VCF/BCF native
            tabix indexing
            No HDF5 step
        Unified CLI
            wasp2-count
            wasp2-map
            wasp2-analyze
        Modern Stack
            Typer
            Rich
            Parquet
        Single-Cell
            AnnData
            scanpy integration
            Millions of cells
```

---

## Episode Reference
- **Episode**: 002 - Building the New Hive
- **Topic**: WASP2 modernization and architecture (2021)
- **Repository**: mcvickerlab/WASP2
