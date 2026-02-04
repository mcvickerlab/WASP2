# Illumination: WASP Mapping Bias Correction
# Episode: 001 - The Origin Swarm

## The Problem: Mapping Bias

When reads contain genetic variants, they may map differently depending on which allele they carry.

```mermaid
graph TD
    subgraph "The Bias Problem"
        R1["Read with REF allele<br/>ACGT<b>A</b>CGT"]
        R2["Read with ALT allele<br/>ACGT<b>G</b>CGT"]

        R1 -->|"Maps perfectly"| M1["✓ High MAPQ<br/>Correct position"]
        R2 -->|"Mismatch penalty"| M2["✗ Lower MAPQ<br/>May mismap or fail"]
    end

    style M1 fill:#90EE90
    style M2 fill:#FFB6C1
```

## The WASP Solution: Allele Swap & Filter

```mermaid
flowchart LR
    subgraph "Step 1: Find Overlapping Reads"
        BAM["BAM File"] --> FIND["Find reads<br/>at het sites"]
        VCF["VCF File"] --> FIND
    end

    subgraph "Step 2: Create Alternate Reads"
        FIND --> ORIG["Original Read<br/>ACGT<b>A</b>CGT"]
        ORIG --> SWAP["Swap allele"]
        SWAP --> ALT["Alternate Read<br/>ACGT<b>G</b>CGT"]
    end

    subgraph "Step 3: Remap Both"
        ORIG --> ALIGN1["Align"]
        ALT --> ALIGN2["Align"]
        ALIGN1 --> POS1["Position 1<br/>MAPQ 60"]
        ALIGN2 --> POS2["Position 2<br/>MAPQ 55"]
    end

    subgraph "Step 4: Compare & Filter"
        POS1 --> COMP{"Same position<br/>& quality?"}
        POS2 --> COMP
        COMP -->|"Yes"| KEEP["✓ KEEP<br/>Unbiased read"]
        COMP -->|"No"| DISCARD["✗ DISCARD<br/>Biased read"]
    end

    style KEEP fill:#90EE90
    style DISCARD fill:#FFB6C1
```

## The Combined Haplotype Test (CHT)

```mermaid
graph TB
    subgraph "Two Sources of Signal"
        RD["Read Depth Signal<br/>(across all individuals)"]
        AI["Allelic Imbalance Signal<br/>(within heterozygotes)"]
    end

    subgraph "Combined Haplotype Test"
        RD --> CHT["Integrate both signals<br/>in likelihood framework"]
        AI --> CHT
        CHT --> BB["Beta-binomial model<br/>handles overdispersion"]
        BB --> LRT["Likelihood Ratio Test"]
        LRT --> PVAL["p-value for QTL"]
    end

    style CHT fill:#87CEEB
    style PVAL fill:#FFD700
```

## The Original WASP Pipeline

```mermaid
flowchart TB
    subgraph "Input"
        VCF["VCF files"]
        BAM["BAM files"]
    end

    subgraph "Preparation"
        VCF --> SNP2H5["snp2h5<br/>(convert to HDF5)"]
        SNP2H5 --> H5["HDF5 database"]
    end

    subgraph "WASP Filtering"
        BAM --> FIND["find_intersecting_snps.py"]
        H5 --> FIND
        FIND --> REMAP["Remap with alternate alleles"]
        REMAP --> FILTER["filter_remapped_reads.py"]
        FILTER --> CLEAN["Filtered BAM<br/>(bias removed)"]
    end

    subgraph "Analysis"
        CLEAN --> COUNT["Count alleles"]
        COUNT --> CHT2["combined_test.py<br/>(CHT)"]
        CHT2 --> QTL["QTL Results"]
    end

    style H5 fill:#FFA07A
    style CLEAN fill:#90EE90
    style QTL fill:#FFD700
```

## Key Insight

```mermaid
graph LR
    A["If a read maps differently<br/>depending on which allele<br/>it carries..."] --> B["...that read is<br/>BIASED<br/>by definition"]
    B --> C["Remove it!"]

    style A fill:#FFB6C1
    style B fill:#FF6347
    style C fill:#90EE90
```

---

## Episode Reference
- **Episode**: 001 - The Origin Swarm
- **Topic**: Original WASP mapping bias correction (2015)
- **Paper**: van de Geijn et al., Nature Methods 2015
