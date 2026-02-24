# WASP2 Container Usage Guide

This guide covers how to use WASP2 containers for local development, HPC clusters, and Nextflow pipelines.

## Container Registries

WASP2 images are available from:

| Registry | Image | Pull Command |
|----------|-------|--------------|
| **DockerHub** | `mcvickerlab/wasp2` | `docker pull mcvickerlab/wasp2:latest` |
| **GitHub Container Registry** | `ghcr.io/mcvickerlab/wasp2` | `docker pull ghcr.io/mcvickerlab/wasp2:latest` |

### Available Tags

- `:latest` - Most recent release
- `:1.3.0` - Specific version
- `:1.3` - Minor version (tracks patches)
- `:main` - Development builds from main branch

## Docker Usage

### Pull and Run

```bash
# Pull the image
docker pull mcvickerlab/wasp2:latest

# Run WASP2 commands
docker run --rm mcvickerlab/wasp2 wasp2-count --help
docker run --rm mcvickerlab/wasp2 wasp2-map --help
docker run --rm mcvickerlab/wasp2 wasp2-analyze --help

# Process files (mount local directory)
docker run --rm -v $(pwd):/data mcvickerlab/wasp2 \
    wasp2-count /data/sample.bam /data/variants.vcf.gz -o /data/counts.tsv
```

### Interactive Shell

```bash
docker run -it --rm -v $(pwd):/data mcvickerlab/wasp2 /bin/bash
```

## Singularity/Apptainer Usage (HPC)

### Pull from Docker Registry

```bash
# Pull and convert to SIF
singularity pull wasp2.sif docker://mcvickerlab/wasp2:latest

# Or from GHCR
singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:latest
```

### Build from Definition File

```bash
# Clone the repository
git clone https://github.com/mcvickerlab/WASP2.git
cd WASP2

# Build the container
singularity build wasp2.sif Singularity.def
```

### Run Commands

```bash
# Run WASP2 commands
singularity exec wasp2.sif wasp2-count --help

# Process files (current directory is auto-bound)
singularity exec wasp2.sif wasp2-count sample.bam variants.vcf.gz -o counts.tsv

# With explicit bindings
singularity exec --bind /scratch:/scratch wasp2.sif \
    wasp2-map make-reads /scratch/input.bam /scratch/variants.vcf
```

### SLURM Job Script Example

```bash
#!/bin/bash
#SBATCH --job-name=wasp2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

module load singularity

CONTAINER=/path/to/wasp2.sif

# Count variants (wasp2-count does not have --threads option)
singularity exec ${CONTAINER} wasp2-count \
    input.bam \
    variants.vcf.gz \
    -o counts.tsv

# WASP mapping filter (supports --threads)
singularity exec ${CONTAINER} wasp2-map make-reads \
    input.bam \
    variants.vcf.gz \
    --threads ${SLURM_CPUS_PER_TASK} \
    --out_dir ./wasp_output
```

## Nextflow Integration

### Configuration

Add to your `nextflow.config`:

```groovy
profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        singularity.cacheDir = "${HOME}/.singularity/cache"
    }
}

process {
    withLabel: 'wasp2' {
        container = 'mcvickerlab/wasp2:latest'
    }
}
```

### Running Pipelines

```bash
# With Docker
nextflow run main.nf -profile docker

# With Singularity
nextflow run main.nf -profile singularity
```

## Building Locally

### Docker Build

```bash
# Clone repository
git clone https://github.com/mcvickerlab/WASP2.git
cd WASP2

# Build image
docker build -t wasp2:local .

# Test the build
docker run --rm wasp2:local wasp2-count --version
docker run --rm wasp2:local python -c "import wasp2_rust; print('OK')"
```

### Manual Build (for maintainers)

Note: Currently only `linux/amd64` is supported.

```bash
# Set up buildx
docker buildx create --name wasp2builder --use

# Build with version argument
docker buildx build \
    --platform linux/amd64 \
    --build-arg VERSION=1.3.0 \
    -t mcvickerlab/wasp2:1.3.0 \
    -t mcvickerlab/wasp2:latest \
    --push .
```

## Container Contents

The WASP2 container includes:

### Python Environment
- Python 3.10+ (container uses 3.11)
- wasp2 package with Rust extension
- Core: pysam, pandas (<2.0), numpy, scipy, polars
- CLI: typer, rich
- Single-cell: anndata, scanpy (optional)

### Rust Components
- Pre-built `wasp2_rust` Python extension
- Compiled with release optimizations

### CLI Tools

Each tool has subcommands for different analysis modes:

- **`wasp2-count`** - Allele counting
  - `count-variants` - Bulk allele counting at heterozygous sites (default)
  - `count-variants-sc` - Single-cell allele counting

- **`wasp2-map`** - WASP mapping filter
  - `make-reads` - Generate reads with swapped alleles for remapping
  - `filter-remapped` - Filter remapped reads using WASP algorithm

- **`wasp2-analyze`** - Statistical analysis
  - `find-imbalance` - Calculate allelic imbalance
  - `find-imbalance-sc` - Single-cell allelic imbalance analysis
  - `compare-imbalance` - Compare imbalance between cell types/groups

### Bioinformatics Tools
- samtools
- bcftools
- bedtools
- tabix

## Troubleshooting

### Permission Issues (Docker)

```bash
# Run as current user
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data mcvickerlab/wasp2 ...
```

### Cache Issues (Singularity)

```bash
# Clear Singularity cache
singularity cache clean

# Use different cache directory
export SINGULARITY_CACHEDIR=/scratch/singularity_cache
```

### Verify Installation

```bash
# Docker
docker run --rm mcvickerlab/wasp2 wasp2-count --version
docker run --rm mcvickerlab/wasp2 python -c "import wasp2_rust; print('Rust extension OK')"

# Singularity
singularity exec wasp2.sif wasp2-count --version
singularity exec wasp2.sif python -c "import wasp2_rust; print('Rust extension OK')"
```

## GitHub Actions Secrets Setup

To enable automated container builds, repository maintainers must configure:

1. **DockerHub Secrets** (Settings → Secrets and variables → Actions):
   - `DOCKERHUB_USERNAME`: Your DockerHub username
   - `DOCKERHUB_TOKEN`: DockerHub access token (Account Settings → Security → Access Tokens)

2. **GitHub Container Registry**: Uses `GITHUB_TOKEN` automatically (no setup needed)

## Related Documentation

- [Nextflow Pipelines](../pipelines/nf-atacseq/README.md)
- [WASP2 Ecosystem](WASP2_ECOSYSTEM.md)
- [GitHub Repository](https://github.com/mcvickerlab/WASP2)
