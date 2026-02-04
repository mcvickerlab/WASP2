# Buzz Report: Building the New Hive
# Episode: 002 | The WASP Chronicles
# Date: 2026-02-03

---

## Opening

[happy buzz]

Welcome to the Hive, fellow worker bees!

I'm the Queen Bee, and this is The WASP's Nest. Today we continue The WASP Chronicles with Episode Two: Building the New Hive.

In our last episode, we explored the original WASP from 2015 - a groundbreaking tool that solved mapping bias. But by 2021, the field had evolved. Single-cell technologies exploded. VCF files became the standard. And a new generation of researchers needed modern tools.

This is the story of WASP2's birth at the McVicker Lab.

---

## The Call to Rebuild

[contemplative buzz]

Let's set the scene. It's late 2021. The original WASP is still widely used, but showing its age:

**The Pain Points:**
- HDF5 conversion required before every analysis
- No support for single-cell data
- CLI tools were scattered Python scripts
- Dependencies becoming harder to manage
- Performance bottlenecks in large-scale studies

[frustrated tone]

Researchers were spending more time wrestling with file formats than doing biology!

**The Opportunity:**
- VCF/BCF had become universal
- Single-cell ATAC-seq and RNA-seq were mainstream
- Modern Python packaging (pyproject.toml, typer, rich) made CLI development elegant
- The core algorithms were still sound - just the interface needed modernization

---

## Foraging: The New Design

[excited waggle]

In December 2021, the McVicker Lab established a new repository: `mcvickerlab/WASP2`. The vision was clear:

### Design Principles

1. **No format conversion** - Read VCF/BCF directly, no HDF5 step
2. **Unified CLI** - One tool, many subcommands (like `git`)
3. **Single-cell native** - First-class support for scATAC and scRNA
4. **Modern packaging** - pip-installable, clean dependencies

### The New Interface

```bash
# The old way (multiple scripts)
python snp2h5.py variants.vcf snps.h5
python find_intersecting_snps.py reads.bam snps.h5 output/
python filter_remapped_reads.py ...

# The new way (unified CLI)
wasp2-count reads.bam variants.vcf.gz -o counts.parquet
wasp2-map reads.bam variants.vcf.gz -o filtered/
wasp2-analyze counts.parquet -o results/
```

[satisfied tone]

Clean. Intuitive. No HDF5 in sight!

---

## Building: The Architecture

[precise tone]

The architects of WASP2 made thoughtful decisions about the new hive's structure:

### Technology Choices

**Typer for CLI** - Modern argument parsing with automatic help generation and shell completion. Each subcommand became a focused tool:
- `wasp2-count` - Allele counting at variant sites
- `wasp2-map` - Unbiased read mapping pipeline
- `wasp2-analyze` - Statistical analysis for QTL discovery

**Rich for Terminal Output** - Beautiful progress bars, colored output, and informative error messages. No more wall-of-text logs!

**AnnData Integration** - The `scanpy` ecosystem's data structure became native. Single-cell researchers could go from WASP2 output directly to downstream analysis.

### Module Organization

```
src/wasp2/
├── counting/        # Allele counting at het sites
├── mapping/         # Read filtering pipeline
├── analysis/        # Statistical tests (CHT, binomial)
├── io/              # VCF, BAM, Parquet readers
└── cli/             # typer-based command line
```

[thoughtful tone]

Pure Python, cleanly organized, well-documented.

---

## Defending: Staying True

[serious tone]

One thing WASP2 never compromised on: the core science.

The mapping bias correction strategy remained unchanged:
1. Find reads overlapping heterozygous variants
2. Swap alleles in the read sequence
3. Remap and compare
4. Filter reads that map differently

The Combined Haplotype Test stayed at the heart of QTL discovery.

[emphatic buzz]

The algorithms that made WASP powerful in 2015? Still powerful in 2021. WASP2 just made them accessible.

---

## Deep Dive: VCF Native

[technical tone]

For the technically curious bees, let's explore the VCF handling.

The original WASP used HDF5 because random access to variants was critical - you need to quickly look up "what variants overlap this read?" HDF5 provided indexed arrays.

WASP2 solved this differently:
1. **VCF indexing via tabix** - `.vcf.gz.tbi` files provide genomic coordinate indexing
2. **pysam's TabixFile** - Fast region queries without format conversion
3. **cyvcf2 backend** - C-accelerated VCF parsing when needed

```python
# Query variants overlapping a read's alignment
with pysam.TabixFile(vcf_path) as vcf:
    for record in vcf.fetch(chrom, start, end):
        # Process variant record
```

This approach meant users could keep their existing VCF files - no conversion pipeline required.

---

## Pollinating: The Ecosystem

[playful buzz]

WASP2 was designed to play nicely with the broader bioinformatics ecosystem:

**Input Formats:**
- BAM/CRAM from any aligner
- VCF/BCF from any variant caller
- Standard FASTQ for remapping

**Output Formats:**
- Parquet for efficient columnar storage
- AnnData for single-cell integration
- TSV for simple downstream processing

**Interoperability:**
- bcftools, samtools compatibility
- scanpy/AnnData ecosystem
- Standard Bioconda packaging

[collaborative tone]

WASP2 didn't reinvent wheels - it connected them.

---

## The Timeline

[narrative tone]

The journey from concept to release:

**December 2021** - Repository established
**2022** - Core counting and mapping modules developed
**2023** - Single-cell support, testing infrastructure
**September 2024** - v1.0.0 official release
**November 2024** - v1.1.0 with Rust acceleration beginning

[pause]

The next chapter would bring a performance revolution. But that's a story for our next episode.

---

## Illumination

See: `illuminations/illumination-002-architecture.md` for the WASP2 architecture diagram showing the module organization and data flow.

---

## Closing

[pause]

And that's the buzz on building the new hive, worker bees!

WASP2 represented a modern reimagining of the original vision. Same proven science, new accessible interface. The McVicker Lab took a decade of lessons learned and built something that felt native to the 2020s research workflow.

Remember:
- Modernization doesn't mean reinvention - the core science remained
- Developer experience matters - unified CLI, no format conversion
- Ecosystem integration accelerates adoption

In our next episode, we'll witness the Rust metamorphosis - when WASP2 learned to fly at lightning speed.

Keep building, keep buzzing!
May your reads map true and your alleles balance.

From the WASP's Nest, this is the Queen Bee.

Buzz out!

---

## Episode Metadata

```yaml
episode:
  number: 2
  title: "Building the New Hive"
  subtitle: "McVicker Lab WASP2"
  series: "The WASP Chronicles"
  date: "2026-02-03"
  duration_estimate: "8-10 minutes"
  source_repo: "https://github.com/mcvickerlab/WASP2"
  timeline:
    established: "2021-12"
    v1_release: "2024-09"
    v1_1_release: "2024-11"
  chapters:
    - name: "The Call"
      topics: ["modernization", "pain points", "opportunity"]
    - name: "Foraging"
      topics: ["design principles", "unified CLI", "no HDF5"]
    - name: "Building"
      topics: ["Typer", "Rich", "AnnData", "module organization"]
    - name: "Defending"
      topics: ["core science preserved", "CHT", "mapping bias"]
    - name: "Deep Dive"
      topics: ["VCF native", "tabix", "pysam"]
    - name: "Pollinating"
      topics: ["ecosystem integration", "format support"]
```
