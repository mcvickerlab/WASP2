# Buzz Report: Building the New Hive
# Episode: 002 | The WASP Chronicles
# Date: 2026-02-03

---

## Opening

Welcome to the Hive, fellow worker bees.

I'm the Queen Bee, and this is The WASP's Nest. Today we continue The WASP Chronicles with Episode Two... Building the New Hive.

In our last episode, we explored the original WASP from 2015... a groundbreaking tool that solved mapping bias. But by 2021, the field had evolved. Single-cell technologies exploded. VCF files became the universal standard. And a new generation of researchers needed modern tools.

This is the story of how WASP2 was born at the McVicker Lab.

---

## The Call to Rebuild

Let's set the scene. It's late 2021 at the Salk Institute. The original WASP is still widely used... but showing its age.

The pain points were real. Researchers had to convert every VCF file to HDF5 format before running any analysis. Single-cell experiments? Not supported. The command-line tools were scattered Python scripts with inconsistent interfaces. Dependencies were becoming harder to manage. And performance bottlenecks were slowing down large-scale studies.

Researchers were spending more time wrestling with file formats... than doing actual biology.

But there was opportunity. VCF and BCF had become universal standards. Single-cell ATAC-seq and RNA-seq were now mainstream. Modern Python packaging... with pyproject.toml, typer, and rich... had made CLI development elegant. The core algorithms were still sound. Only the interface needed modernization.

---

## Foraging: The New Design

Aaron Ho, working with the McVicker Lab, established a new repository... mcvickerlab WASP2. The vision was clear from day one.

The design principles were straightforward. First... no format conversion. Read VCF and BCF files directly. Eliminate the HDF5 step entirely. Second... a unified CLI. One tool with many subcommands, like git. Third... single-cell native support. First-class handling for scATAC and scRNA experiments. Fourth... modern packaging. A simple pip install. Clean dependencies. No headaches.

Here's what the transformation looked like in practice. The old way required multiple scripts... snp2h5 dot py to convert variants... find intersecting snps dot py to identify overlaps... filter remapped reads dot py for the filtering step. Multiple commands, multiple outputs, multiple opportunities for confusion.

The new way is elegantly simple. wasp2-count for counting alleles at variant sites. wasp2-map for the mapping bias correction pipeline. wasp2-analyze for detecting allelic imbalance. Clean. Intuitive. No HDF5 in sight.

---

## Building: The Architecture

The architects of WASP2 made thoughtful choices about the new hive's structure.

For the command-line interface, they chose Typer. Modern argument parsing with automatic help generation and shell completion. Each subcommand became a focused tool. wasp2-count handles allele counting at heterozygous variant sites. wasp2-map provides the unbiased read mapping pipeline. wasp2-analyze runs statistical analysis for detecting allelic imbalance. And wasp2-ipscore enables QTL scoring workflows.

For terminal output, they integrated Rich. Beautiful progress bars, colored output, and informative error messages. No more walls of text flooding the terminal.

For single-cell support, they built native AnnData integration. The scanpy ecosystem's data structure became a first-class citizen. Single-cell researchers could take WASP2 output and flow directly into downstream analysis.

The module organization reflects this clarity. The counting module handles allele counting at heterozygous sites. The mapping module manages the read filtering pipeline. The analysis module contains the statistical models... specifically the beta-binomial distribution for detecting allelic imbalance. And the I/O module supports VCF, BCF, and even the high-performance PGEN format.

Pure Python... cleanly organized... well-documented.

---

## Defending: The Statistical Heart

One thing WASP2 never compromised on... the core science.

The mapping bias correction strategy remained unchanged from the original. Find reads overlapping heterozygous variants. Swap the alleles in the read sequence. Remap both versions. Filter out any reads that map differently. Simple. Principled. Effective.

But the statistical analysis evolved. While the original WASP used the Combined Haplotype Test... WASP2 took a different approach. The new analysis module centers on the beta-binomial distribution.

Here's why this matters. When you count alleles at a heterozygous site, you expect roughly fifty-fifty between reference and alternate. But biological and technical variation create overdispersion... more variance than a simple binomial would predict. The beta-binomial model captures this elegantly with two parameters. Mu represents the mean imbalance probability. Rho captures the dispersion.

WASP2 fits these parameters using likelihood optimization, then runs a likelihood ratio test. The null hypothesis... no allelic imbalance, mu equals 0.5. The alternative... imbalance exists. The test statistic follows a chi-squared distribution... giving you a p-value you can trust.

The model supports both phased and unphased genotypes. For phased data, the optimization is direct. For unphased data, a clever dynamic programming approach averages over possible phase configurations.

This is the scientific heart of WASP2. Robust statistical testing... properly accounting for overdispersion... with principled inference.

---

## Deep Dive: VCF Native

For the technically curious bees... let's explore the VCF handling innovation.

The original WASP used HDF5 because random access to variants was critical. You need to quickly look up which variants overlap each read. HDF5 provided indexed arrays for this.

WASP2 solved this problem differently. VCF indexing via tabix provides genomic coordinate indexing through the tbi files. Pysam's TabixFile class enables fast region queries without any format conversion. And for maximum speed, the cyvcf2 backend offers C-accelerated VCF parsing... roughly seven times faster than pure Python.

But WASP2 went further. Beyond VCF, the BCF format... the binary version of VCF... offers another seven-fold speedup through native binary parsing. And for the ultimate performance, PGEN format support via Pgenlib delivers a stunning twenty-five times speedup over standard VCF.

Users can keep their existing files... no conversion pipeline required. Just choose the format that matches your performance needs.

---

## Pollinating: The Ecosystem

WASP2 was designed to play nicely with the broader bioinformatics ecosystem.

For inputs... BAM or CRAM files from any aligner. VCF, BCF, or PGEN from any variant caller or imputation pipeline. Standard FASTQ for the remapping step.

For outputs... TSV files for simple downstream processing. Parquet for efficient columnar storage and fast queries. And AnnData in H5AD format for seamless single-cell integration.

The interoperability is deliberate. Standard bcftools and samtools compatibility. Integration with the scanpy and AnnData ecosystem. Bioconda packaging for easy installation.

WASP2 didn't reinvent wheels... it connected them.

---

## The Timeline

The journey from concept to release tells a story of steady progress.

December 2021... the repository was established. Through 2022... the core counting and mapping modules took shape. In 2023... single-cell support arrived alongside robust testing infrastructure. September 2024 marked the v1.0.0 official release. November 2024 brought v1.1.0... and the beginning of Rust acceleration.

That performance revolution... that's a story for our next episode.

---

## Closing

And that's the buzz on building the new hive, worker bees.

WASP2 represented a modern reimagining of the original vision. Same proven science for mapping bias correction. New accessible interface for modern workflows. The McVicker Lab took a decade of lessons learned and built something that feels native to 2020s research.

The key insights from this chapter... Modernization doesn't mean reinvention. The core science remained. Developer experience matters... unified CLI, no format conversion, clean outputs. And ecosystem integration accelerates adoption.

In our next episode... we'll witness the Rust metamorphosis. When WASP2 learned to fly at lightning speed.

Keep building... keep buzzing. May your reads map true and your alleles balance.

From the WASP's Nest... this is the Queen Bee.

Buzz out.

---

## Episode Metadata

```yaml
episode:
  number: 2
  title: "Building the New Hive"
  subtitle: "McVicker Lab WASP2"
  series: "The WASP Chronicles"
  date: "2026-02-03"
  duration_estimate: "10-12 minutes"
  source_repo: "https://github.com/mcvickerlab/WASP2"
  authors:
    - "Aaron Ho - Creator of WASP2"
    - "Jeff Jaureguy - Developer and maintainer"
    - "McVicker Lab, Salk Institute"
  timeline:
    established: "2021-12"
    v1_release: "2024-09"
    v1_1_release: "2024-11"
  technical_highlights:
    - "Beta-binomial model for allelic imbalance (NOT CHT)"
    - "VCF/BCF/PGEN native support (no HDF5)"
    - "Single-cell via AnnData/H5AD"
    - "Unified CLI: wasp2-count, wasp2-map, wasp2-analyze, wasp2-ipscore"
  chapters:
    - name: "The Call"
      topics: ["modernization", "pain points", "opportunity"]
    - name: "Foraging"
      topics: ["design principles", "unified CLI", "no HDF5"]
    - name: "Building"
      topics: ["Typer", "Rich", "AnnData", "module organization"]
    - name: "Defending"
      topics: ["beta-binomial model", "likelihood ratio test", "phased/unphased"]
    - name: "Deep Dive"
      topics: ["VCF native", "BCF 7x", "PGEN 25x", "pysam", "cyvcf2"]
    - name: "Pollinating"
      topics: ["ecosystem integration", "format support", "AnnData output"]
```
