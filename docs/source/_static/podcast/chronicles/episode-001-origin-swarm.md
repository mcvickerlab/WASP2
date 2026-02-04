# Buzz Report: The Origin Swarm
# Episode: 001 | The WASP Chronicles
# Date: 2026-02-03

---

## Opening

[happy buzz]

Welcome to the Hive, fellow worker bees!

I'm the Queen Bee, and this is The WASP's Nest - bringing you something special today. Instead of our usual release notes, we're going back to the beginning. This is Episode One of The WASP Chronicles, where we trace the lineage of our hive.

Today's Buzz Report takes us back to 2015, when the first WASP was born.

---

## The Problem: Mapping Bias

[serious tone]

Picture this, worker bees: you're a researcher trying to understand which version of a gene is more active. You sequence RNA from cells, map those reads to the genome, and count how many come from each allele.

Simple, right? **Wrong.**

Here's the sting: reads carrying the reference allele map differently than reads carrying the alternate allele. If your read has a variant that doesn't match the reference genome, the aligner might:
- Map it to the wrong place
- Give it a lower quality score
- Fail to map it entirely

This creates **systematic bias** toward the reference allele. And when you're looking for allele-specific expression? That bias looks exactly like the biological signal you're hunting for.

[frustrated buzz]

False positives everywhere! Real signals getting buried!

---

## The Foraging: A Clever Solution

[excited waggle]

In 2015, a team of brilliant researchers at Stanford and the University of Chicago forged a solution. Bryce van de Geijn, Graham McVicker, Yoav Gilad, and Jonathan Pritchard published their landmark paper in *Nature Methods*:

**"WASP: allele-specific software for robust molecular quantitative trait locus discovery"**

Their approach was elegantly simple:

### The WASP Read Filtering Strategy

1. **Find reads overlapping variants** - identify which reads touch heterozygous sites
2. **Swap the alleles** - create an alternate version of each read with the other allele
3. **Remap both versions** - send both through the aligner
4. **Filter discordant reads** - if they don't map to the same place with the same quality, throw them out

[impressed tone]

The genius? Any read that maps differently depending on which allele it carries is biased by definition. By removing these reads, you eliminate the bias at its source!

---

## Building: The Combined Haplotype Test

[precise tone]

But wait - there's more! The original WASP didn't just fix mapping bias. It introduced a powerful statistical test called the **Combined Haplotype Test (CHT)**.

Traditional approaches tested either:
- **Read depth** - does a genetic variant affect total expression?
- **Allelic imbalance** - among heterozygotes, is one allele more expressed?

The CHT combined both signals into a single test:

[technical tone]

The test integrates across individuals, combining:
- Total read counts at the gene level
- Allele-specific read counts at heterozygous sites within the gene
- Proper handling of overdispersion using a beta-binomial model

This gave substantially more power to detect expression QTLs than either approach alone.

---

## The Original Architecture

[contemplative buzz]

The 2015 WASP was built for its era:

**Technology Stack:**
- Python 3.x with C extensions (77.7% Python, 19.5% C)
- HDF5 format for variant storage (via PyTables)
- NumPy and SciPy for numerical computation
- pysam for BAM file handling

**The Tools:**
- `snp2h5` - Convert VCF files to HDF5 format
- `find_intersecting_snps.py` - Find reads overlapping variants
- `filter_remapped_reads.py` - Remove biased reads after remapping
- `combined_test.py` - Run the CHT for QTL discovery

[nostalgic tone]

The HDF5 requirement was pragmatic for 2015 - it offered fast random access to millions of variants. But it also meant users had to convert their VCF files before running the pipeline.

---

## Deep Dive: The Science

[scholarly tone]

For the bioinformaticians in the hive, let's go deeper.

The key insight was modeling read mapping as a stochastic process. Given a heterozygous site with alleles A and B, a read carrying allele A might have mapping probability P_A, while the same read with allele B has probability P_B.

If P_A != P_B, that read is biased. By simulating the alternate allele and testing empirically, WASP avoided the need to model aligner behavior analytically.

The CHT used a likelihood ratio test:
- Null: no genetic effect (expression independent of genotype)
- Alternative: genetic effect present (QTL exists)

The test statistic followed a chi-squared distribution under the null, with overdispersion handled by the beta-binomial model for allelic counts.

---

## Illumination

See: `illuminations/illumination-001-wasp-mapping.md` for the visual diagram of the WASP mapping bias correction strategy.

---

## The Impact

[proud buzz]

The original WASP made a lasting mark:

- **529 commits** over 4+ years of development
- **111 stars** on GitHub (github.com/bmvdgeijn/WASP)
- **Last release**: v0.3.4 (April 2019)
- **Cited by**: hundreds of eQTL and ASE studies

But perhaps most importantly, it established the fundamental approach that all subsequent allele-specific analysis tools would build upon.

---

## Closing

[pause]

And that's the buzz on where it all began, worker bees!

The original WASP showed us that mapping bias isn't just a nuisance - it's a fundamental problem that requires a principled solution. By swapping alleles and filtering discordant reads, van de Geijn and colleagues gave the field a tool that remains influential a decade later.

Remember:
- Mapping bias is real and can masquerade as biological signal
- The WASP filtering strategy removes bias at its source
- Combining read depth and allelic imbalance increases power

In our next episode, we'll see how the McVicker Lab took these foundational ideas and built something new.

Keep building, keep buzzing!
May your reads map true and your alleles balance.

From the WASP's Nest, this is the Queen Bee.

Buzz out!

---

## Episode Metadata

```yaml
episode:
  number: 1
  title: "The Origin Swarm"
  subtitle: "Original WASP (2015)"
  series: "The WASP Chronicles"
  date: "2026-02-03"
  duration_estimate: "8-10 minutes"
  source_paper:
    title: "WASP: allele-specific software for robust molecular quantitative trait locus discovery"
    authors: ["van de Geijn B", "McVicker G", "Gilad Y", "Pritchard JK"]
    journal: "Nature Methods"
    year: 2015
    pmid: 26366987
    doi: "10.1038/nmeth.3582"
  source_repo: "https://github.com/bmvdgeijn/WASP"
  chapters:
    - name: "The Problem"
      topics: ["mapping bias", "allele-specific analysis", "false positives"]
    - name: "Foraging"
      topics: ["WASP filtering", "allele swapping", "read remapping"]
    - name: "Building"
      topics: ["Combined Haplotype Test", "beta-binomial", "QTL detection"]
    - name: "Deep Dive"
      topics: ["statistical model", "likelihood ratio test"]
    - name: "Impact"
      topics: ["citations", "field influence"]
```
