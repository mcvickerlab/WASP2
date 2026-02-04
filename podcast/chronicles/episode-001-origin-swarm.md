# Buzz Report: The Origin Swarm
# Episode: 001 | The WASP Chronicles
# Date: 2026-02-03

---

## Opening

Welcome to the Hive, fellow worker bees.

I'm the Queen Bee, and this is The WASP's Nest. Today we're bringing you something special. Instead of our usual release notes, we're going back to the beginning. This is Episode One of The WASP Chronicles... where we trace the lineage of our hive.

Today's Buzz Report takes us back to 2015... when the first WASP was born.

---

## The Problem: Mapping Bias

Picture this, worker bees. You're a researcher trying to understand which version of a gene is more active. You sequence RNA from cells, map those reads to the genome, and count how many come from each allele.

Simple, right?... Wrong.

Here's the sting. Reads carrying the reference allele map differently than reads carrying the alternate allele. If your read has a variant that doesn't match the reference genome, the aligner might map it to the wrong place... give it a lower quality score... or fail to map it entirely.

This creates systematic bias toward the reference allele. And when you're looking for allele-specific expression?... That bias looks exactly like the biological signal you're hunting for.

False positives everywhere. Real signals getting buried.

---

## The Foraging: A Clever Solution

In 2015, a team of brilliant researchers at Stanford and the University of Chicago forged a solution. Bryce van de Geijn, Graham McVicker, Yoav Gilad, and Jonathan Pritchard published their landmark paper in Nature Methods.

The title... "WASP: allele-specific software for robust molecular quantitative trait locus discovery."

Their approach was elegantly simple. The WASP Read Filtering Strategy works in four steps.

First... find reads overlapping variants. Identify which reads touch heterozygous sites.

Second... swap the alleles. Create an alternate version of each read with the other allele.

Third... remap both versions. Send both through the aligner.

Fourth... filter discordant reads. If they don't map to the same place with the same quality... throw them out.

The genius of this approach is clear. Any read that maps differently depending on which allele it carries is biased by definition. By removing these reads... you eliminate the bias at its source.

---

## Building: The Combined Haplotype Test

But wait... there's more. The original WASP didn't just fix mapping bias. It introduced a powerful statistical test called the Combined Haplotype Test... or CHT.

Traditional approaches tested either read depth... does a genetic variant affect total expression?... or allelic imbalance... among heterozygotes, is one allele more expressed?

The CHT combined both signals into a single test.

The test integrates across individuals, combining total read counts at the gene level... allele-specific read counts at heterozygous sites within the gene... and proper handling of overdispersion using a beta-binomial model.

This gave substantially more power to detect expression QTLs than either approach alone.

---

## The Original Architecture

The 2015 WASP was built for its era.

The technology stack included Python 3.x with C extensions... about 77 percent Python and 19 percent C. HDF5 format for variant storage via PyTables. NumPy and SciPy for numerical computation. And pysam for BAM file handling.

The tools were straightforward. snp2h5 converted VCF files to HDF5 format. find_intersecting_snps.py found reads overlapping variants. filter_remapped_reads.py removed biased reads after remapping. And combined_test.py ran the CHT for QTL discovery.

The HDF5 requirement was pragmatic for 2015... it offered fast random access to millions of variants. But it also meant users had to convert their VCF files before running the pipeline.

---

## Deep Dive: The Science

For the bioinformaticians in the hive... let's go deeper.

The key insight was modeling read mapping as a stochastic process. Given a heterozygous site with alleles A and B, a read carrying allele A might have mapping probability P_A... while the same read with allele B has probability P_B.

If P_A is not equal to P_B... that read is biased. By simulating the alternate allele and testing empirically, WASP avoided the need to model aligner behavior analytically.

The CHT used a likelihood ratio test. The null hypothesis states no genetic effect... expression is independent of genotype. The alternative hypothesis states a genetic effect is present... a QTL exists.

The test statistic follows a chi-squared distribution under the null... with overdispersion handled by the beta-binomial model for allelic counts.

---

## The Impact

The original WASP made a lasting mark.

529 commits over four-plus years of development. 111 stars on GitHub at github.com slash bmvdgeijn slash WASP. Last release v0.3.4 in April 2019. And cited by hundreds of eQTL and ASE studies worldwide.

But perhaps most importantly... it established the fundamental approach that all subsequent allele-specific analysis tools would build upon.

---

## Closing

And that's the buzz on where it all began, worker bees.

The original WASP showed us that mapping bias isn't just a nuisance... it's a fundamental problem that requires a principled solution. By swapping alleles and filtering discordant reads, van de Geijn and colleagues gave the field a tool that remains influential a decade later.

The key takeaways from this episode. Mapping bias is real and can masquerade as biological signal. The WASP filtering strategy removes bias at its source. And combining read depth and allelic imbalance increases statistical power.

In our next episode... we'll see how the McVicker Lab took these foundational ideas and built something new.

Keep building... keep buzzing. May your reads map true and your alleles balance.

From the WASP's Nest... this is the Queen Bee.

Buzz out.

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
  note: "The original WASP used the Combined Haplotype Test (CHT). WASP2 replaced CHT with a beta-binomial model for allelic imbalance detection."
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
