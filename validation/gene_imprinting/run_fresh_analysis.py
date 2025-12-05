#!/usr/bin/env python3
"""Fresh WASP2-Rust gene imprinting analysis with INDEL support."""
import sys
import gzip
sys.path.insert(0, '/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp')

from wasp2_rust import BamCounter
from datetime import datetime

BAM = "/iblm/netapp/data3/aho/alignment/GM12878_rna_v2/GM12878_merged.sorted.bam"
VCF = "/iblm/netapp/data3/jjaureguy/gvl_files/wasp2/WASP2_extensive_evaluation/WASP2_current/cvpc/WASP2-exp/validation/gene_imprinting/results/imprinted_het_variants.vcf.gz"

print("=" * 60)
print("WASP2-Rust Gene Imprinting Analysis - FRESH RUN")
print(f"Time: {datetime.now()}")
print("=" * 60)
print(f"\nBAM: {BAM}")
print(f"VCF: {VCF}")

# Parse VCF
regions = []
variant_info = []

with gzip.open(VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]

        # Format: (chrom, pos, ref, alt) - 1-based position
        regions.append((chrom, pos, ref, alt))
        is_indel = len(ref) != len(alt)
        variant_info.append({
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'is_indel': is_indel
        })

print(f"\nTotal variants: {len(regions)}")
snp_count = sum(1 for v in variant_info if not v['is_indel'])
indel_count = sum(1 for v in variant_info if v['is_indel'])
print(f"SNPs: {snp_count}")
print(f"INDELs: {indel_count}")

# Run BamCounter
print("\nRunning WASP2-Rust BamCounter...")
counter = BamCounter(BAM)
counts = counter.count_alleles(regions, min_qual=20, threads=4)
print("Done!")

# Analyze results
snp_with_reads = []
indel_with_reads = []

for i, (ref_cnt, alt_cnt, other_cnt) in enumerate(counts):
    total = ref_cnt + alt_cnt
    if total > 0:
        if variant_info[i]['is_indel']:
            indel_with_reads.append((variant_info[i], ref_cnt, alt_cnt))
        else:
            snp_with_reads.append((variant_info[i], ref_cnt, alt_cnt))

print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)
print(f"SNPs with reads: {len(snp_with_reads)}")
print(f"INDELs with reads: {len(indel_with_reads)}")

snp_total = sum(r + a for _, r, a in snp_with_reads)
indel_total = sum(r + a for _, r, a in indel_with_reads)
total = snp_total + indel_total

print(f"\nTotal reads at SNPs: {snp_total}")
print(f"Total reads at INDELs: {indel_total}")
print(f"Grand total: {total}")

if total > 0:
    print(f"\n*** INDEL contribution: {100*indel_total/total:.1f}% of reads ***")
    print(f"*** INDEL variants: {len(indel_with_reads)} sites ({100*len(indel_with_reads)/(len(snp_with_reads)+len(indel_with_reads)):.1f}% of covered variants) ***")

print("\n" + "=" * 60)
print("INDEL DETAILS")
print("=" * 60)
for info, ref_cnt, alt_cnt in indel_with_reads:
    ratio = ref_cnt / (ref_cnt + alt_cnt) if (ref_cnt + alt_cnt) > 0 else 0.5
    print(f"  {info['chrom']}:{info['pos']} {info['ref']}>{info['alt']}: ref={ref_cnt}, alt={alt_cnt}, ratio={ratio:.2f}")

print("\n" + "=" * 60)
print(f"COMPLETED: {datetime.now()}")
print("=" * 60)
