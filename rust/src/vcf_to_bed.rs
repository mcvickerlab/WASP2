//! VCF to BED conversion using noodles
//!
//! Replaces bcftools subprocess with pure Rust implementation for VCF files.
//! BCF files fall back to bcftools due to noodles API complexity.
//!
//! # Performance
//! Expected 5-6x speedup over bcftools subprocess due to:
//! - No process spawn overhead
//! - No pipe overhead
//! - Streaming output with large buffers
//!
//! # Output Format (matches bcftools query)
//! ```text
//! chrom  start  end  ref  alt  genotype
//! chr1   12345  12346  A   G    A|G
//! ```

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// ============================================================================
// Configuration
// ============================================================================

/// Configuration for VCF → BED conversion
#[derive(Debug, Clone)]
pub struct VcfToBedConfig {
    /// Sample names to extract (None = all samples)
    pub samples: Option<Vec<String>>,
    /// Only output heterozygous sites
    pub het_only: bool,
    /// Include indels (not just SNPs)
    pub include_indels: bool,
    /// Maximum indel length (abs(len(ref) - len(alt)))
    pub max_indel_len: usize,
    /// Include genotype column in output
    pub include_genotypes: bool,
}

impl Default for VcfToBedConfig {
    fn default() -> Self {
        Self {
            samples: None,
            het_only: true,
            include_indels: false,
            max_indel_len: 10,
            include_genotypes: true,
        }
    }
}

// ============================================================================
// Genotype Classification
// ============================================================================

/// Genotype classification (matches Python Genotype enum)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    HomRef,  // 0/0, 0|0
    Het,     // 0/1, 1/0, 0|1, 1|0
    HomAlt,  // 1/1, 1|1
    Missing, // ./., .|.
}

// ============================================================================
// Main Entry Point
// ============================================================================

/// Convert VCF to BED format
///
/// Auto-detects VCF vs BCF from file extension.
/// Supports plain VCF and gzipped VCF (.vcf.gz) - BCF returns error.
///
/// # Arguments
/// * `vcf_path` - Input VCF file
/// * `bed_path` - Output BED file
/// * `config` - Conversion configuration
///
/// # Returns
/// Number of variants written, or error for unsupported formats
pub fn vcf_to_bed<P: AsRef<Path>>(
    vcf_path: P,
    bed_path: P,
    config: &VcfToBedConfig,
) -> Result<usize> {
    let vcf_path = vcf_path.as_ref();
    let path_str = vcf_path.to_string_lossy().to_lowercase();

    // Determine format from extension
    let is_bcf = path_str.ends_with(".bcf") || path_str.ends_with(".bcf.gz");
    let is_gzipped = path_str.ends_with(".gz") || path_str.ends_with(".bgz");

    eprintln!(
        "  VCF to BED: {} (bcf={}, gzip={})",
        vcf_path.display(),
        is_bcf,
        is_gzipped
    );

    if is_bcf {
        // BCF not supported in Rust - caller should fall back to bcftools
        return Err(anyhow::anyhow!(
            "BCF format not supported in Rust, use bcftools fallback"
        ));
    } else if is_gzipped {
        vcf_to_bed_vcf_gz(vcf_path, bed_path.as_ref(), config)
    } else {
        vcf_to_bed_vcf_plain(vcf_path, bed_path.as_ref(), config)
    }
}

// ============================================================================
// Plain VCF (uncompressed)
// ============================================================================

fn vcf_to_bed_vcf_plain(
    vcf_path: &Path,
    bed_path: &Path,
    config: &VcfToBedConfig,
) -> Result<usize> {
    let file = File::open(vcf_path).context("Failed to open VCF file")?;
    let reader = BufReader::with_capacity(1024 * 1024, file);

    vcf_to_bed_from_reader(reader, bed_path, config)
}

// ============================================================================
// Gzipped VCF (.vcf.gz, .vcf.bgz)
// ============================================================================

fn vcf_to_bed_vcf_gz(vcf_path: &Path, bed_path: &Path, config: &VcfToBedConfig) -> Result<usize> {
    let file = File::open(vcf_path).context("Failed to open VCF.gz file")?;

    // Try BGZF first (standard for indexed VCF)
    let reader = bgzf::Reader::new(file);
    let buf_reader = BufReader::with_capacity(1024 * 1024, reader);

    vcf_to_bed_from_reader(buf_reader, bed_path, config)
}

// ============================================================================
// Generic VCF Reader (works with plain and gzipped)
// ============================================================================

fn vcf_to_bed_from_reader<R: BufRead>(
    reader: R,
    bed_path: &Path,
    config: &VcfToBedConfig,
) -> Result<usize> {
    let mut vcf_reader = vcf::io::Reader::new(reader);

    let header = vcf_reader
        .read_header()
        .context("Failed to read VCF header")?;

    // Get sample indices
    // When include_genotypes=False, we only need to check one sample to determine
    // if a variant passes filters (het_only, etc.) - no need to output duplicates
    let all_sample_indices = get_sample_indices_from_header(&header, &config.samples)?;
    let sample_indices = if !config.include_genotypes && all_sample_indices.len() > 1 {
        // Only use first sample when not outputting genotypes to avoid duplicates
        vec![all_sample_indices[0]]
    } else {
        all_sample_indices
    };

    eprintln!(
        "  Processing {} samples: {:?}",
        sample_indices.len(),
        config.samples.as_ref().unwrap_or(&vec!["all".to_string()])
    );

    let out_file = File::create(bed_path).context("Failed to create output BED file")?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);

    let mut variant_count = 0;
    let mut total_records = 0;

    for result in vcf_reader.records() {
        let record = result.context("Failed to read VCF record")?;
        total_records += 1;

        if let Some(count) =
            process_vcf_record(&record, &header, &sample_indices, config, &mut writer)?
        {
            variant_count += count;
        }
    }

    writer.flush()?;
    eprintln!(
        "  Processed {} records, wrote {} variants to BED",
        total_records, variant_count
    );

    Ok(variant_count)
}

// ============================================================================
// Record Processing (VCF)
// ============================================================================

fn process_vcf_record<W: Write>(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_indices: &[usize],
    config: &VcfToBedConfig,
    writer: &mut W,
) -> Result<Option<usize>> {
    use vcf::variant::record::AlternateBases;

    // Get reference bases - vcf::Record returns &str directly
    let ref_allele = record.reference_bases().to_string();

    // Get alternate bases
    let alt_bases = record.alternate_bases();

    // Collect all ALT alleles
    let alt_alleles: Vec<String> = alt_bases
        .iter()
        .filter_map(|r| r.ok().map(|a| a.to_string()))
        .collect();

    if alt_alleles.is_empty() {
        return Ok(None); // No valid ALT alleles
    }

    // Get chromosome and position
    let chrom = record.reference_sequence_name();
    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get(), // 1-based
        _ => return Ok(None),
    };
    let pos0 = pos - 1; // 0-based for BED

    // Calculate end position (BED end is exclusive)
    let end = pos0 + ref_allele.len();

    // Special case: when not filtering by het and not including genotypes,
    // output all variants regardless of sample genotypes (like bcftools --drop-genotypes)
    if !config.het_only && !config.include_genotypes {
        let mut written = 0;
        for alt_allele in alt_alleles.iter() {
            // Check SNP vs indel
            let is_snp = ref_allele.len() == 1 && alt_allele.len() == 1;
            if !is_snp && !config.include_indels {
                continue;
            }

            // Check indel length
            if !is_snp {
                let len_diff = (ref_allele.len() as i32 - alt_allele.len() as i32).abs() as usize;
                if len_diff > config.max_indel_len {
                    continue;
                }
            }

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                chrom, pos0, end, ref_allele, alt_allele
            )?;
            written += 1;
        }
        return Ok(Some(written));
    }

    // Process each sample for het filtering or genotype output
    let samples = record.samples();
    let mut written = 0;

    for &sample_idx in sample_indices {
        // Get genotype indices for this sample
        let (gt_indices, is_phased) = get_genotype_indices(&samples, header, sample_idx)?;

        if gt_indices.is_empty() || gt_indices.iter().any(|&i| i.is_none()) {
            continue; // Skip missing genotypes
        }

        let gt_indices: Vec<usize> = gt_indices.iter().filter_map(|&i| i).collect();

        // For multi-allelic sites, we output each heterozygous ALT allele separately
        // This matches bcftools -g het behavior
        for (alt_idx, alt_allele) in alt_alleles.iter().enumerate() {
            let alt_index = alt_idx + 1; // ALT indices are 1-based (0 = REF)

            // Check if this sample is heterozygous for this specific ALT
            // Het means one allele is REF (0) and one is this ALT
            let has_ref = gt_indices.iter().any(|&i| i == 0);
            let has_this_alt = gt_indices.iter().any(|&i| i == alt_index);
            let is_het_for_this_alt = has_ref && has_this_alt;

            // Also handle het between two different ALTs (e.g., 1/2)
            // In this case, we should still output each ALT allele
            let num_different_alleles = gt_indices
                .iter()
                .collect::<std::collections::HashSet<_>>()
                .len();
            let is_het_multi_alt = num_different_alleles > 1 && has_this_alt;

            let is_het = is_het_for_this_alt || is_het_multi_alt;

            // Filter het-only
            if config.het_only && !is_het {
                continue;
            }

            // Check SNP vs indel for this specific ALT
            let is_snp = ref_allele.len() == 1 && alt_allele.len() == 1;
            if !is_snp && !config.include_indels {
                continue; // Skip indels if not requested
            }

            // Check indel length
            if !is_snp {
                let len_diff = (ref_allele.len() as i32 - alt_allele.len() as i32).abs() as usize;
                if len_diff > config.max_indel_len {
                    continue;
                }
            }

            // Build genotype string (e.g., "A|G")
            let gt_string =
                build_genotype_string(&ref_allele, &alt_alleles, &gt_indices, is_phased);

            // Write BED line
            if config.include_genotypes {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    chrom, pos0, end, ref_allele, alt_allele, gt_string
                )?;
            } else {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    chrom, pos0, end, ref_allele, alt_allele
                )?;
            }

            written += 1;
        }
    }

    Ok(Some(written))
}

/// Get genotype indices from sample (returns allele indices like [0, 1] for 0/1)
fn get_genotype_indices(
    samples: &vcf::record::Samples,
    header: &vcf::Header,
    sample_idx: usize,
) -> Result<(Vec<Option<usize>>, bool)> {
    use vcf::variant::record::samples::keys::key::GENOTYPE as GT_KEY;
    use vcf::variant::record::samples::Sample as SampleTrait;

    // Get sample at index
    let sample = match samples.iter().nth(sample_idx) {
        Some(s) => s,
        None => return Ok((vec![], false)),
    };

    // Try to get GT field from sample
    let gt_value = match sample.get(header, GT_KEY) {
        Some(Ok(Some(v))) => v,
        _ => return Ok((vec![], false)),
    };

    // Convert value to string using Debug and parse manually
    let gt_string = format!("{:?}", gt_value);
    let gt_clean = extract_genotype_string(&gt_string);

    // Check for missing genotype
    if gt_clean.contains('.') {
        return Ok((vec![None], false));
    }

    // Parse genotype - format is "0|1", "0/1", etc.
    let is_phased = gt_clean.contains('|');

    let indices: Vec<Option<usize>> = gt_clean
        .split(|c| c == '|' || c == '/')
        .map(|s| s.parse().ok())
        .collect();

    Ok((indices, is_phased))
}

/// Build genotype string from allele indices (e.g., [0, 1] -> "A|G")
fn build_genotype_string(
    ref_allele: &str,
    alt_alleles: &[String],
    gt_indices: &[usize],
    is_phased: bool,
) -> String {
    let allele_strs: Vec<String> = gt_indices
        .iter()
        .map(|&idx| {
            if idx == 0 {
                ref_allele.to_string()
            } else if idx <= alt_alleles.len() {
                alt_alleles[idx - 1].clone()
            } else {
                idx.to_string() // Fallback
            }
        })
        .collect();

    allele_strs.join(if is_phased { "|" } else { "/" })
}

// ============================================================================
// Genotype String Extraction
// ============================================================================

/// Extract genotype string from Debug format
/// Handles formats like: Genotype(Genotype("0|1")), String("0|1"), "0|1"
fn extract_genotype_string(debug_str: &str) -> String {
    // Find the innermost quoted string
    if let Some(start) = debug_str.rfind('"') {
        if let Some(end) = debug_str[..start].rfind('"') {
            return debug_str[end + 1..start].to_string();
        }
    }

    // Fallback: try to find pattern like 0|1 or 0/1
    for part in debug_str.split(|c: char| !c.is_ascii_digit() && c != '|' && c != '/' && c != '.') {
        let trimmed = part.trim();
        if !trimmed.is_empty() && (trimmed.contains('|') || trimmed.contains('/')) {
            return trimmed.to_string();
        }
    }

    // If all else fails, return as-is
    debug_str.to_string()
}

// ============================================================================
// Sample Index Lookup
// ============================================================================

fn get_sample_indices_from_header(
    header: &vcf::Header,
    requested: &Option<Vec<String>>,
) -> Result<Vec<usize>> {
    let sample_names = header.sample_names();

    match requested {
        Some(names) => {
            let mut indices = Vec::with_capacity(names.len());
            for name in names {
                let idx = sample_names.iter().position(|s| s == name).ok_or_else(|| {
                    anyhow::anyhow!(
                        "Sample '{}' not found in VCF. Available: {:?}",
                        name,
                        sample_names.iter().take(5).collect::<Vec<_>>()
                    )
                })?;
                indices.push(idx);
            }
            Ok(indices)
        }
        None => Ok((0..sample_names.len()).collect()),
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::NamedTempFile;

    fn create_test_vcf() -> NamedTempFile {
        let mut vcf = NamedTempFile::new().unwrap();
        writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf, "##contig=<ID=chr1,length=1000000>").unwrap();
        writeln!(
            vcf,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            vcf,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"
        )
        .unwrap();
        writeln!(vcf, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(vcf, "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t1|1").unwrap(); // HomAlt - should be filtered
        writeln!(vcf, "chr1\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(vcf, "chr1\t400\t.\tAT\tA\t.\t.\t.\tGT\t0|1").unwrap(); // Deletion - skipped by default
        vcf.flush().unwrap();
        vcf
    }

    #[test]
    fn test_vcf_to_bed_het_only() {
        let vcf = create_test_vcf();
        let bed = NamedTempFile::new().unwrap();

        let config = VcfToBedConfig {
            samples: Some(vec!["SAMPLE1".to_string()]),
            het_only: true,
            include_indels: false,
            max_indel_len: 10,
            include_genotypes: true,
        };

        let count = vcf_to_bed(vcf.path(), bed.path(), &config).unwrap();

        // Should have 2 het SNPs (pos 100 and 300), skipping homalt at 200 and indel at 400
        assert_eq!(count, 2);

        // Read output
        let content = std::fs::read_to_string(bed.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr1\t99\t100\tA\tG"));
        assert!(lines[1].starts_with("chr1\t299\t300\tG\tA"));
    }

    #[test]
    fn test_vcf_to_bed_with_indels() {
        let vcf = create_test_vcf();
        let bed = NamedTempFile::new().unwrap();

        let config = VcfToBedConfig {
            samples: Some(vec!["SAMPLE1".to_string()]),
            het_only: true,
            include_indels: true,
            max_indel_len: 10,
            include_genotypes: true,
        };

        let count = vcf_to_bed(vcf.path(), bed.path(), &config).unwrap();

        // Should have 3 het variants (2 SNPs + 1 deletion)
        assert_eq!(count, 3);
    }

    #[test]
    fn test_vcf_to_bed_all_genotypes() {
        let vcf = create_test_vcf();
        let bed = NamedTempFile::new().unwrap();

        let config = VcfToBedConfig {
            samples: Some(vec!["SAMPLE1".to_string()]),
            het_only: false, // Include all genotypes
            include_indels: false,
            max_indel_len: 10,
            include_genotypes: true,
        };

        let count = vcf_to_bed(vcf.path(), bed.path(), &config).unwrap();

        // Should have 3 SNPs (het at 100, homalt at 200, het at 300)
        assert_eq!(count, 3);
    }

    /// Test that multi-allelic heterozygous sites are properly included
    /// This was the root cause of the 2,167 missing variants in WASP2-Rust
    #[test]
    fn test_vcf_to_bed_multiallelic() {
        let mut vcf = NamedTempFile::new().unwrap();
        writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf, "##contig=<ID=chr1,length=1000000>").unwrap();
        writeln!(
            vcf,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            vcf,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"
        )
        .unwrap();
        // Biallelic het (baseline)
        writeln!(vcf, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        // Multi-allelic: C -> A,T with het for first ALT (0|1 = het C/A)
        writeln!(vcf, "chr1\t200\t.\tC\tA,T\t.\t.\t.\tGT\t0|1").unwrap();
        // Multi-allelic: G -> A,C with het for second ALT (0|2 = het G/C)
        writeln!(vcf, "chr1\t300\t.\tG\tA,C\t.\t.\t.\tGT\t0|2").unwrap();
        // Multi-allelic: het between two ALTs (1|2 = het A/T)
        writeln!(vcf, "chr1\t400\t.\tT\tA,G\t.\t.\t.\tGT\t1|2").unwrap();
        // Multi-allelic: hom ref (0|0) - should be filtered by het_only
        writeln!(vcf, "chr1\t500\t.\tA\tG,C\t.\t.\t.\tGT\t0|0").unwrap();
        vcf.flush().unwrap();

        let bed = NamedTempFile::new().unwrap();

        let config = VcfToBedConfig {
            samples: Some(vec!["SAMPLE1".to_string()]),
            het_only: true,
            include_indels: false,
            max_indel_len: 10,
            include_genotypes: true,
        };

        let count = vcf_to_bed(vcf.path(), bed.path(), &config).unwrap();

        // Should include:
        // - pos 100: 1 het SNP (biallelic)
        // - pos 200: 1 het for ALT A (0|1)
        // - pos 300: 1 het for ALT C (0|2)
        // - pos 400: 2 hets for ALT A and ALT G (1|2 is het for both)
        // Total: 5 het entries
        assert_eq!(count, 5);

        // Read output and verify
        let content = std::fs::read_to_string(bed.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 5);

        // Verify multi-allelic sites are present
        assert!(
            lines.iter().any(|l| l.contains("chr1\t199\t200\tC\tA")),
            "Missing multi-allelic het 0|1 for A"
        );
        assert!(
            lines.iter().any(|l| l.contains("chr1\t299\t300\tG\tC")),
            "Missing multi-allelic het 0|2 for C"
        );
    }
}
