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
        return Err(anyhow::anyhow!("BCF format not supported in Rust, use bcftools fallback"));
    } else if is_gzipped {
        vcf_to_bed_vcf_gz(vcf_path, bed_path.as_ref(), config)
    } else {
        vcf_to_bed_vcf_plain(vcf_path, bed_path.as_ref(), config)
    }
}

// ============================================================================
// Plain VCF (uncompressed)
// ============================================================================

fn vcf_to_bed_vcf_plain(vcf_path: &Path, bed_path: &Path, config: &VcfToBedConfig) -> Result<usize> {
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

    let header = vcf_reader.read_header().context("Failed to read VCF header")?;

    // Get sample indices
    let sample_indices = get_sample_indices_from_header(&header, &config.samples)?;

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

        if let Some(count) = process_vcf_record(&record, &header, &sample_indices, config, &mut writer)? {
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

    // Check biallelic (exactly one ALT allele)
    let alt_count = alt_bases.iter().count();
    if alt_count != 1 {
        return Ok(None); // Skip multi-allelic
    }

    // Get the single ALT allele
    let alt_result = alt_bases.iter().next().unwrap();
    let alt_allele = match alt_result {
        Ok(alt) => alt.to_string(),
        Err(_) => return Ok(None), // Skip if can't parse ALT
    };

    // Check SNP vs indel
    let is_snp = ref_allele.len() == 1 && alt_allele.len() == 1;
    if !is_snp && !config.include_indels {
        return Ok(None); // Skip indels if not requested
    }

    // Check indel length
    if !is_snp {
        let len_diff = (ref_allele.len() as i32 - alt_allele.len() as i32).abs() as usize;
        if len_diff > config.max_indel_len {
            return Ok(None);
        }
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

    // Process each sample
    let samples = record.samples();
    let mut written = 0;

    for &sample_idx in sample_indices {
        // Get genotype for this sample
        let (genotype, gt_string) = get_genotype_from_samples(&samples, header, sample_idx, &ref_allele, &alt_allele)?;

        // Filter het-only
        if config.het_only && genotype != Genotype::Het {
            continue;
        }

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

    Ok(Some(written))
}

// ============================================================================
// Genotype Extraction - Simple String-Based Approach
// ============================================================================

fn get_genotype_from_samples(
    samples: &vcf::record::Samples,
    header: &vcf::Header,
    sample_idx: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<(Genotype, String)> {
    use vcf::variant::record::samples::keys::key::GENOTYPE as GT_KEY;
    use vcf::variant::record::samples::Sample as SampleTrait;

    // Get sample at index
    let sample = match samples.iter().nth(sample_idx) {
        Some(s) => s,
        None => return Ok((Genotype::Missing, "./.".to_string())),
    };

    // Try to get GT field from sample
    let gt_value = match sample.get(header, GT_KEY) {
        Some(Ok(Some(v))) => v,
        _ => return Ok((Genotype::Missing, "./.".to_string())),
    };

    // Convert value to string using Debug and parse manually
    // Debug format is like: Genotype(Genotype("0|1"))
    let gt_string = format!("{:?}", gt_value);

    // Extract the actual genotype string from debug format
    // Format: Genotype(Genotype("0|1")) or String("0|1") or just "0|1"
    let gt_clean = extract_genotype_string(&gt_string);

    parse_genotype_string(&gt_clean, ref_allele, alt_allele)
}

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

/// Parse a genotype string (like "0|1") into classification and allele representation
fn parse_genotype_string(
    gt_string: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<(Genotype, String)> {
    // Check for missing genotype
    if gt_string.contains('.') {
        return Ok((Genotype::Missing, "./.".to_string()));
    }

    // Parse genotype - format is "0|1", "0/1", etc.
    // Determine if phased (|) or unphased (/)
    let is_phased = gt_string.contains('|');

    let allele_indices: Vec<usize> = gt_string
        .split(|c| c == '|' || c == '/')
        .filter_map(|s| s.parse().ok())
        .collect();

    if allele_indices.is_empty() {
        return Ok((Genotype::Missing, "./.".to_string()));
    }

    // Count alt alleles
    let alt_count = allele_indices.iter().filter(|&&idx| idx > 0).count();

    let genotype = match alt_count {
        0 => Genotype::HomRef,
        n if n == allele_indices.len() => Genotype::HomAlt,
        _ => Genotype::Het,
    };

    // Build allele-based string representation (like A|G)
    let alleles = [ref_allele, alt_allele];
    let allele_strs: Vec<String> = allele_indices
        .iter()
        .map(|&idx| {
            if idx < alleles.len() {
                alleles[idx].to_string()
            } else {
                idx.to_string() // Fallback for multi-allelic
            }
        })
        .collect();

    let gt_output = allele_strs.join(if is_phased { "|" } else { "/" });
    Ok((genotype, gt_output))
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
                let idx = sample_names
                    .iter()
                    .position(|s| s == name)
                    .ok_or_else(|| {
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
        writeln!(vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
        writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1").unwrap();
        writeln!(vcf, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(vcf, "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t1|1").unwrap();  // HomAlt - should be filtered
        writeln!(vcf, "chr1\t300\t.\tG\tA\t.\t.\t.\tGT\t0|1").unwrap();
        writeln!(vcf, "chr1\t400\t.\tAT\tA\t.\t.\t.\tGT\t0|1").unwrap();  // Deletion - skipped by default
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
            het_only: false,  // Include all genotypes
            include_indels: false,
            max_indel_len: 10,
            include_genotypes: true,
        };

        let count = vcf_to_bed(vcf.path(), bed.path(), &config).unwrap();

        // Should have 3 SNPs (het at 100, homalt at 200, het at 300)
        assert_eq!(count, 3);
    }
}
