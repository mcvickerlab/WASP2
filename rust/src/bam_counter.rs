use pyo3::prelude::*;
use pyo3::types::PyList;
use rust_htslib::{bam, bam::Read as BamRead, bam::ext::BamRecordExtensions};
use std::collections::HashSet;
use std::path::Path;

/// BAM allele counter using rust-htslib (same backend as pysam)
#[pyclass]
pub struct BamCounter {
    bam_path: String,
}

#[derive(Debug)]
struct Region {
    chrom: String,
    pos: u32,  // 1-based position from Python
    ref_base: char,
    alt_base: char,
}

#[pymethods]
impl BamCounter {
    #[new]
    fn new(bam_path: String) -> PyResult<Self> {
        // Verify BAM file exists
        if !Path::new(&bam_path).exists() {
            return Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
                format!("BAM file not found: {}", bam_path)
            ));
        }

        Ok(BamCounter { bam_path })
    }

    /// Count alleles at SNP positions
    ///
    /// Args:
    ///     regions: List of (chrom, pos, ref, alt) tuples
    ///     min_qual: Minimum base quality (default: 20)
    ///
    /// Returns:
    ///     List of (ref_count, alt_count, other_count) tuples
    #[pyo3(signature = (regions, min_qual=20))]
    fn count_alleles(
        &self,
        py: Python,
        regions: &PyList,
        min_qual: u8,
    ) -> PyResult<Vec<(u32, u32, u32)>> {
        // Parse Python regions
        let mut rust_regions = Vec::new();
        for item in regions.iter() {
            let tuple = item.downcast::<pyo3::types::PyTuple>()?;
            let chrom: String = tuple.get_item(0)?.extract()?;
            let pos: u32 = tuple.get_item(1)?.extract()?;
            let ref_base: String = tuple.get_item(2)?.extract()?;
            let alt_base: String = tuple.get_item(3)?.extract()?;

            rust_regions.push(Region {
                chrom,
                pos,
                ref_base: ref_base.chars().next().unwrap(),
                alt_base: alt_base.chars().next().unwrap(),
            });
        }

        // Release GIL for parallel processing
        py.allow_threads(|| {
            self.count_alleles_impl(&rust_regions, min_qual)
        })
    }
}

impl BamCounter {
    fn count_alleles_impl(
        &self,
        regions: &[Region],
        min_qual: u8,
    ) -> PyResult<Vec<(u32, u32, u32)>> {
        // Open indexed BAM reader
        let mut bam = bam::IndexedReader::from_path(&self.bam_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to open BAM: {}", e)
            ))?;

        // Read header (needed for fetch to work properly)
        let _header = bam::Header::from_template(bam.header());

        // Initialize results
        let mut results = vec![(0u32, 0u32, 0u32); regions.len()];

        // Track which reads we've seen (deduplication across ALL regions)
        // Important: Python keeps one set for all SNPs, so same read only counted once
        let mut seen_reads: HashSet<Vec<u8>> = HashSet::new();

        // Process each region
        for (idx, region) in regions.iter().enumerate() {
            // DO NOT clear seen_reads - must persist across all SNPs!

            // Fetch reads overlapping this position (0-based for pysam compatibility)
            // Python uses 1-based positions, BAM uses 0-based
            let start = region.pos.saturating_sub(1) as i64;
            let end = region.pos as i64;

            // Fetch using region name (rust-htslib accepts string region)
            if let Err(_) = bam.fetch((&region.chrom as &str, start, end)) {
                // Skip if fetch fails (e.g., chromosome not in BAM)
                continue;
            }

            // Iterate through overlapping reads
            for result in bam.records() {
                let record = match result {
                    Ok(r) => r,
                    Err(_) => continue,
                };

                // Skip unmapped reads
                if record.is_unmapped() {
                    continue;
                }

                // Deduplication: skip if we've seen this read
                let qname = record.qname().to_vec();
                if seen_reads.contains(&qname) {
                    continue;
                }
                seen_reads.insert(qname);

                // Find the base at our SNP position using aligned_pairs
                if let Some(base) = get_base_at_position(&record, region.pos - 1, min_qual) {
                    // Count the allele
                    if base == region.ref_base {
                        results[idx].0 += 1;
                    } else if base == region.alt_base {
                        results[idx].1 += 1;
                    } else {
                        results[idx].2 += 1;
                    }
                }
            }
        }

        Ok(results)
    }
}

/// Get base at genomic position, accounting for CIGAR operations
fn get_base_at_position(
    record: &bam::Record,
    target_pos: u32,  // 0-based genomic position
    min_qual: u8,
) -> Option<char> {
    // Get read sequence and qualities
    let seq = record.seq();
    let qual = record.qual();

    // Use aligned_pairs to get CIGAR-aware position mapping
    let aligned_pairs = record.aligned_pairs();

    // Find the query position that aligns to our target reference position
    for pair in aligned_pairs {
        let qpos = pair[0];
        let rpos = pair[1];

        // Check if this is a valid match (not a deletion/insertion)
        if qpos >= 0 && rpos >= 0 && rpos == target_pos as i64 {
            // Found the alignment - check quality
            if qual[qpos as usize] < min_qual {
                return None;
            }

            // Get the base (using array indexing)
            let base = match seq[qpos as usize] {
                b'A' => 'A',
                b'C' => 'C',
                b'G' => 'G',
                b'T' => 'T',
                b'N' => 'N',
                _ => return None,
            };
            return Some(base);
        }
    }

    None
}
