use pyo3::prelude::*;
use pyo3::types::PyList;
use rust_htslib::{bam, bam::Read as BamRead, bam::ext::BamRecordExtensions};
use std::collections::HashSet;
use std::path::Path;

/// BAM allele counter using rust-htslib with batched fetching
#[pyclass]
pub struct BamCounter {
    bam_path: String,
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    pos: u32,  // 1-based position from Python
    ref_base: char,
    alt_base: char,
}

/// Genomic window containing multiple SNPs
struct GenomicWindow {

// PyO3 expands #[pymethods] into impl blocks that trigger non_local_definitions warnings;
// suppress the noise until we restructure.
#[allow(non_local_definitions)]
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

    /// Count alleles at SNP positions using batched fetching
    ///
    /// Args:
    ///     regions: List of (chrom, pos, ref, alt) tuples
    ///     min_qual: Minimum base quality (default: 0 for WASP2 compatibility)
    ///
    /// Returns:
    ///     List of (ref_count, alt_count, other_count) tuples
    #[pyo3(signature = (regions, min_qual=0))]
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

        // Initialize results
        let mut results = vec![(0u32, 0u32, 0u32); regions.len()];

        // Process per chromosome with per-chrom dedup (matches Python baseline)
        let mut chrom_regions: Vec<&Region> = regions.iter().collect();
        chrom_regions.sort_by(|a, b| a.chrom.cmp(&b.chrom));

        let mut current_chrom: Option<String> = None;
        let mut seen: HashSet<Vec<u8>> = HashSet::new();

        for (idx, region) in regions.iter().enumerate() {
            // Reset seen set when chromosome changes
            if current_chrom.as_ref() != Some(&region.chrom) {
                current_chrom = Some(region.chrom.clone());
                seen.clear();
            }

            // Fetch reads overlapping this SNP (htslib uses 0-based, half-open)
            let start = (region.pos - 1) as i64;
            let end = region.pos as i64; // exclusive
            if bam.fetch((&region.chrom as &str, start, end)).is_err() {
                continue;
            }

            for result in bam.records() {
                let record = match result {
                    Ok(r) => r,
                    Err(_) => continue,
                };

                if record.is_unmapped() {
                    continue;
                }

                let qname = record.qname().to_vec();
                if seen.contains(&qname) {
                    continue;
                }

                if let Some(base) = get_base_at_position(&record, region.pos - 1, min_qual) {
                    seen.insert(qname);
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
/// Matches WASP2 behavior: NO quality filtering by default
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
            // Optional quality filtering (min_qual=0 means no filtering like WASP2)
            if min_qual > 0 && qual[qpos as usize] < min_qual {
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

#[cfg(test)]
mod tests {
    use super::{BamCounter, Region};

    #[test]
    fn groups_regions_into_chrom_sorted_windows() {
        let counter = BamCounter { bam_path: "dummy.bam".to_string() };
        let regions = vec![
            Region { chrom: "chr1".into(), pos: 10, ref_base: 'A', alt_base: 'G' },
            Region { chrom: "chr1".into(), pos: 20, ref_base: 'C', alt_base: 'T' },
            // This one is >100kb away, should start a new window
            Region { chrom: "chr1".into(), pos: 150_500, ref_base: 'G', alt_base: 'A' },
            // Different chromosome, always new window
            Region { chrom: "chr2".into(), pos: 5, ref_base: 'T', alt_base: 'C' },
        ];

        let windows = counter.group_into_windows(&regions);
        assert_eq!(windows.len(), 3, "expected three windows (two on chr1, one on chr2)");

        // Windows should be sorted by chrom then pos
        assert_eq!(windows[0].chrom, "chr1");
        assert_eq!(windows[1].chrom, "chr1");
        assert_eq!(windows[2].chrom, "chr2");

        // First window has the two close chr1 SNPs
        assert_eq!(windows[0].snps.len(), 2);
        // Second window has the distant chr1 SNP
        assert_eq!(windows[1].snps.len(), 1);
        // Third window has the chr2 SNP
        assert_eq!(windows[2].snps.len(), 1);
    }
}
