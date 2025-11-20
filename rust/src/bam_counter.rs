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
    chrom: String,
    start: i64,
    end: i64,
    snps: Vec<(usize, Region)>,  // (original_index, region)
}

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

        // Group regions into genomic windows
        let windows = self.group_into_windows(regions);

        // Track which reads we've seen (deduplication across ALL windows)
        let mut seen_reads: HashSet<Vec<u8>> = HashSet::new();

        // Process each window with ONE fetch per window
        for window in windows {
            // Fetch all reads in this window at once
            if let Err(_) = bam.fetch((&window.chrom as &str, window.start, window.end)) {
                // Skip if fetch fails (e.g., chromosome not in BAM)
                continue;
            }

            // Process all reads in this window
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

                // Check this read against ALL SNPs in the window
                let mut counted_any = false;
                for (idx, region) in &window.snps {
                    // Check if read overlaps this SNP position
                    if let Some(base) = get_base_at_position(&record, region.pos - 1, min_qual) {
                        counted_any = true;
                        // Count the allele
                        if base == region.ref_base {
                            results[*idx].0 += 1;
                        } else if base == region.alt_base {
                            results[*idx].1 += 1;
                        } else {
                            results[*idx].2 += 1;
                        }
                    }
                }

                // Only add to seen_reads if we actually counted this read at least once
                // This prevents marking reads as "seen" if they don't overlap any SNPs in this window
                if counted_any {
                    seen_reads.insert(qname);
                }
            }
        }

        Ok(results)
    }

    /// Group regions into genomic windows to reduce fetch() calls
    fn group_into_windows(&self, regions: &[Region]) -> Vec<GenomicWindow> {
        const WINDOW_SIZE: i64 = 100_000; // 100Kb windows

        let mut windows: Vec<GenomicWindow> = Vec::new();
        let mut indexed_regions: Vec<(usize, Region)> = regions
            .iter()
            .enumerate()
            .map(|(i, r)| (i, r.clone()))
            .collect();

        // Sort by chromosome and position
        indexed_regions.sort_by(|a, b| {
            a.1.chrom.cmp(&b.1.chrom)
                .then(a.1.pos.cmp(&b.1.pos))
        });

        let mut current_window: Option<GenomicWindow> = None;

        for (idx, region) in indexed_regions {
            let pos_i64 = region.pos as i64;

            match &mut current_window {
                Some(window) if window.chrom == region.chrom && pos_i64 < window.end => {
                    // Add to current window
                    window.snps.push((idx, region.clone()));
                    // Extend window if needed
                    if pos_i64 > window.end {
                        window.end = pos_i64 + 1;
                    }
                }
                Some(window) => {
                    // Different chromosome or too far away, finalize current window
                    windows.push(GenomicWindow {
                        chrom: window.chrom.clone(),
                        start: window.start,
                        end: window.end,
                        snps: window.snps.clone(),
                    });

                    // Start new window
                    current_window = Some(GenomicWindow {
                        chrom: region.chrom.clone(),
                        start: (pos_i64 - 1).max(0),
                        end: pos_i64 + WINDOW_SIZE,
                        snps: vec![(idx, region)],
                    });
                }
                None => {
                    // First window
                    current_window = Some(GenomicWindow {
                        chrom: region.chrom.clone(),
                        start: (pos_i64 - 1).max(0),
                        end: pos_i64 + WINDOW_SIZE,
                        snps: vec![(idx, region)],
                    });
                }
            }
        }

        // Add last window
        if let Some(window) = current_window {
            windows.push(window);
        }

        windows
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
