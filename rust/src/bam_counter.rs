use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use std::collections::HashMap;
use std::path::Path;

/// BAM allele counter using noodles
#[pyclass]
pub struct BamCounter {
    bam_path: String,
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

        // Parse regions
        let rust_regions: Vec<Region> = regions
            .iter()
            .map(|item| {
                let tuple = item.downcast::<PyTuple>()?;

                let chrom: String = tuple.get_item(0)?.extract()?;
                let pos: u32 = tuple.get_item(1)?.extract()?;
                let ref_str: String = tuple.get_item(2)?.extract()?;
                let alt_str: String = tuple.get_item(3)?.extract()?;

                Ok(Region {
                    chrom,
                    pos,
                    ref_base: ref_str.chars().next().unwrap_or('N'),
                    alt_base: alt_str.chars().next().unwrap_or('N'),
                })
            })
            .collect::<PyResult<Vec<_>>>()?;

        // Release GIL for parallel processing
        py.allow_threads(|| {
            self.count_alleles_impl(&rust_regions, min_qual)
        })
    }

    /// Get BAM file path
    fn get_bam_path(&self) -> String {
        self.bam_path.clone()
    }
}

#[derive(Clone, Debug)]
struct Region {
    chrom: String,
    pos: u32,
    ref_base: char,
    alt_base: char,
}

impl BamCounter {
    fn count_alleles_impl(
        &self,
        regions: &[Region],
        min_qual: u8,
    ) -> PyResult<Vec<(u32, u32, u32)>> {
        use noodles::bam;
        use std::fs::File;

        // Open BAM file
        let file = File::open(&self.bam_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to open BAM: {}", e)
            ))?;

        let mut reader = bam::io::Reader::new(file);

        // Read header
        let _header = reader.read_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to read header: {}", e)
            ))?;

        // Create a hashmap of position -> (region_idx, region)
        let mut pos_map: HashMap<(String, u32), Vec<(usize, &Region)>> = HashMap::new();
        for (idx, region) in regions.iter().enumerate() {
            pos_map
                .entry((region.chrom.clone(), region.pos))
                .or_insert_with(Vec::new)
                .push((idx, region));
        }

        // Initialize results
        let mut results = vec![(0u32, 0u32, 0u32); regions.len()];

        // Iterate through all reads
        for result in reader.records() {
            let record = result.map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Failed to read record: {}", e)
                )
            })?;

            // Get alignment info
            let Some(_ref_id) = record.reference_sequence_id() else {
                continue;
            };
            let Some(Ok(align_start_pos)) = record.alignment_start() else {
                continue;
            };
            let align_start = align_start_pos.get() as u32;
            let align_end = align_start + record.sequence().len() as u32;

            // Check if any of our regions overlap this read
            for ((_chrom, pos), region_list) in &pos_map {
                // Note: chromosome filtering would require header lookup to map ref_id to name
                // For now, we'll check position overlap
                if *pos >= align_start && *pos < align_end {
                    // This read overlaps our position
                    for (idx, region) in region_list {
                        if let Some(base) = get_base_at_pos(&record, *pos, min_qual) {
                            if base == region.ref_base {
                                results[*idx].0 += 1;
                            } else if base == region.alt_base {
                                results[*idx].1 += 1;
                            } else {
                                results[*idx].2 += 1;
                            }
                        }
                    }
                }
            }
        }

        Ok(results)
    }
}

/// Extract base from read at specific genomic position
fn get_base_at_pos(
    record: &noodles::bam::Record,
    target_pos: u32,
    min_qual: u8,
) -> Option<char> {
    use noodles::sam::alignment::record::QualityScores;

    // Get alignment start position
    let align_start = record.alignment_start()?.ok()?.get() as u32;

    // Calculate offset in read (0-based)
    if target_pos < align_start {
        return None;
    }

    let offset = (target_pos - align_start) as usize;

    // Get sequence and quality scores
    let seq = record.sequence();
    let qual = record.quality_scores();

    // Check bounds
    if offset >= seq.len() {
        return None;
    }

    // Check quality - iterate to find the quality score at offset
    let qual_vec: Vec<u8> = qual.iter().collect();
    if offset < qual_vec.len() {
        if qual_vec[offset] < min_qual {
            return None;
        }
    } else {
        return None;
    }

    // Get base - iterate to find the base at offset
    let bases: Vec<u8> = seq.iter().collect();
    if offset < bases.len() {
        Some(bases[offset] as char)
    } else {
        None
    }
}
