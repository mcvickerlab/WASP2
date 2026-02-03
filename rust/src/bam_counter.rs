use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rayon::prelude::*;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Read as BamRead};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;

/// BAM allele counter using rust-htslib with batched fetching
#[pyclass]
pub struct BamCounter {
    bam_path: String,
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    pos: u32, // 1-based position from Python
    ref_allele: String,  // Full reference allele (supports INDELs)
    alt_allele: String,  // Full alternate allele (supports INDELs)
}

impl Region {
    /// Returns true if this variant is a simple SNP (single base change)
    fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() == 1
    }
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
                format!("BAM file not found: {}", bam_path),
            ));
        }

        Ok(BamCounter { bam_path })
    }

    /// Count alleles at variant positions (SNPs and INDELs) using batched fetching
    ///
    /// Args:
    ///     regions: List of (chrom, pos, ref, alt) tuples - supports both SNPs and INDELs
    ///     min_qual: Minimum base quality (default: 0 for WASP2 compatibility)
    ///     threads: Number of worker threads (default: 1). Use >1 to enable Rayon parallelism per chromosome.
    ///
    /// Returns:
    ///     List of (ref_count, alt_count, other_count) tuples
    #[pyo3(signature = (regions, min_qual=0, threads=1))]
    fn count_alleles(
        &self,
        py: Python<'_>,
        regions: &Bound<'_, PyList>,
        min_qual: u8,
        threads: usize,
    ) -> PyResult<Vec<(u32, u32, u32)>> {
        // Parse Python regions (supports both SNPs and INDELs)
        let mut rust_regions = Vec::new();
        for item in regions.iter() {
            let tuple = item.cast::<PyTuple>()?;
            let chrom: String = tuple.get_item(0)?.extract()?;
            let pos: u32 = tuple.get_item(1)?.extract()?;
            let ref_allele: String = tuple.get_item(2)?.extract()?;
            let alt_allele: String = tuple.get_item(3)?.extract()?;

            // Validate alleles are non-empty
            if ref_allele.is_empty() {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("Empty ref_allele for variant at {}:{}", chrom, pos)
                ));
            }
            if alt_allele.is_empty() {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("Empty alt_allele for variant at {}:{}", chrom, pos)
                ));
            }

            rust_regions.push(Region {
                chrom,
                pos,
                ref_allele,
                alt_allele,
            });
        }

        // Release GIL for parallel processing
        py.detach(|| self.count_alleles_impl(&rust_regions, min_qual, threads))
    }
}

impl BamCounter {
    fn count_alleles_impl(
        &self,
        regions: &[Region],
        min_qual: u8,
        threads: usize,
    ) -> PyResult<Vec<(u32, u32, u32)>> {
        // Initialize results
        let mut results = vec![(0u32, 0u32, 0u32); regions.len()];

        // Group regions by chromosome while preserving encounter order
        let grouped = self.group_regions_by_chrom(regions);
        let debug_sites = parse_debug_sites();

        // Process chromosomes in parallel if threads > 1
        if threads > 1 {
            // Set thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Failed to create thread pool: {}",
                        e
                    ))
                })?
                .install(|| {
                    // Process chromosomes in parallel
                    let partial_results: Result<Vec<_>, _> = grouped
                        .par_iter()
                        .map(|(chrom, chrom_regions)| {
                            self.process_chromosome_reads(
                                chrom,
                                chrom_regions,
                                min_qual,
                                &debug_sites,
                            )
                        })
                        .collect();

                    // Merge results
                    for partial in partial_results? {
                        for (idx, (r, a, o)) in partial {
                            let entry = &mut results[idx];
                            entry.0 += r;
                            entry.1 += a;
                            entry.2 += o;
                        }
                    }
                    Ok::<(), PyErr>(())
                })?;
        } else {
            // Single-threaded path
            for (chrom, chrom_regions) in grouped {
                let partial =
                    self.process_chromosome_reads(&chrom, &chrom_regions, min_qual, &debug_sites)?;
                for (idx, (r, a, o)) in partial {
                    let entry = &mut results[idx];
                    entry.0 += r;
                    entry.1 += a;
                    entry.2 += o;
                }
            }
        }

        Ok(results)
    }

    /// Process a single chromosome by scanning reads once, honoring encounter order per read
    fn process_chromosome_reads(
        &self,
        chrom: &str,
        regions: &[(usize, Region)],
        min_qual: u8,
        debug_sites: &FxHashMap<(String, u32), usize>,
    ) -> PyResult<FxHashMap<usize, (u32, u32, u32)>> {
        let mut bam = bam::IndexedReader::from_path(&self.bam_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM: {}", e))
        })?;

        let mut seen_reads: FxHashSet<Vec<u8>> = FxHashSet::default();
        let total_variants: usize = regions.len();
        let mut counts: FxHashMap<usize, (u32, u32, u32)> = FxHashMap::default();
        counts.reserve(total_variants);

        // Track skipped records for logging (prevents silent failures)
        let mut skipped_records: u32 = 0;
        const MAX_SKIP_WARNINGS: u32 = 5;

        // Build position -> variant list, preserving encounter order
        let mut pos_map: FxHashMap<u32, Vec<(usize, Region)>> = FxHashMap::default();
        let mut min_pos: u32 = u32::MAX;
        let mut max_pos: u32 = 0;
        for (idx, region) in regions.iter() {
            pos_map
                .entry(region.pos)
                .or_insert_with(Vec::new)
                .push((*idx, region.clone()));
            if region.pos < min_pos {
                min_pos = region.pos;
            }
            if region.pos > max_pos {
                max_pos = region.pos;
            }
        }

        if pos_map.is_empty() {
            return Ok(counts);
        }

        // Fetch the span covering all SNPs on this chromosome
        let start = if min_pos == 0 {
            0
        } else {
            (min_pos - 1) as i64
        };
        let end = max_pos.saturating_add(1) as i64;
        if let Err(e) = bam.fetch((chrom, start, end)) {
            eprintln!(
                "[WARN] Failed to fetch {}:{}-{}: {}. Skipping {} variants.",
                chrom, start, end, e, regions.len()
            );
            return Ok(counts);
        }

        // For each read, assign to the earliest SNP in encounter order that it overlaps
        let mut read_iter = bam.records();
        while let Some(res) = read_iter.next() {
            let record = match res {
                Ok(r) => r,
                Err(e) => {
                    skipped_records += 1;
                    if skipped_records <= MAX_SKIP_WARNINGS {
                        eprintln!(
                            "[WARN] Skipped corrupted BAM record on {}: {}",
                            chrom, e
                        );
                    }
                    continue;
                }
            };
            if record.is_unmapped()
                || record.is_secondary()
                || record.is_supplementary()
                || record.is_duplicate()
            {
                continue;
            }
            let qname = record.qname().to_vec();
            if seen_reads.contains(&qname) {
                continue;
            }

            // Find earliest-overlap SNP by encounter index
            let mut best: Option<(usize, &Region, usize, u32)> = None; // (encounter_idx, region, qpos, pos1)
            for pair in record.aligned_pairs() {
                let qpos = pair[0];
                let rpos = pair[1];
                if qpos < 0 || rpos < 0 {
                    continue;
                }
                let pos1 = (rpos as u32).saturating_add(1);
                if let Some(list) = pos_map.get(&pos1) {
                    for (enc_idx, region) in list {
                        if let Some((best_idx, _, _, _)) = best {
                            if *enc_idx >= best_idx {
                                continue;
                            }
                        }
                        best = Some((*enc_idx, region, qpos as usize, pos1));
                    }
                }
            }

            if let Some((enc_idx, region, qpos, pos1)) = best {
                let quals = record.qual();
                if min_qual > 0 {
                    if qpos >= quals.len() || quals[qpos] < min_qual {
                        continue;
                    }
                }

                let entry_counts = counts.entry(enc_idx).or_insert((0, 0, 0));
                let read_allele: String;

                if region.is_snp() {
                    // Fast path for SNPs: single base comparison
                    read_allele = match record.seq()[qpos] {
                        b'A' => "A".to_string(),
                        b'C' => "C".to_string(),
                        b'G' => "G".to_string(),
                        b'T' => "T".to_string(),
                        b'N' => "N".to_string(),
                        _ => continue,
                    };
                } else {
                    // INDEL path: extract sequence span from read
                    // For INDELs, we need to determine which allele the read supports
                    // by comparing the read sequence to expected ref/alt alleles
                    let seq = record.seq();
                    let seq_len = seq.len();

                    // Extract enough bases to compare against both alleles
                    let max_allele_len = std::cmp::max(region.ref_allele.len(), region.alt_allele.len());
                    let end_pos = std::cmp::min(qpos + max_allele_len, seq_len);

                    if qpos >= seq_len {
                        continue;
                    }

                    // Build read sequence string from the relevant positions
                    let mut seq_string = String::with_capacity(end_pos - qpos);
                    for i in qpos..end_pos {
                        let base = match seq[i] {
                            b'A' => 'A',
                            b'C' => 'C',
                            b'G' => 'G',
                            b'T' => 'T',
                            _ => 'N',
                        };
                        seq_string.push(base);
                    }
                    read_allele = seq_string;
                }

                // Compare read allele to ref/alt
                if read_allele == region.ref_allele {
                    entry_counts.0 += 1;
                } else if read_allele == region.alt_allele {
                    entry_counts.1 += 1;
                } else if region.is_snp() {
                    // For SNPs, a mismatch to both is "other"
                    entry_counts.2 += 1;
                } else {
                    // For INDELs, check if read starts with either allele
                    // This handles cases where read extends beyond the variant
                    if read_allele.starts_with(&region.ref_allele) {
                        entry_counts.0 += 1;
                    } else if read_allele.starts_with(&region.alt_allele) {
                        entry_counts.1 += 1;
                    } else {
                        entry_counts.2 += 1;
                    }
                }
                seen_reads.insert(qname.clone());

                if let Some(limit) = debug_sites.get(&(chrom.to_string(), pos1)) {
                    if *limit > 0
                        && entry_counts.0 + entry_counts.1 + entry_counts.2 <= *limit as u32
                    {
                        eprintln!(
                            "[DEBUG VARIANT] {}:{} read={} flags(unmap/sec/supp/dup)={}/{}/{}/{} qpos={} read_seq={} -> idx={} ref={} alt={} snp={}",
                            chrom,
                            pos1,
                            String::from_utf8_lossy(&qname),
                            record.is_unmapped(),
                            record.is_secondary(),
                            record.is_supplementary(),
                            record.is_duplicate(),
                            qpos,
                            read_allele,
                            enc_idx,
                            region.ref_allele,
                            region.alt_allele,
                            region.is_snp()
                        );
                    }
                }
            }
        }

        // Log summary if records were skipped (prevents silent data loss)
        if skipped_records > 0 {
            eprintln!(
                "[WARN] Skipped {} corrupted BAM record(s) on {} (shown first {})",
                skipped_records, chrom, MAX_SKIP_WARNINGS.min(skipped_records)
            );
        }

        Ok(counts)
    }

    /// Group regions by chromosome while preserving encounter order
    fn group_regions_by_chrom(&self, regions: &[Region]) -> Vec<(String, Vec<(usize, Region)>)> {
        let mut grouped: Vec<Vec<(usize, Region)>> = Vec::new();
        let mut chrom_order: Vec<String> = Vec::new();
        let mut chrom_index: FxHashMap<String, usize> = FxHashMap::default();

        for (idx, region) in regions.iter().enumerate() {
            if let Some(&i) = chrom_index.get(&region.chrom) {
                grouped[i].push((idx, region.clone()));
            } else {
                let i = grouped.len();
                chrom_index.insert(region.chrom.clone(), i);
                chrom_order.push(region.chrom.clone());
                grouped.push(vec![(idx, region.clone())]);
            }
        }

        chrom_order.into_iter().zip(grouped).collect()
    }
}

/// Get base at genomic position, accounting for CIGAR operations
/// Matches WASP2 behavior: NO quality filtering by default
#[allow(dead_code)]
fn get_base_at_position(
    record: &bam::Record,
    target_pos: u32, // 0-based genomic position
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

/// Parse optional debug sites from env var WASP2_DEBUG_SNP (format: chr:pos or chr:pos:limit, comma-separated)
fn parse_debug_sites() -> FxHashMap<(String, u32), usize> {
    let mut map = FxHashMap::default();
    if let Ok(val) = std::env::var("WASP2_DEBUG_SNP") {
        for tok in val.split(',') {
            let tok = tok.trim();
            if tok.is_empty() {
                continue;
            }
            let parts: Vec<&str> = tok.split(':').collect();
            if parts.len() < 2 {
                continue;
            }
            let chrom = parts[0].to_string();
            if let Ok(pos) = parts[1].parse::<u32>() {
                let limit = if parts.len() >= 3 {
                    parts[2].parse::<usize>().unwrap_or(10)
                } else {
                    10
                };
                map.insert((chrom, pos), limit);
            }
        }
    }
    map
}
#[cfg(test)]
mod tests {
    use super::{BamCounter, Region};

    #[test]
    fn groups_regions_by_chrom_preserving_order() {
        let counter = BamCounter {
            bam_path: "dummy.bam".to_string(),
        };
        let regions = vec![
            Region {
                chrom: "chr1".into(),
                pos: 10,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
            Region {
                chrom: "chr1".into(),
                pos: 20,
                ref_allele: "C".into(),
                alt_allele: "T".into(),
            },
            Region {
                chrom: "chr2".into(),
                pos: 5,
                ref_allele: "T".into(),
                alt_allele: "C".into(),
            },
        ];

        let grouped = counter.group_regions_by_chrom(&regions);
        assert_eq!(grouped.len(), 2, "expected two chromosome groups");
        assert_eq!(grouped[0].0, "chr1");
        assert_eq!(grouped[1].0, "chr2");
        assert_eq!(grouped[0].1.len(), 2);
        assert_eq!(grouped[1].1.len(), 1);
        // Order preserved
        assert_eq!(grouped[0].1[0].1.pos, 10);
        assert_eq!(grouped[0].1[1].1.pos, 20);
    }
}
