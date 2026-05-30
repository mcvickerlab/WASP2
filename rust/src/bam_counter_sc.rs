use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read as BamRead};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::HashMap;
use std::path::Path;

use crate::cigar_utils::find_query_position;

/// Per-cell (single-cell) BAM allele counter using rust-htslib with batched fetching.
///
/// Mirrors the bulk `BamCounter::process_chromosome_reads` structure, but produces
/// per-(SNP, barcode) allele counts instead of aggregate per-SNP counts. The Python
/// reference being mirrored is `count_bc_snp_alleles` in
/// `src/counting/count_alleles_sc.py`.
///
/// Parity invariants (must hold vs. the Python reference):
///   * per-chromosome qname dedup (a read is counted at most once per chromosome)
///   * each read is assigned to the EARLIEST-overlap SNP by encounter index
///   * SNP single-base classification via `record.seq()[qpos]` (A/C/G/T/N)
///   * barcode read from the `CB` aux tag; reads with no CB are skipped
///   * barcodes absent from `bc_dict` are skipped (matches Python)
///   * `min_qual` defaults to 0 (no quality filtering)
///   * counts are u32
#[pyclass]
pub struct BamCounterSC {
    bam_path: String,
}

#[derive(Debug, Clone)]
struct ScRegion {
    snp_idx: u32, // stable SNP index (the `index` column from the Python snp_df)
    chrom: String,
    pos: u32,           // 1-based position from Python
    ref_allele: String, // Full reference allele
    alt_allele: String, // Full alternate allele
}

impl ScRegion {
    /// Returns true if this variant is a simple SNP (single base change)
    fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() == 1
    }
}

// PyO3 expands #[pymethods] into impl blocks that trigger non_local_definitions warnings;
// suppress the noise until we restructure.
#[allow(non_local_definitions)]
#[pymethods]
impl BamCounterSC {
    #[new]
    fn new(bam_path: String) -> PyResult<Self> {
        if !Path::new(&bam_path).exists() {
            return Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
                format!("BAM file not found: {}", bam_path),
            ));
        }

        Ok(BamCounterSC { bam_path })
    }

    /// Count alleles at variant positions per cell barcode using batched fetching.
    ///
    /// Args:
    ///     regions: List of (snp_idx, chrom, pos, ref, alt) tuples
    ///     bc_dict: Dict mapping cell barcode string -> integer index
    ///     min_qual: Minimum base quality (default: 0 for WASP2 compatibility)
    ///     threads: Number of worker threads (default: 1). >1 enables Rayon per-chromosome parallelism.
    ///
    /// Returns:
    ///     List of (snp_idx, bc_idx, ref_count, alt_count, other_count) COO tuples (nonzero only)
    #[pyo3(signature = (regions, bc_dict, min_qual=0, threads=1))]
    fn count_bc_alleles(
        &self,
        py: Python<'_>,
        regions: &Bound<'_, PyList>,
        bc_dict: &Bound<'_, PyDict>,
        min_qual: u8,
        threads: usize,
    ) -> PyResult<Vec<(u32, u32, u32, u32, u32)>> {
        // Parse Python regions
        let mut rust_regions = Vec::with_capacity(regions.len());
        for item in regions.iter() {
            let tuple = item.cast::<PyTuple>()?;
            let snp_idx: u32 = tuple.get_item(0)?.extract()?;
            let chrom: String = tuple.get_item(1)?.extract()?;
            let pos: u32 = tuple.get_item(2)?.extract()?;
            let ref_allele: String = tuple.get_item(3)?.extract()?;
            let alt_allele: String = tuple.get_item(4)?.extract()?;

            if ref_allele.is_empty() {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Empty ref_allele for variant at {}:{}",
                    chrom, pos
                )));
            }
            if alt_allele.is_empty() {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Empty alt_allele for variant at {}:{}",
                    chrom, pos
                )));
            }

            rust_regions.push(ScRegion {
                snp_idx,
                chrom,
                pos,
                ref_allele,
                alt_allele,
            });
        }

        // Parse barcode dict into a plain Rust HashMap so it can cross the GIL boundary.
        let mut bc_map: HashMap<String, usize> = HashMap::with_capacity(bc_dict.len());
        for (k, v) in bc_dict.iter() {
            let bc: String = k.extract()?;
            let idx: usize = v.extract()?;
            bc_map.insert(bc, idx);
        }

        py.detach(|| self.count_bc_alleles_impl(&rust_regions, &bc_map, min_qual, threads))
    }
}

impl BamCounterSC {
    fn count_bc_alleles_impl(
        &self,
        regions: &[ScRegion],
        bc_dict: &HashMap<String, usize>,
        min_qual: u8,
        threads: usize,
    ) -> PyResult<Vec<(u32, u32, u32, u32, u32)>> {
        // Group regions by chromosome while preserving encounter order
        let grouped = self.group_regions_by_chrom(regions);

        // Per-chromosome partial COO maps keyed by (snp_idx, bc_idx)
        let mut combined: FxHashMap<(u32, usize), (u32, u32, u32)> = FxHashMap::default();

        if threads > 1 {
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
                    let partial_results: Result<Vec<_>, _> = grouped
                        .par_iter()
                        .map(|(chrom, chrom_regions)| {
                            self.process_chromosome_reads(chrom, chrom_regions, bc_dict, min_qual)
                        })
                        .collect();

                    for partial in partial_results? {
                        for (key, (r, a, o)) in partial {
                            let entry = combined.entry(key).or_insert((0, 0, 0));
                            entry.0 += r;
                            entry.1 += a;
                            entry.2 += o;
                        }
                    }
                    Ok::<(), PyErr>(())
                })?;
        } else {
            for (chrom, chrom_regions) in &grouped {
                let partial =
                    self.process_chromosome_reads(chrom, chrom_regions, bc_dict, min_qual)?;
                for (key, (r, a, o)) in partial {
                    let entry = combined.entry(key).or_insert((0, 0, 0));
                    entry.0 += r;
                    entry.1 += a;
                    entry.2 += o;
                }
            }
        }

        // Emit COO: only nonzero entries (per-allele counts are always >0 by construction here,
        // since we only insert on an actual classification).
        let mut out: Vec<(u32, u32, u32, u32, u32)> = Vec::with_capacity(combined.len());
        for ((snp_idx, bc_idx), (r, a, o)) in combined {
            out.push((snp_idx, bc_idx as u32, r, a, o));
        }

        Ok(out)
    }

    /// Process a single chromosome SNP-major with a per-SNP narrow fetch.
    ///
    /// This is a direct mirror of the Python reference `count_bc_snp_alleles`
    /// (`src/counting/count_alleles_sc.py`):
    ///
    /// ```text
    /// read_set = set()                       # global across this chromosome
    /// for idx, pos, ref, alt in snp_list:    # SNP-major, INPUT ORDER
    ///     for read in bam.fetch(chrom, pos-1, pos):   # narrow 1bp window
    ///         if read.query_name in read_set: continue
    ///         read_bc = read.get_tag("CB") else continue
    ///         if read_bc not in bc_dict: continue
    ///         qpos = find_read_aln_pos(read, pos-1)   # matches_only semantics
    ///         if qpos is None: continue               # NO read_set insert (deletion case)
    ///         base = seq[qpos]
    ///         classify ref/alt/other; counts[(idx, bc)] += 1
    ///         read_set.add(read.query_name)           # mark seen ONLY after success
    /// ```
    ///
    /// Because the SNP loop is in input order and a qname is added to `read_set`
    /// only after a successful classification, each qname is assigned to the
    /// EARLIEST classifiable SNP across all of its alignments — identical to
    /// Python. The earlier read-major scan + per-qname best-map is no longer
    /// needed and is removed. Crucially, when `find_query_position` returns
    /// `None` (deletion / non-match), we `continue` WITHOUT inserting the qname
    /// into `seen`, so a later SNP can still classify it (byte-identity).
    ///
    /// Returns a map keyed by (snp_idx, bc_idx) -> (ref, alt, other).
    fn process_chromosome_reads(
        &self,
        chrom: &str,
        regions: &[(usize, ScRegion)],
        bc_dict: &HashMap<String, usize>,
        min_qual: u8,
    ) -> PyResult<FxHashMap<(u32, usize), (u32, u32, u32)>> {
        let mut bam = bam::IndexedReader::from_path(&self.bam_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM: {}", e))
        })?;

        // Multithreaded BGZF decode. Defaults to 2 background decode threads, which
        // accelerates the load-imbalance tail (when only the largest chromosomes
        // remain and most Rayon workers are idle): whole-genome 60.4s -> 40.5s in
        // benchmarking, bit-identical. Set WASP2_RUST_BGZF_THREADS=0 to disable.
        // Only changes decode threading, never record order/bytes.
        let bgzf_threads: usize = std::env::var("WASP2_RUST_BGZF_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2);
        if bgzf_threads > 0 {
            let _ = bam.set_threads(bgzf_threads);
        }

        let mut counts: FxHashMap<(u32, usize), (u32, u32, u32)> = FxHashMap::default();

        if regions.is_empty() {
            return Ok(counts);
        }

        // Global per-chromosome qname dedup set, mirroring Python's `read_set`.
        let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();

        // Reuse ONE Record buffer across all fetches on this chromosome to avoid
        // per-read allocations. `record` persists across `bam.fetch(...)` calls.
        let mut record = bam::record::Record::new();

        let mut skipped_records: u32 = 0;
        const MAX_SKIP_WARNINGS: u32 = 5;

        // SNP-major: iterate variants in INPUT ORDER (the encounter index is the
        // stable order from Python's snp_list). `_enc_idx` is unused now because
        // earliest-SNP assignment is implied by the SNP loop order itself.
        for (_enc_idx, region) in regions.iter() {
            // region.pos is 1-based; fetch the 0-based half-open [pos-1, pos) window.
            let start = (region.pos.saturating_sub(1)) as i64;
            let end = region.pos as i64;
            if let Err(e) = bam.fetch((chrom, start, end)) {
                eprintln!(
                    "[WARN] Failed to fetch {}:{}-{}: {}. Skipping variant.",
                    chrom, start, end, e
                );
                continue;
            }

            // 0-based reference position used for query-position lookup.
            let target_ref_pos = start;

            while let Some(res) = bam.read(&mut record) {
                match res {
                    Ok(_) => {}
                    Err(e) => {
                        skipped_records += 1;
                        if skipped_records <= MAX_SKIP_WARNINGS {
                            eprintln!("[WARN] Skipped corrupted BAM record on {}: {}", chrom, e);
                        }
                        continue;
                    }
                }
                if record.is_unmapped() {
                    continue;
                }

                // Borrowed qname (&[u8], no allocation); dedup BEFORE any allocation.
                let qname = record.qname();
                // Already counted this qname at an earlier SNP on this chromosome.
                if seen.contains(qname) {
                    continue;
                }

                // Read the CB (cell barcode) aux tag; skip reads without one.
                // Borrowed &str (no allocation); map to index then drop the borrow.
                let cb = match record.aux(b"CB") {
                    Ok(rust_htslib::bam::record::Aux::String(s)) => s,
                    _ => continue,
                };
                // Map barcode -> index via bc_dict; skip barcodes not present (matches Python).
                // bc_idx is a copied usize, so the `cb`/qname borrows end here.
                let bc_idx = match bc_dict.get(cb) {
                    Some(i) => *i,
                    None => continue,
                };

                // matches_only query-position lookup; None for deletions / non-match.
                // CRITICAL: on None we continue WITHOUT inserting into `seen`.
                let qpos = match find_query_position(&record, target_ref_pos) {
                    Some(q) => q,
                    None => continue,
                };

                let quals = record.qual();
                if min_qual > 0 && (qpos >= quals.len() || quals[qpos] < min_qual) {
                    continue;
                }

                let seq = record.seq();
                if qpos >= seq.len() {
                    continue;
                }

                // Single-base classification via seq[qpos], mirroring Python's
                // `seq[qpos] == ref / == alt / else other`. scATAC SNPs are
                // single-base ref/alt, so this is the byte-identical path.
                let base: &[u8] = match seq[qpos] {
                    b'A' => b"A",
                    b'C' => b"C",
                    b'G' => b"G",
                    b'T' => b"T",
                    b'N' => b"N",
                    _ => b"",
                };

                let entry = counts.entry((region.snp_idx, bc_idx)).or_insert((0, 0, 0));
                if base == region.ref_allele.as_bytes() {
                    entry.0 += 1;
                } else if base == region.alt_allele.as_bytes() {
                    entry.1 += 1;
                } else {
                    entry.2 += 1;
                }

                // Mark seen ONLY after a successful classification (owning copy on success).
                seen.insert(record.qname().to_vec());
            }
        }

        if skipped_records > 0 {
            eprintln!(
                "[WARN] Skipped {} corrupted BAM record(s) on {} (shown first {})",
                skipped_records,
                chrom,
                MAX_SKIP_WARNINGS.min(skipped_records)
            );
        }

        Ok(counts)
    }

    /// Group regions by chromosome while preserving encounter order.
    fn group_regions_by_chrom(
        &self,
        regions: &[ScRegion],
    ) -> Vec<(String, Vec<(usize, ScRegion)>)> {
        let mut grouped: Vec<Vec<(usize, ScRegion)>> = Vec::new();
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

#[cfg(test)]
mod tests {
    use super::{BamCounterSC, ScRegion};

    #[test]
    fn groups_sc_regions_by_chrom_preserving_order() {
        let counter = BamCounterSC {
            bam_path: "dummy.bam".to_string(),
        };
        let regions = vec![
            ScRegion {
                snp_idx: 0,
                chrom: "chr1".into(),
                pos: 10,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
            ScRegion {
                snp_idx: 1,
                chrom: "chr1".into(),
                pos: 20,
                ref_allele: "C".into(),
                alt_allele: "T".into(),
            },
            ScRegion {
                snp_idx: 2,
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
        // Encounter order preserved; snp_idx travels with the region.
        assert_eq!(grouped[0].1[0].1.pos, 10);
        assert_eq!(grouped[0].1[0].1.snp_idx, 0);
        assert_eq!(grouped[0].1[1].1.pos, 20);
        assert_eq!(grouped[0].1[1].1.snp_idx, 1);
        assert_eq!(grouped[1].1[0].1.snp_idx, 2);
    }

    #[test]
    fn sc_region_is_snp_classification() {
        let snp = ScRegion {
            snp_idx: 0,
            chrom: "chr1".into(),
            pos: 1,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
        };
        assert!(snp.is_snp());

        let indel = ScRegion {
            snp_idx: 1,
            chrom: "chr1".into(),
            pos: 1,
            ref_allele: "A".into(),
            alt_allele: "AG".into(),
        };
        assert!(!indel.is_snp());
    }
}
