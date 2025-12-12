use pyo3::prelude::*;
use rust_htslib::bam::{self, Read, Writer};
use rustc_hash::FxHashSet;

/// Buffered record info for paired-read handling
struct BufferedRead {
    pos: i64,
    mpos: i64,
}

/// Parsed WASP name components
/// Supports both old format (4 parts) and new extended format (6 parts)
#[derive(Debug, Clone)]
struct WaspNameInfo {
    orig_name: String,
    pos1: i64,
    pos2: i64,
    seq_num: i64,
    total_seqs: i64,
    /// Trim combination ID (hap_idx * 1000 + combo_idx), 0 for old format
    trim_combo_id: u16,
    /// Total number of trim combinations, 1 for old format
    total_combos: u16,
    /// Expected position shift tolerance per mate (absolute delta of indels)
    delta1: i64,
    delta2: i64,
}

/// Parse WASP-encoded name into components
/// Supports both old format: {name}_WASP_{pos1}_{pos2}_{seq}_{total}
/// And new format: {name}_WASP_{pos1}_{pos2}_{seq}_{total}_{trim_combo}_{total_combos}
fn parse_wasp_name(qname: &[u8]) -> Option<WaspNameInfo> {
    let split_idx = qname.windows(6).position(|w| w == b"_WASP_")?;

    let orig_name = std::str::from_utf8(&qname[..split_idx]).ok()?.to_string();
    let suffix = &qname[split_idx + 6..];
    let parts: Vec<&[u8]> = suffix.split(|b| *b == b'_').collect();

    let parse_i64 =
        |bytes: &[u8]| -> Option<i64> { std::str::from_utf8(bytes).ok()?.parse::<i64>().ok() };

    if parts.len() >= 8 {
        // New extended format with trim combinations and per-mate deltas
        Some(WaspNameInfo {
            orig_name,
            pos1: parse_i64(parts[0])?,
            pos2: parse_i64(parts[1])?,
            seq_num: parse_i64(parts[2])?,
            total_seqs: parse_i64(parts[3])?,
            trim_combo_id: parse_i64(parts[4])? as u16,
            total_combos: parse_i64(parts[5])? as u16,
            delta1: parse_i64(parts[6])?.abs(),
            delta2: parse_i64(parts[7])?.abs(),
        })
    } else if parts.len() >= 6 {
        // New extended format with trim combinations
        Some(WaspNameInfo {
            orig_name,
            pos1: parse_i64(parts[0])?,
            pos2: parse_i64(parts[1])?,
            seq_num: parse_i64(parts[2])?,
            total_seqs: parse_i64(parts[3])?,
            trim_combo_id: parse_i64(parts[4])? as u16,
            total_combos: parse_i64(parts[5])? as u16,
            delta1: 0,
            delta2: 0,
        })
    } else if parts.len() >= 4 {
        // Old format - backwards compatibility
        Some(WaspNameInfo {
            orig_name,
            pos1: parse_i64(parts[0])?,
            pos2: parse_i64(parts[1])?,
            seq_num: parse_i64(parts[2])?,
            total_seqs: parse_i64(parts[3])?,
            trim_combo_id: 0,
            total_combos: 1,
            delta1: 0,
            delta2: 0,
        })
    } else {
        None
    }
}

/// WASP-aware remap filter:
/// - Reads the remapped BAM with `_WASP_`-encoded names
/// - Buffers records until both mates of a pair arrive (like Python's paired_read_gen)
/// - Keeps pairs that returned to their original positions and saw all expected copies
/// - Writes a filtered BAM from the original `to_remap_bam` containing only kept read names
/// Returns (kept_reads, removed_moved, removed_missing)
#[pyfunction]
#[pyo3(signature = (to_remap_bam, remapped_bam, remap_keep_bam, keep_read_file=None, threads=1, same_locus_slop=0, expected_sidecar=None))]
pub fn filter_bam_wasp(
    to_remap_bam: String,
    remapped_bam: String,
    remap_keep_bam: String,
    keep_read_file: Option<String>,
    threads: usize,
    same_locus_slop: i64,
    expected_sidecar: Option<String>,
) -> PyResult<(u64, u64, u64)> {
    // Allow env override when Python binding lacks expected_sidecar kwarg
    let expected_sidecar = expected_sidecar.or_else(|| {
        std::env::var("WASP2_EXPECTED_SIDECAR")
            .ok()
            .map(|s| if s.is_empty() { None } else { Some(s) })
            .flatten()
    });

    // Optional sidecar of expected positions keyed by full qname
    let expected_map: Option<std::collections::HashMap<String, (i64, i64)>> =
        if let Some(sidecar_path) = expected_sidecar.as_ref() {
            let mut map = std::collections::HashMap::new();
            let content = std::fs::read_to_string(sidecar_path).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Failed to read sidecar {}: {}",
                    sidecar_path, e
                ))
            })?;
            for line in content.lines() {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() < 5 {
                    continue;
                }
                let q = parts[0].to_string();
                let p1 = parts[1].parse::<i64>().unwrap_or(0);
                let p2 = parts[2].parse::<i64>().unwrap_or(0);
                map.insert(q, (p1, p2));
            }
            Some(map)
        } else {
            None
        };

    // Track expected positions and remaining remapped copies
    let mut keep_set: FxHashSet<String> = FxHashSet::default();
    let mut pos_map: std::collections::HashMap<String, (i64, i64)> =
        std::collections::HashMap::new();
    let mut remaining: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    let mut removed_moved: u64 = 0;
    // Track per-read dynamic slop (variant-aware deltas)
    let mut slop_map: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    // Optional expected positions from sidecar
    let mut expected_pos_map: std::collections::HashMap<String, (i64, i64)> =
        std::collections::HashMap::new();

    // Buffer for incomplete pairs: keyed by full qname (with WASP suffix)
    // This mimics Python's paired_read_gen which buffers until both mates arrive
    let mut read_buffer: std::collections::HashMap<String, BufferedRead> =
        std::collections::HashMap::new();

    let mut remapped_reader = bam::Reader::from_path(&remapped_bam).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open remapped BAM: {}", e))
    })?;
    if threads > 1 {
        let _ = remapped_reader.set_threads(threads);
    }

    for rec_res in remapped_reader.records() {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };
        if rec.is_unmapped()
            || !rec.is_proper_pair()
            || rec.is_secondary()
            || rec.is_supplementary()
        {
            continue;
        }

        let qname = rec.qname();

        // Parse WASP name using the new function (handles both old and extended formats)
        let wasp_info = match parse_wasp_name(qname) {
            Some(info) => info,
            None => continue,
        };

        let name = wasp_info.orig_name.clone();
        let pos1 = wasp_info.pos1;
        let pos2 = wasp_info.pos2;
        let total = wasp_info.total_seqs;
        let dyn_slop = if same_locus_slop > 0 {
            same_locus_slop
        } else {
            wasp_info.delta1.max(wasp_info.delta2)
        };

        // Full query name for buffering (includes WASP suffix to distinguish haplotype copies)
        let full_qname = match std::str::from_utf8(qname) {
            Ok(s) => s.to_owned(),
            Err(_) => continue,
        };

        // Buffer records until both mates arrive (like Python's paired_read_gen)
        let rec_pos = rec.pos();
        let mate_pos = rec.mpos();

        if !read_buffer.contains_key(&full_qname) {
            // First mate of this pair - buffer it and continue
            read_buffer.insert(
                full_qname,
                BufferedRead {
                    pos: rec_pos,
                    mpos: mate_pos,
                },
            );
            continue;
        }

        // Second mate arrived - now we have a complete pair, process it
        let _first_read = read_buffer.remove(&full_qname).unwrap();

        // Initialize tracking for this original read name if not seen
        if !pos_map.contains_key(&name) {
            pos_map.insert(name.clone(), (pos1, pos2));
            remaining.insert(name.clone(), total);
            keep_set.insert(name.clone());
            slop_map.insert(name.clone(), dyn_slop);
            if let Some(ref m) = expected_map {
                if let Some((e1, e2)) = m.get(&full_qname) {
                    expected_pos_map.insert(full_qname.clone(), (*e1, *e2));
                }
            }
        } else if !keep_set.contains(&name) {
            // Already marked as failed
            continue;
        }

        // Count down expected copies - once per PAIR (not per record)
        if let Some(rem) = remaining.get_mut(&name) {
            *rem -= 1;
        }

        // Check if the remapped position matches original coordinates (mate order agnostic)
        // For indels, allow slop tolerance to handle micro-homology shifts
        if let Some((expect_pos, expect_mate)) = pos_map.get(&name) {
            // Prefer expected positions from sidecar (variant-aware), else use slop
            if let Some((e1, e2)) = expected_pos_map.get(&full_qname) {
                // Require remap to land on expected coords (mate-order agnostic)
                if !((rec_pos == *e1 && mate_pos == *e2) || (rec_pos == *e2 && mate_pos == *e1)) {
                    keep_set.remove(&name);
                    removed_moved += 1;
                    continue;
                }
            } else {
                let slop = *slop_map.get(&name).unwrap_or(&same_locus_slop);
                let matches = if slop == 0 {
                    // Strict matching for SNPs
                    (rec_pos == *expect_pos && mate_pos == *expect_mate)
                        || (rec_pos == *expect_mate && mate_pos == *expect_pos)
                } else {
                    // Allow slop tolerance for indels
                    let pos_diff1 = (rec_pos - *expect_pos).abs();
                    let mate_diff1 = (mate_pos - *expect_mate).abs();
                    let pos_diff2 = (rec_pos - *expect_mate).abs();
                    let mate_diff2 = (mate_pos - *expect_pos).abs();

                    (pos_diff1 <= slop && mate_diff1 <= slop)
                        || (pos_diff2 <= slop && mate_diff2 <= slop)
                };

                if !matches {
                    keep_set.remove(&name);
                    removed_moved += 1;
                    continue;
                }
            }
        }

        // Drop bookkeeping if all expected pairs seen
        if let Some(rem) = remaining.get(&name) {
            if *rem <= 0 {
                remaining.remove(&name);
                pos_map.remove(&name);
            }
        }
    }

    // Remove reads with missing counts
    let missing_count = remaining.len() as u64;
    removed_moved += missing_count;
    if missing_count > 0 {
        for name in remaining.keys() {
            keep_set.remove(name);
        }
    }

    // Persist keep list if requested
    if let Some(path) = keep_read_file.as_ref() {
        let mut file = std::fs::File::create(path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to create keep_read_file: {}",
                e
            ))
        })?;
        for name in keep_set.iter() {
            use std::io::Write;
            file.write_all(name.as_bytes())
                .and_then(|_| file.write_all(b"\n"))
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                        "Failed to write keep_read_file: {}",
                        e
                    ))
                })?;
        }
    }

    // Write filtered BAM from original to_remap input
    let mut to_reader = bam::Reader::from_path(&to_remap_bam).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open to_remap BAM: {}", e))
    })?;
    if threads > 1 {
        let _ = to_reader.set_threads(threads);
    }
    let header = bam::Header::from_template(to_reader.header());
    let mut writer =
        Writer::from_path(&remap_keep_bam, &header, bam::Format::Bam).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to create remap_keep_bam: {}",
                e
            ))
        })?;
    if threads > 1 {
        let _ = writer.set_threads(threads);
    }

    let mut kept_written: u64 = 0;
    for rec_res in to_reader.records() {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };
        let qname = match std::str::from_utf8(rec.qname()) {
            Ok(s) => s,
            Err(_) => continue,
        };
        if keep_set.contains(qname) {
            writer.write(&rec).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Write failed: {}", e))
            })?;
            kept_written += 1;
        }
    }

    Ok((kept_written, removed_moved, missing_count))
}
