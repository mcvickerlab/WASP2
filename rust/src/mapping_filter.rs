use pyo3::prelude::*;
use rust_htslib::bam::{self, Read, Writer};
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader};

/// Marker for buffered first mate (position data is obtained from the BAM record directly)
struct BufferedRead;

struct ExpectedPos {
    pos1: i64,
    pos2: i64,
    slop: i64,
}

/// Minimal parsed WASP name components needed for filtering.
///
/// Supports:
/// - Old format: `{name}_WASP_{pos1}_{pos2}_{seq}_{total}`
/// - New format: `{name}_WASP_{pos1}_{pos2}_{seq}_{total}_{trim_combo}_{total_combos}`
/// - New+delta:  `{name}_WASP_{pos1}_{pos2}_{seq}_{total}_{trim_combo}_{total_combos}_{d1}_{d2}`
#[derive(Debug, Clone, Copy)]
struct WaspNameInfo<'a> {
    orig_name: &'a [u8],
    pos1: i64,
    pos2: i64,
    total_seqs: i64,
    /// Expected position shift tolerance per mate (absolute delta of indels)
    delta1: i64,
    delta2: i64,
}

fn parse_i64_ascii(bytes: &[u8]) -> Option<i64> {
    if bytes.is_empty() {
        return None;
    }
    let mut idx = 0;
    let mut neg = false;
    if bytes[0] == b'-' {
        neg = true;
        idx = 1;
    } else if bytes[0] == b'+' {
        idx = 1;
    }
    if idx >= bytes.len() {
        return None;
    }
    let mut val: i64 = 0;
    let mut seen_digit = false;
    for &b in &bytes[idx..] {
        if !(b'0'..=b'9').contains(&b) {
            break;
        }
        seen_digit = true;
        val = val.checked_mul(10)? + (b - b'0') as i64;
    }
    if !seen_digit {
        return None;
    }
    Some(if neg { -val } else { val })
}

/// Parse WASP-encoded name into components
/// Supports both old format: {name}_WASP_{pos1}_{pos2}_{seq}_{total}
/// And new format: {name}_WASP_{pos1}_{pos2}_{seq}_{total}_{trim_combo}_{total_combos}
fn parse_wasp_name(qname: &[u8]) -> Option<WaspNameInfo<'_>> {
    let split_idx = qname.windows(6).position(|w| w == b"_WASP_")?;

    let orig_name = &qname[..split_idx];
    let suffix = &qname[split_idx + 6..];
    let mut parts = suffix.split(|&b| b == b'_');

    let pos1 = parse_i64_ascii(parts.next()?)?;
    let pos2 = parse_i64_ascii(parts.next()?)?;
    // seq_num is not needed by the filter
    let _seq_num = parts.next()?;
    let total_seqs = parse_i64_ascii(parts.next()?)?;

    // Optional fields
    let _trim_combo = parts.next();
    let _total_combos = parts.next();
    let delta1 = parts
        .next()
        .and_then(parse_i64_ascii)
        .map(|v| v.abs())
        .unwrap_or(0);
    let delta2 = parts
        .next()
        .and_then(parse_i64_ascii)
        .map(|v| v.abs())
        .unwrap_or(0);

    Some(WaspNameInfo {
        orig_name,
        pos1,
        pos2,
        total_seqs,
        delta1,
        delta2,
    })
}

/// Check if remapped positions match expected positions (mate-order agnostic)
fn positions_match(rec_pos: i64, mate_pos: i64, exp_pos1: i64, exp_pos2: i64, slop: i64) -> bool {
    if slop < 0 {
        eprintln!(
            "[WARN] positions_match: negative slop ({}), clamping to 0",
            slop
        );
    }
    let slop = slop.max(0);
    if slop == 0 {
        (rec_pos == exp_pos1 && mate_pos == exp_pos2)
            || (rec_pos == exp_pos2 && mate_pos == exp_pos1)
    } else {
        let pos_diff1 = (rec_pos - exp_pos1).abs();
        let mate_diff1 = (mate_pos - exp_pos2).abs();
        let pos_diff2 = (rec_pos - exp_pos2).abs();
        let mate_diff2 = (mate_pos - exp_pos1).abs();

        (pos_diff1 <= slop && mate_diff1 <= slop) || (pos_diff2 <= slop && mate_diff2 <= slop)
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
            .filter(|s| !s.is_empty())
    });

    // Optional sidecar of expected positions keyed by full qname.
    // Stored as bytes to avoid per-read UTF-8/String allocations in the hot loop.
    let expected_map: Option<FxHashMap<Vec<u8>, (i64, i64)>> =
        if let Some(sidecar_path) = expected_sidecar.as_ref() {
            let file = std::fs::File::open(sidecar_path).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Failed to open sidecar {}: {}",
                    sidecar_path, e
                ))
            })?;
            let mut reader = BufReader::new(file);
            let mut buf: Vec<u8> = Vec::new();
            let mut map: FxHashMap<Vec<u8>, (i64, i64)> = FxHashMap::default();

            loop {
                buf.clear();
                let n = reader.read_until(b'\n', &mut buf).map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                        "Failed to read sidecar {}: {}",
                        sidecar_path, e
                    ))
                })?;
                if n == 0 {
                    break;
                }
                if buf.ends_with(b"\n") {
                    buf.pop();
                    if buf.ends_with(b"\r") {
                        buf.pop();
                    }
                }

                let mut parts = buf.split(|&b| b == b'\t');
                let q = match parts.next() {
                    Some(v) if !v.is_empty() => v,
                    _ => continue,
                };
                let p1 = match parts.next().and_then(parse_i64_ascii) {
                    Some(v) => v,
                    None => continue,
                };
                let p2 = match parts.next().and_then(parse_i64_ascii) {
                    Some(v) => v,
                    None => continue,
                };
                // Keep compatibility with older sidecars: require at least 5 columns (q, p1, p2, ...)
                if parts.next().is_none() || parts.next().is_none() {
                    continue;
                }
                map.insert(q.to_vec(), (p1, p2));
            }
            Some(map)
        } else {
            None
        };

    // Track expected positions and remaining remapped copies
    let mut keep_set: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut pos_map: FxHashMap<Vec<u8>, ExpectedPos> = FxHashMap::default();
    let mut remaining: FxHashMap<Vec<u8>, i64> = FxHashMap::default();
    let mut removed_moved: u64 = 0;

    // Buffer for incomplete pairs: keyed by full qname (with WASP suffix)
    // This mimics Python's paired_read_gen which buffers until both mates arrive
    let mut read_buffer: FxHashMap<Vec<u8>, BufferedRead> = FxHashMap::default();

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

        let name = wasp_info.orig_name;
        let pos1 = wasp_info.pos1;
        let pos2 = wasp_info.pos2;
        let total = wasp_info.total_seqs;
        let dyn_slop = if same_locus_slop > 0 {
            same_locus_slop
        } else {
            wasp_info.delta1.max(wasp_info.delta2)
        };

        // Buffer records until both mates arrive (like Python's paired_read_gen)
        let rec_pos = rec.pos();
        let mate_pos = rec.mpos();

        if !read_buffer.contains_key(qname) {
            // First mate of this pair - buffer it and continue
            read_buffer.insert(qname.to_vec(), BufferedRead);
            continue;
        }

        // Second mate arrived - now we have a complete pair, process it
        let _first_read = read_buffer.remove(qname).unwrap();

        // Initialize tracking for this original read name if not seen
        if !pos_map.contains_key(name) {
            let owned_name = name.to_vec();
            pos_map.insert(
                owned_name.clone(),
                ExpectedPos {
                    pos1,
                    pos2,
                    slop: dyn_slop,
                },
            );
            remaining.insert(owned_name.clone(), total);
            keep_set.insert(owned_name);
        } else if !keep_set.contains(name) {
            // Already marked as failed
            continue;
        }

        // Count down expected copies - once per PAIR (not per record)
        if let Some(rem) = remaining.get_mut(name) {
            *rem -= 1;
        }

        // Check if the remapped position matches original coordinates (mate order agnostic)
        // For indels, allow slop tolerance to handle micro-homology shifts
        if let Some(expect) = pos_map.get(name) {
            // Prefer expected positions from sidecar (variant-aware), else use slop
            let matches = expected_map
                .as_ref()
                .and_then(|m| m.get(qname))
                .map(|(e1, e2)| positions_match(rec_pos, mate_pos, *e1, *e2, 0))
                .unwrap_or_else(|| {
                    positions_match(rec_pos, mate_pos, expect.pos1, expect.pos2, expect.slop)
                });

            if !matches {
                keep_set.remove(name);
                removed_moved += 1;
                continue;
            }
        }

        // Drop bookkeeping if all expected pairs seen
        if let Some(rem) = remaining.get(name) {
            if *rem <= 0 {
                remaining.remove(name);
                pos_map.remove(name);
            }
        }
    }

    // Remove reads with missing counts (tracked separately from moved reads)
    let missing_count = remaining.len() as u64;
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
            file.write_all(name)
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
        if keep_set.contains(rec.qname()) {
            writer.write(&rec).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Write failed: {}", e))
            })?;
            kept_written += 1;
        }
    }

    Ok((kept_written, removed_moved, missing_count))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_i64_ascii() {
        assert_eq!(parse_i64_ascii(b"123"), Some(123));
        assert_eq!(parse_i64_ascii(b"-123"), Some(-123));
        assert_eq!(parse_i64_ascii(b"+123"), Some(123));
        assert_eq!(parse_i64_ascii(b"123/1"), Some(123));
        assert_eq!(parse_i64_ascii(b"/1"), None);
        assert_eq!(parse_i64_ascii(b""), None);
        assert_eq!(parse_i64_ascii(b"abc"), None);
    }

    #[test]
    fn test_parse_wasp_name_old_format_with_mate_suffix() {
        let qname = b"readX_WASP_100_200_1_10/1";
        let info = parse_wasp_name(qname).unwrap();
        assert_eq!(info.orig_name, b"readX");
        assert_eq!(info.pos1, 100);
        assert_eq!(info.pos2, 200);
        assert_eq!(info.total_seqs, 10);
        assert_eq!(info.delta1, 0);
        assert_eq!(info.delta2, 0);
    }

    #[test]
    fn test_parse_wasp_name_extended_without_delta() {
        let qname = b"readX_WASP_100_200_1_10_5_6/1";
        let info = parse_wasp_name(qname).unwrap();
        assert_eq!(info.orig_name, b"readX");
        assert_eq!(info.pos1, 100);
        assert_eq!(info.pos2, 200);
        assert_eq!(info.total_seqs, 10);
        assert_eq!(info.delta1, 0);
        assert_eq!(info.delta2, 0);
    }

    #[test]
    fn test_positions_match_exact_forward() {
        assert!(positions_match(100, 200, 100, 200, 0));
    }

    #[test]
    fn test_positions_match_exact_swapped() {
        assert!(positions_match(200, 100, 100, 200, 0));
    }

    #[test]
    fn test_positions_match_exact_mismatch() {
        assert!(!positions_match(100, 201, 100, 200, 0));
    }

    #[test]
    fn test_positions_match_slop_within() {
        assert!(positions_match(101, 199, 100, 200, 2));
    }

    #[test]
    fn test_positions_match_slop_swapped_within() {
        assert!(positions_match(199, 101, 100, 200, 2));
    }

    #[test]
    fn test_positions_match_slop_at_boundary() {
        assert!(positions_match(102, 198, 100, 200, 2));
    }

    #[test]
    fn test_positions_match_slop_past_boundary() {
        assert!(!positions_match(103, 200, 100, 200, 2));
    }

    #[test]
    fn test_positions_match_slop_mixed_fail() {
        // One within slop, the other outside
        assert!(!positions_match(101, 210, 100, 200, 2));
    }

    #[test]
    fn test_positions_match_negative_positions() {
        assert!(positions_match(-5, 10, -5, 10, 0));
        assert!(positions_match(10, -5, -5, 10, 0));
    }

    #[test]
    fn test_positions_match_negative_slop_clamped() {
        // Negative slop is clamped to 0, so exact match is required
        assert!(positions_match(100, 200, 100, 200, -5));
        assert!(!positions_match(101, 200, 100, 200, -5));
    }

    #[test]
    fn test_parse_wasp_name_extended_with_delta() {
        let qname = b"readX_WASP_100_200_1_10_5_6_2_3/1";
        let info = parse_wasp_name(qname).unwrap();
        assert_eq!(info.orig_name, b"readX");
        assert_eq!(info.pos1, 100);
        assert_eq!(info.pos2, 200);
        assert_eq!(info.total_seqs, 10);
        assert_eq!(info.delta1, 2);
        assert_eq!(info.delta2, 3);
    }
}
