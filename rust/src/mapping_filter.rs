use pyo3::prelude::*;
use rust_htslib::bam::{self, Read, Writer};
use rustc_hash::FxHashSet;

/// WASP-aware remap filter:
/// - Reads the remapped BAM with `_WASP_`-encoded names
/// - Keeps pairs that returned to their original positions and saw all expected copies
/// - Writes a filtered BAM from the original `to_remap_bam` containing only kept read names
/// Returns (kept_reads, removed_moved, removed_missing)
#[pyfunction]
#[pyo3(signature = (to_remap_bam, remapped_bam, remap_keep_bam, keep_read_file=None, threads=1, same_locus_slop=0))]
pub fn filter_bam_wasp(
    to_remap_bam: String,
    remapped_bam: String,
    remap_keep_bam: String,
    keep_read_file: Option<String>,
    threads: usize,
    same_locus_slop: i64,
) -> PyResult<(u64, u64, u64)> {
    // Track expected positions and remaining remapped copies
    let mut keep_set: FxHashSet<String> = FxHashSet::default();
    let mut pos_map: std::collections::HashMap<String, (i64, i64)> = std::collections::HashMap::new();
    let mut remaining: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    let mut removed_moved: u64 = 0;

    let mut remapped_reader = bam::Reader::from_path(&remapped_bam)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open remapped BAM: {}", e)))?;
    if threads > 1 {
        let _ = remapped_reader.set_threads(threads);
    }

    for rec_res in remapped_reader.records() {
        let rec = match rec_res {
            Ok(r) => r,
            Err(_) => continue,
        };
        if rec.is_unmapped() || !rec.is_proper_pair() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }

        let qname = rec.qname();
        let split_idx = qname.windows(6).position(|w| w == b"_WASP_");
        let split_idx = match split_idx {
            Some(i) => i,
            None => continue,
        };

        let orig_name = &qname[..split_idx];
        let suffix = &qname[split_idx + 6..];
        let parts: Vec<&[u8]> = suffix.split(|b| *b == b'_').collect();
        if parts.len() < 4 {
            continue;
        }

        let parse_i64 = |bytes: &[u8]| -> Option<i64> {
            std::str::from_utf8(bytes).ok()?.parse::<i64>().ok()
        };

        let pos1 = match parse_i64(parts[0]) {
            Some(v) => v,
            None => continue,
        };
        let pos2 = match parse_i64(parts[1]) {
            Some(v) => v,
            None => continue,
        };
        let total = match parse_i64(parts[3]) {
            Some(v) => v,
            None => continue,
        };

        let name = match std::str::from_utf8(orig_name) {
            Ok(s) => s.to_owned(),
            Err(_) => continue,
        };

        if !pos_map.contains_key(&name) {
            pos_map.insert(name.clone(), (pos1, pos2));
            remaining.insert(name.clone(), total);
            keep_set.insert(name.clone());
        } else if !keep_set.contains(&name) {
            continue;
        }

        // Count down expected copies
        if let Some(rem) = remaining.get_mut(&name) {
            *rem -= 1;
        }

        // Check if the remapped position matches original coordinates (mate order agnostic)
        // For indels, allow slop tolerance to handle micro-homology shifts
        let rec_pos = rec.pos();
        let mate_pos = rec.mpos();
        if let Some((expect_pos, expect_mate)) = pos_map.get(&name) {
            let matches = if same_locus_slop == 0 {
                // Strict matching for SNPs
                (rec_pos == *expect_pos && mate_pos == *expect_mate)
                    || (rec_pos == *expect_mate && mate_pos == *expect_pos)
            } else {
                // Allow slop tolerance for indels
                let pos_diff1 = (rec_pos - *expect_pos).abs();
                let mate_diff1 = (mate_pos - *expect_mate).abs();
                let pos_diff2 = (rec_pos - *expect_mate).abs();
                let mate_diff2 = (mate_pos - *expect_pos).abs();

                (pos_diff1 <= same_locus_slop && mate_diff1 <= same_locus_slop)
                    || (pos_diff2 <= same_locus_slop && mate_diff2 <= same_locus_slop)
            };

            if !matches {
                keep_set.remove(&name);
                remaining.remove(&name);
                removed_moved += 1;
                continue;
            }
        }

        // Drop bookkeeping if done
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
        let mut file = std::fs::File::create(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create keep_read_file: {}", e)))?;
        for name in keep_set.iter() {
            use std::io::Write;
            file.write_all(name.as_bytes())
                .and_then(|_| file.write_all(b"\n"))
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write keep_read_file: {}", e)))?;
        }
    }

    // Write filtered BAM from original to_remap input
    let mut to_reader = bam::Reader::from_path(&to_remap_bam)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open to_remap BAM: {}", e)))?;
    if threads > 1 {
        let _ = to_reader.set_threads(threads);
    }
    let header = bam::Header::from_template(to_reader.header());
    let mut writer = Writer::from_path(&remap_keep_bam, &header, bam::Format::Bam)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create remap_keep_bam: {}", e)))?;
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
            writer
                .write(&rec)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Write failed: {}", e)))?;
            kept_written += 1;
        }
    }

    Ok((kept_written, removed_moved, missing_count))
}
