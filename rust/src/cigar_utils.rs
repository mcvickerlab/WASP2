//! CIGAR-aware position mapping utilities for INDEL support
//!
#![allow(dead_code)] // Utility functions for future optimization paths
//!
//! This module provides efficient reference-to-query position mapping using
//! rust-htslib's `aligned_pairs_full()` API, which matches pysam's
//! `get_aligned_pairs(matches_only=False)`.
//!
//! # Key Concepts
//!
//! When a read has insertions or deletions in its CIGAR string, the simple
//! arithmetic `query_pos = ref_pos - read_start` is WRONG. We need to account
//! for CIGAR operations that consume reference vs query bases differently.
//!
//! ## CIGAR Operations
//! - M/=/X: consume both ref and query (1:1 mapping)
//! - I: consume query only (insertion in read)
//! - D/N: consume ref only (deletion/skip in read)
//! - S: consume query only (soft clip)
//! - H: consume neither (hard clip)
//!
//! ## Position Mapping for Indels
//!
//! For a deletion in the read (ref bases with no query bases), we need TWO mappings:
//! - `ref2query_left`: maps ref_pos to the LAST query position BEFORE the deletion
//! - `ref2query_right`: maps ref_pos to the FIRST query position AFTER the deletion
//!
//! This allows proper slicing: use left for variant start, right for variant end.
//!
//! # Performance
//!
//! - `aligned_pairs_full()` is O(n) where n = alignment length
//! - Building maps is O(n) with two passes
//! - Single position lookup via `find_query_position()` is O(k) where k = CIGAR ops
//!
//! For reads with few variants, targeted lookup is faster than building full maps.

use anyhow::Result;
use rust_htslib::bam::{self, ext::BamRecordExtensions};
use rustc_hash::FxHashMap;

/// Position mapping result for a reference position
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryPosition {
    /// Exact match: ref position maps to this query position
    Mapped(usize),
    /// Deletion: ref position is deleted, use flanking positions
    Deleted { left_flank: usize, right_flank: usize },
    /// Not covered: ref position is outside the alignment
    NotCovered,
}

/// Build reference-to-query position mappings using rust-htslib's aligned_pairs_full
///
/// This is the Rust equivalent of Python's:
/// ```python
/// pairs = read.get_aligned_pairs(matches_only=False)
/// ```
///
/// # Returns
/// Two HashMaps:
/// - `ref2query_left`: For each ref position, the nearest LEFT query position
/// - `ref2query_right`: For each ref position, the nearest RIGHT query position
///
/// For matched positions, both maps return the same value.
/// For deletions, left gives the position BEFORE, right gives the position AFTER.
///
/// # Performance
/// O(n) where n = alignment length. Builds ~n entries in each map.
/// Consider using `find_query_position()` for single lookups.
pub fn build_ref2query_maps(
    read: &bam::Record,
) -> (FxHashMap<i64, usize>, FxHashMap<i64, usize>) {
    let mut ref2query_left: FxHashMap<i64, usize> = FxHashMap::default();
    let mut ref2query_right: FxHashMap<i64, usize> = FxHashMap::default();

    // Collect aligned pairs: [Option<query_pos>, Option<ref_pos>]
    // - Both Some: matched base
    // - query=Some, ref=None: insertion
    // - query=None, ref=Some: deletion
    let pairs: Vec<[Option<i64>; 2]> = read.aligned_pairs_full().collect();

    if pairs.is_empty() {
        return (ref2query_left, ref2query_right);
    }

    // Forward pass: build left mapping
    let mut last_query_pos: Option<usize> = None;
    for pair in &pairs {
        let query_pos = pair[0];
        let ref_pos = pair[1];

        if let Some(rp) = ref_pos {
            if let Some(qp) = query_pos {
                // Matched base
                ref2query_left.insert(rp, qp as usize);
                last_query_pos = Some(qp as usize);
            } else {
                // Deletion: use last known query position (left flank)
                if let Some(lqp) = last_query_pos {
                    ref2query_left.insert(rp, lqp);
                }
            }
        } else if let Some(qp) = query_pos {
            // Insertion: just update last_query_pos
            last_query_pos = Some(qp as usize);
        }
    }

    // Backward pass: build right mapping
    let mut next_query_pos: Option<usize> = None;
    for pair in pairs.iter().rev() {
        let query_pos = pair[0];
        let ref_pos = pair[1];

        if let Some(rp) = ref_pos {
            if let Some(qp) = query_pos {
                // Matched base
                ref2query_right.insert(rp, qp as usize);
                next_query_pos = Some(qp as usize);
            } else {
                // Deletion: use next known query position (right flank)
                if let Some(nqp) = next_query_pos {
                    ref2query_right.insert(rp, nqp);
                }
            }
        } else if let Some(qp) = query_pos {
            // Insertion: just update next_query_pos
            next_query_pos = Some(qp as usize);
        }
    }

    (ref2query_left, ref2query_right)
}

/// Find query position for a single reference position by walking CIGAR
///
/// This is more efficient than building full maps when you only need 1-4 lookups.
///
/// # Arguments
/// * `read` - BAM record
/// * `target_ref_pos` - Reference position to find (0-based)
///
/// # Returns
/// - `Some(query_pos)` if the position is mapped
/// - `None` if the position is in a deletion or outside alignment
///
/// # Performance
/// O(k) where k = number of CIGAR operations (typically <10)
pub fn find_query_position(read: &bam::Record, target_ref_pos: i64) -> Option<usize> {
    use rust_htslib::bam::record::Cigar;

    let cigar = read.cigar();
    let mut query_pos: usize = 0;
    let mut ref_pos = read.pos();

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                // Check if target is in this match block
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + (*len as i64) {
                    let offset = (target_ref_pos - ref_pos) as usize;
                    return Some(query_pos + offset);
                }
                query_pos += *len as usize;
                ref_pos += *len as i64;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                // Only query advances
                query_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                // Only reference advances - position is in deletion
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + (*len as i64) {
                    return None; // Position is deleted
                }
                ref_pos += *len as i64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                // No advancement
            }
        }
    }

    None // Position not found
}

/// Find query position with flanking information for deletions
///
/// Enhanced version that returns flanking positions for deleted bases.
///
/// # Returns
/// - `QueryPosition::Mapped(pos)` - exact mapping
/// - `QueryPosition::Deleted { left, right }` - position is deleted, use flanks
/// - `QueryPosition::NotCovered` - position outside alignment
pub fn find_query_position_with_flanks(
    read: &bam::Record,
    target_ref_pos: i64,
) -> QueryPosition {
    use rust_htslib::bam::record::Cigar;

    let cigar = read.cigar();
    let mut query_pos: usize = 0;
    let mut ref_pos = read.pos();
    let mut last_query_pos: usize = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + (*len as i64) {
                    let offset = (target_ref_pos - ref_pos) as usize;
                    return QueryPosition::Mapped(query_pos + offset);
                }
                query_pos += *len as usize;
                ref_pos += *len as i64;
                last_query_pos = query_pos.saturating_sub(1);
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                query_pos += *len as usize;
                last_query_pos = query_pos.saturating_sub(1);
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if target_ref_pos >= ref_pos && target_ref_pos < ref_pos + (*len as i64) {
                    // Position is in deletion - return flanking positions
                    return QueryPosition::Deleted {
                        left_flank: last_query_pos,
                        right_flank: query_pos, // Next query position after deletion
                    };
                }
                ref_pos += *len as i64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    QueryPosition::NotCovered
}

/// Apply allele substitution to a sequence with CIGAR awareness
///
/// This handles:
/// - SNPs: simple base replacement
/// - Deletions: remove bases from sequence
/// - Insertions: add bases to sequence
///
/// # Arguments
/// * `seq` - Original read sequence
/// * `qual` - Original quality scores
/// * `ref_start` - Variant reference start position (0-based)
/// * `ref_end` - Variant reference end position (exclusive, 0-based)
/// * `ref_allele` - Reference allele string
/// * `alt_allele` - Alternate allele to substitute
/// * `ref2query_left` - Left position mapping (for variant start)
/// * `ref2query_right` - Right position mapping (for variant end)
///
/// # Returns
/// (new_sequence, new_quality) with substitution applied
pub fn apply_cigar_aware_substitution(
    seq: &[u8],
    qual: &[u8],
    ref_start: i64,
    ref_end: i64,
    ref_allele: &str,
    alt_allele: &str,
    ref2query_left: &FxHashMap<i64, usize>,
    ref2query_right: &FxHashMap<i64, usize>,
) -> Result<(Vec<u8>, Vec<u8>)> {
    // Get query positions using appropriate mappings
    let query_start = ref2query_left
        .get(&ref_start)
        .copied()
        .ok_or_else(|| anyhow::anyhow!("Ref position {} not in left map", ref_start))?;

    // For end position, we want the position AFTER the last ref base
    // ref_end is exclusive, so we look up ref_end - 1 and add 1
    let query_end = ref2query_right
        .get(&(ref_end - 1))
        .map(|&p| p + 1)
        .ok_or_else(|| anyhow::anyhow!("Ref position {} not in right map", ref_end - 1))?;

    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();

    // Build new sequence
    let mut new_seq = Vec::with_capacity(seq.len() + alt_len.saturating_sub(ref_len));
    let mut new_qual = Vec::with_capacity(qual.len() + alt_len.saturating_sub(ref_len));

    // Part before variant
    new_seq.extend_from_slice(&seq[..query_start]);
    new_qual.extend_from_slice(&qual[..query_start]);

    // Substitute allele
    new_seq.extend_from_slice(alt_allele.as_bytes());

    // Handle quality scores for the substituted region
    if alt_len == ref_len {
        // Same length: use original qualities
        if query_end <= qual.len() {
            new_qual.extend_from_slice(&qual[query_start..query_end]);
        }
    } else if alt_len < ref_len {
        // Deletion: truncate qualities
        let qual_to_copy = alt_len.min(query_end.saturating_sub(query_start));
        if query_start + qual_to_copy <= qual.len() {
            new_qual.extend_from_slice(&qual[query_start..query_start + qual_to_copy]);
        }
    } else {
        // Insertion: copy original quals + fill extra with default Q30
        let orig_qual_len = query_end.saturating_sub(query_start).min(qual.len() - query_start);
        if query_start + orig_qual_len <= qual.len() {
            new_qual.extend_from_slice(&qual[query_start..query_start + orig_qual_len]);
        }
        let extra_needed = alt_len.saturating_sub(orig_qual_len);
        new_qual.extend(std::iter::repeat(30u8).take(extra_needed));
    }

    // Part after variant
    if query_end < seq.len() {
        new_seq.extend_from_slice(&seq[query_end..]);
    }
    if query_end < qual.len() {
        new_qual.extend_from_slice(&qual[query_end..]);
    }

    Ok((new_seq, new_qual))
}

/// Check if any variants in a list are indels (different ref/alt lengths)
pub fn has_indels(variants: &[(i64, i64, &str, &str)]) -> bool {
    variants.iter().any(|(_, _, ref_allele, alt_allele)| {
        ref_allele.len() != alt_allele.len()
    })
}

/// Segment a sequence based on variant positions
///
/// Returns segments suitable for haplotype generation:
/// - Even indices (0, 2, 4, ...): non-variant regions
/// - Odd indices (1, 3, 5, ...): variant regions to be swapped
///
/// # Arguments
/// * `seq` - Original sequence
/// * `qual` - Original quality scores
/// * `variant_positions` - List of (query_start, query_end) positions
///
/// # Returns
/// (seq_segments, qual_segments) where segments alternate between
/// non-variant and variant regions
pub fn segment_sequence(
    seq: &[u8],
    qual: &[u8],
    variant_positions: &[(usize, usize)],
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut seq_segments = Vec::new();
    let mut qual_segments = Vec::new();
    let mut last_end = 0;

    for &(start, end) in variant_positions {
        // Non-variant segment before this variant
        seq_segments.push(seq[last_end..start].to_vec());
        qual_segments.push(qual[last_end..start].to_vec());

        // Variant segment
        seq_segments.push(seq[start..end].to_vec());
        qual_segments.push(qual[start..end].to_vec());

        last_end = end;
    }

    // Final non-variant segment
    seq_segments.push(seq[last_end..].to_vec());
    qual_segments.push(qual[last_end..].to_vec());

    (seq_segments, qual_segments)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_position_enum() {
        let mapped = QueryPosition::Mapped(42);
        let deleted = QueryPosition::Deleted {
            left_flank: 10,
            right_flank: 11,
        };
        let not_covered = QueryPosition::NotCovered;

        assert_eq!(mapped, QueryPosition::Mapped(42));
        assert_eq!(
            deleted,
            QueryPosition::Deleted {
                left_flank: 10,
                right_flank: 11
            }
        );
        assert_eq!(not_covered, QueryPosition::NotCovered);
    }

    #[test]
    fn test_has_indels_snp_only() {
        let variants = vec![
            (100, 101, "A", "G"),
            (200, 201, "C", "T"),
        ];
        let variants_ref: Vec<(i64, i64, &str, &str)> = variants
            .iter()
            .map(|(s, e, r, a)| (*s as i64, *e as i64, *r, *a))
            .collect();
        assert!(!has_indels(&variants_ref));
    }

    #[test]
    fn test_has_indels_with_deletion() {
        let variants = vec![
            (100, 101, "A", "G"),      // SNP
            (200, 203, "ACG", "A"),    // Deletion
        ];
        let variants_ref: Vec<(i64, i64, &str, &str)> = variants
            .iter()
            .map(|(s, e, r, a)| (*s as i64, *e as i64, *r, *a))
            .collect();
        assert!(has_indels(&variants_ref));
    }

    #[test]
    fn test_has_indels_with_insertion() {
        let variants = vec![
            (100, 101, "A", "ACGT"),   // Insertion
        ];
        let variants_ref: Vec<(i64, i64, &str, &str)> = variants
            .iter()
            .map(|(s, e, r, a)| (*s as i64, *e as i64, *r, *a))
            .collect();
        assert!(has_indels(&variants_ref));
    }

    #[test]
    fn test_segment_sequence() {
        let seq = b"AAAAABBBBBCCCCC";
        let qual = vec![30u8; 15];
        let positions = vec![(5, 10)]; // Variant at positions 5-10

        let (seq_segs, qual_segs) = segment_sequence(seq, &qual, &positions);

        assert_eq!(seq_segs.len(), 3); // before, variant, after
        assert_eq!(seq_segs[0], b"AAAAA"); // before
        assert_eq!(seq_segs[1], b"BBBBB"); // variant
        assert_eq!(seq_segs[2], b"CCCCC"); // after

        assert_eq!(qual_segs.len(), 3);
        assert_eq!(qual_segs[0].len(), 5);
        assert_eq!(qual_segs[1].len(), 5);
        assert_eq!(qual_segs[2].len(), 5);
    }

    #[test]
    fn test_segment_sequence_multiple_variants() {
        let seq = b"AAABBBCCCDDDEEE";
        let qual = vec![30u8; 15];
        let positions = vec![(3, 6), (9, 12)]; // Two variants

        let (seq_segs, _qual_segs) = segment_sequence(seq, &qual, &positions);

        assert_eq!(seq_segs.len(), 5); // before, var1, between, var2, after
        assert_eq!(seq_segs[0], b"AAA");
        assert_eq!(seq_segs[1], b"BBB");
        assert_eq!(seq_segs[2], b"CCC");
        assert_eq!(seq_segs[3], b"DDD");
        assert_eq!(seq_segs[4], b"EEE");
    }
}
