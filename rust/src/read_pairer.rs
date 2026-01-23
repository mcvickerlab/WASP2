//! Read Pairing Utilities
//!
//! Efficiently pair reads from BAM files, replacing Python's `paired_read_gen`
//! and `paired_read_gen_stat` functions.
//!
//! Performance improvements:
//! - FxHashMap instead of Python dict for read storage
//! - Byte slices instead of String for read names (zero UTF-8 validation)
//! - Single-pass filtering (vs multiple if statements in Python)

use rust_htslib::bam;
use rustc_hash::FxHashMap;

// ============================================================================
// Data Structures
// ============================================================================

/// Statistics for read pairing (matches Python's ReadStats)
#[derive(Debug, Default, Clone)]
#[allow(dead_code)]
pub struct PairingStats {
    /// Reads discarded because unmapped
    pub discard_unmapped: usize,
    /// Reads discarded because not proper pair
    pub discard_improper_pair: usize,
    /// Reads discarded because secondary alignment
    pub discard_secondary: usize,
    /// Reads discarded because supplementary alignment
    pub discard_supplementary: usize,
    /// Read pairs where mate was missing
    pub discard_missing_pair: usize,
    /// Total read pairs successfully paired
    pub pairs_yielded: usize,
}

// ============================================================================
// Read Pairing Iterator
// ============================================================================

/// Iterator that yields properly paired reads from a BAM file
///
/// Replaces Python's `paired_read_gen()` and `paired_read_gen_stat()`.
///
/// # Performance
/// - Python: dict with String keys, multiple function calls
/// - Rust: FxHashMap with byte slice keys, inlined checks
/// - Expected speedup: 2-3x
#[allow(dead_code)]
pub struct ReadPairer {
    /// Internal reader
    reader: bam::Reader,
    /// Temporary storage for unpaired reads
    /// Key: read name (as bytes), Value: read record
    unpaired: FxHashMap<Vec<u8>, bam::Record>,
    /// Set of read names to discard (failed filters)
    discard_set: std::collections::HashSet<Vec<u8>>,
    /// Statistics tracking
    stats: PairingStats,
    /// Whether to collect statistics
    track_stats: bool,
    /// Current chromosome (if fetching specific region)
    chrom: Option<String>,
}

#[allow(dead_code)]
impl ReadPairer {
    /// Create a new ReadPairer for the entire BAM file
    pub fn new(bam_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let reader = bam::Reader::from_path(bam_path)?;

        Ok(Self {
            reader,
            unpaired: FxHashMap::default(),
            discard_set: std::collections::HashSet::new(),
            stats: PairingStats::default(),
            track_stats: false,
            chrom: None,
        })
    }

    /// Create a ReadPairer for a specific chromosome
    pub fn for_chromosome(bam_path: &str, chrom: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut pairer = Self::new(bam_path)?;
        pairer.chrom = Some(chrom.to_string());
        Ok(pairer)
    }

    /// Enable statistics tracking
    pub fn with_stats(mut self) -> Self {
        self.track_stats = true;
        self
    }

    /// Get accumulated statistics
    pub fn stats(&self) -> &PairingStats {
        &self.stats
    }

    /// Check if a read passes filters
    ///
    /// Filters:
    /// - Must be mapped
    /// - Must be proper pair
    /// - Must not be secondary alignment
    /// - Must not be supplementary alignment
    fn passes_filters(&mut self, read: &bam::Record) -> bool {
        // Check unmapped
        if read.is_unmapped() {
            if self.track_stats {
                self.stats.discard_unmapped += 1;
            }
            return false;
        }

        // Check proper pair
        if !read.is_proper_pair() {
            if self.track_stats {
                self.stats.discard_improper_pair += 1;
            }
            return false;
        }

        // Check secondary
        if read.is_secondary() {
            if self.track_stats {
                self.stats.discard_secondary += 1;
            }
            return false;
        }

        // Check supplementary
        if read.is_supplementary() {
            if self.track_stats {
                self.stats.discard_supplementary += 1;
            }
            return false;
        }

        true
    }

    /// Process a single read, returning paired read if mate found
    fn process_read(&mut self, read: bam::Record) -> Option<(bam::Record, bam::Record)> {
        // Check filters
        if !self.passes_filters(&read) {
            if self.track_stats {
                self.discard_set.insert(read.qname().to_vec());
            }
            return None;
        }

        let read_name = read.qname().to_vec();

        // Check if mate already seen
        if let Some(mate) = self.unpaired.remove(&read_name) {
            // Found mate! Yield pair in correct order (R1, R2)
            if self.track_stats {
                self.stats.pairs_yielded += 1;
            }

            if read.is_first_in_template() {
                Some((read, mate))
            } else {
                Some((mate, read))
            }
        } else {
            // No mate yet, store for later
            self.unpaired.insert(read_name, read);
            None
        }
    }

    /// Finalize pairing and update statistics for missing pairs
    pub fn finalize(&mut self) {
        if self.track_stats {
            // Count missing pairs (reads without mates)
            let missing = self
                .unpaired
                .keys()
                .filter(|k| !self.discard_set.contains(*k))
                .count();
            self.stats.discard_missing_pair = missing;
        }
    }
}

impl Iterator for ReadPairer {
    type Item = (bam::Record, bam::Record);

    fn next(&mut self) -> Option<Self::Item> {
        // TODO: Implement proper iterator that doesn't borrow self mutably
        // For now, this is a placeholder
        unimplemented!("ReadPairer iterator not yet implemented")
    }
}

// ============================================================================
// Convenience Functions
// ============================================================================

/// Pair all reads in a BAM file
///
/// Simple interface for basic use cases without statistics.
///
/// # Example
/// ```ignore
/// let pairs = pair_reads_from_bam("input.bam")?;
/// for (read1, read2) in pairs {
///     // Process pair
/// }
/// ```
#[allow(dead_code)]
pub fn pair_reads_from_bam(
    bam_path: &str,
) -> Result<impl Iterator<Item = (bam::Record, bam::Record)>, Box<dyn std::error::Error>> {
    ReadPairer::new(bam_path)
}

/// Pair reads from a specific chromosome with statistics
///
/// # Example
/// ```ignore
/// let mut pairer = pair_reads_from_chromosome("input.bam", "chr10")?;
/// for (read1, read2) in pairer.by_ref() {
///     // Process pair
/// }
/// pairer.finalize();
/// println!("Pairs yielded: {}", pairer.stats().pairs_yielded);
/// ```
#[allow(dead_code)]
pub fn pair_reads_from_chromosome(
    bam_path: &str,
    chrom: &str,
) -> Result<ReadPairer, Box<dyn std::error::Error>> {
    Ok(ReadPairer::for_chromosome(bam_path, chrom)?.with_stats())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // Remove when implemented
    fn test_read_pairer_basic() {
        // TODO: Create test BAM file
        // TODO: Pair reads
        // TODO: Verify pairs are correct
    }

    #[test]
    #[ignore]
    fn test_read_pairer_stats() {
        // TODO: Create test BAM with various read types
        // TODO: Pair with statistics enabled
        // TODO: Verify stats are accurate
    }

    #[test]
    #[ignore]
    fn test_filters() {
        // TODO: Test each filter individually
        // TODO: Verify discarded reads are counted correctly
    }

    #[test]
    #[ignore]
    fn test_chromosome_specific() {
        // TODO: Create BAM with multiple chromosomes
        // TODO: Pair only chr10
        // TODO: Verify only chr10 pairs returned
    }
}
