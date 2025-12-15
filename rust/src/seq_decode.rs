use rust_htslib::bam;

// Matches rust-htslib's internal decode table for BAM 4-bit base encoding.
// See: rust-htslib bam/record.rs `DECODE_BASE`.
const DECODE_BASE: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Decode a BAM record's 4-bit encoded sequence into `out`.
///
/// This avoids the heavy `read.seq().as_bytes()` allocation by reusing `out`.
pub fn decode_seq_into(read: &bam::Record, out: &mut Vec<u8>) {
    let seq = read.seq();
    let len = seq.len();
    let encoded = seq.encoded;

    out.clear();
    out.resize(len, 0);

    // Decode two bases per packed byte (high then low nibble).
    for (i, packed) in encoded.iter().copied().enumerate() {
        let pos = i * 2;
        if pos >= len {
            break;
        }
        let hi = (packed >> 4) as usize;
        out[pos] = DECODE_BASE[hi];
        let pos2 = pos + 1;
        if pos2 < len {
            let lo = (packed & 0x0F) as usize;
            out[pos2] = DECODE_BASE[lo];
        }
    }
}

/// Copy a BAM record's qualities into `out` (reusing the allocation).
pub fn copy_qual_into(read: &bam::Record, out: &mut Vec<u8>) {
    let qual = read.qual();
    out.clear();
    out.extend_from_slice(qual);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn make_record(seq: &[u8], qual: &[u8]) -> bam::Record {
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        let mut rec = bam::Record::new();
        rec.set(b"q1", Some(&cigar), seq, qual);
        rec.set_pos(100);
        rec
    }

    #[test]
    fn decode_seq_into_matches_rust_htslib() {
        let seq = b"ACGTNACGTN";
        let qual = vec![10u8; seq.len()];
        let rec = make_record(seq, &qual);

        let mut buf = Vec::new();
        decode_seq_into(&rec, &mut buf);
        assert_eq!(buf, rec.seq().as_bytes());

        // Reuse the buffer with a different length.
        let rec2 = make_record(b"NNNN", &[1, 2, 3, 4]);
        decode_seq_into(&rec2, &mut buf);
        assert_eq!(buf, rec2.seq().as_bytes());
    }

    #[test]
    fn copy_qual_into_matches_rust_htslib() {
        let seq = b"ACGTN";
        let qual = vec![0u8, 1, 2, 40, 41];
        let rec = make_record(seq, &qual);

        let mut buf = Vec::new();
        copy_qual_into(&rec, &mut buf);
        assert_eq!(buf, rec.qual().to_vec());
    }
}
