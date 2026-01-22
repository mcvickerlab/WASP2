use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use rust_htslib::bam::{self, header::HeaderRecord, Format, Header, Writer};
use std::collections::HashMap;

/// Create a synthetic BAM file for benchmarking
fn create_test_bam(path: &str, n_reads: usize, include_wasp_suffix: bool) -> std::io::Result<()> {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", &"1.6");
    hd.push_tag(b"SO", &"coordinate");
    header.push_record(&hd);

    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"chr1");
    sq.push_tag(b"LN", &"248956422");
    header.push_record(&sq);

    let mut pg = HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", &"test");
    pg.push_tag(b"PN", &"benchmark");
    pg.push_tag(b"VN", &"1.0");
    header.push_record(&pg);
    let mut writer = Writer::from_path(path, &header, Format::Bam).unwrap();

    for i in 0..n_reads {
        let mut record = bam::Record::new();

        // Create read name with WASP suffix for remapped BAM
        let qname = if include_wasp_suffix {
            format!(
                "read_{}_WASP_{}_{}_{}_2",
                i,
                1000 + i * 100,
                1300 + i * 100,
                i % 10
            )
        } else {
            format!("read_{}", i)
        };

        record.set_qname(qname.as_bytes());
        record.set_tid(0); // chr1
        record.set_pos(1000 + i as i64 * 100);
        record.set_mpos(1300 + i as i64 * 100);
        record.set_mapq(60);
        record.set_flags(99); // Proper pair, first read
        record.set_insert_size(450);

        let seq = b"ATCGATCGATCGATCGATCGATCG";
        let qual = vec![30u8; seq.len()];

        // Set qname/seq/qual/cigar together (rust-htslib 0.44 API)
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(seq.len() as u32)]);
        record.set(qname.as_bytes(), Some(&cigar), seq, &qual);

        writer.write(&record).unwrap();
    }

    Ok(())
}

/// Benchmark the WASP name parsing (hottest part)
fn bench_qname_parsing(c: &mut Criterion) {
    let test_names: Vec<&[u8]> = vec![
        b"read_1_WASP_1000_1300_5_10".as_ref(),
        b"read_2_WASP_2000_2300_3_8".as_ref(),
        b"read_3_WASP_3000_3300_7_12".as_ref(),
        b"very_long_read_name_12345_WASP_4000_4300_2_15".as_ref(),
    ];

    c.bench_function("qname_wasp_parse", |b| {
        b.iter(|| {
            for qname in &test_names {
                // Simulate the WASP parsing from mapping_filter.rs
                let split_idx = black_box(qname).windows(6).position(|w| w == b"_WASP_");
                if let Some(idx) = split_idx {
                    let suffix = &qname[idx + 6..];
                    let parts: Vec<&[u8]> = suffix.split(|b| *b == b'_').collect();

                    if parts.len() >= 4 {
                        // Parse positions
                        let _ = std::str::from_utf8(parts[0])
                            .ok()
                            .and_then(|s| s.parse::<i64>().ok());
                        let _ = std::str::from_utf8(parts[1])
                            .ok()
                            .and_then(|s| s.parse::<i64>().ok());
                        let _ = std::str::from_utf8(parts[3])
                            .ok()
                            .and_then(|s| s.parse::<i64>().ok());
                    }
                }
            }
        });
    });
}

/// Benchmark position comparison logic
fn bench_position_matching(c: &mut Criterion) {
    let test_cases = vec![
        ((1000i64, 1300i64), (1000i64, 1300i64), 0i64), // Exact match
        ((1000i64, 1300i64), (1002i64, 1298i64), 5i64), // Within slop
        ((1000i64, 1300i64), (1010i64, 1310i64), 5i64), // Outside slop
    ];

    c.bench_function("position_matching", |b| {
        b.iter(|| {
            for (rec_pos, expect_pos, slop) in &test_cases {
                let (rec_p, rec_m) = rec_pos;
                let (exp_p, exp_m) = expect_pos;

                let _ = if *slop == 0 {
                    (*rec_p == *exp_p && *rec_m == *exp_m) || (*rec_p == *exp_m && *rec_m == *exp_p)
                } else {
                    let pos_diff1 = (*rec_p - *exp_p).abs();
                    let mate_diff1 = (*rec_m - *exp_m).abs();
                    let pos_diff2 = (*rec_p - *exp_m).abs();
                    let mate_diff2 = (*rec_m - *exp_p).abs();

                    (pos_diff1 <= *slop && mate_diff1 <= *slop)
                        || (pos_diff2 <= *slop && mate_diff2 <= *slop)
                };
            }
        });
    });
}

/// Benchmark HashMap operations (keeping track of read names)
fn bench_hashmap_operations(c: &mut Criterion) {
    use rustc_hash::FxHashSet;

    let mut group = c.benchmark_group("hashmap_ops");

    for size in [100, 1000, 10000].iter() {
        group.bench_with_input(BenchmarkId::new("insert_lookup", size), size, |b, &size| {
            b.iter(|| {
                let mut keep_set: FxHashSet<String> = FxHashSet::default();
                let mut pos_map: HashMap<String, (i64, i64)> = HashMap::new();

                for i in 0..size {
                    let name = format!("read_{}", i);
                    keep_set.insert(name.clone());
                    pos_map.insert(name, (1000 + i as i64, 1300 + i as i64));
                }

                // Lookup
                for i in 0..size {
                    let name = format!("read_{}", i);
                    let _ = black_box(keep_set.contains(&name));
                    let _ = black_box(pos_map.get(&name));
                }
            });
        });
    }
    group.finish();
}

/// Benchmark String allocation in hot loop
fn bench_string_allocation(c: &mut Criterion) {
    let qname_bytes = b"read_12345";

    let mut group = c.benchmark_group("string_alloc");

    group.bench_function("string_from_utf8_owned", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = black_box(std::str::from_utf8(qname_bytes).ok().map(|s| s.to_owned()));
            }
        });
    });

    group.bench_function("string_from_utf8_borrowed", |b| {
        b.iter(|| {
            for _ in 0..1000 {
                let _ = black_box(std::str::from_utf8(qname_bytes).ok());
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_qname_parsing,
    bench_position_matching,
    bench_hashmap_operations,
    bench_string_allocation
);
criterion_main!(benches);
