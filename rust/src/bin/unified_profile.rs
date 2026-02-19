use anyhow::{Context, Result};
use std::path::PathBuf;
use wasp2_rust::{unified_make_reads, unified_make_reads_parallel, UnifiedConfig};

fn parse_arg(flag: &str) -> Option<String> {
    let mut args = std::env::args();
    while let Some(a) = args.next() {
        if a == flag {
            return args.next();
        }
    }
    None
}

fn parse_usize(flag: &str, default: usize) -> usize {
    parse_arg(flag)
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() -> Result<()> {
    let bam = parse_arg("--bam").context("Missing --bam")?;
    let bed = parse_arg("--bed").context("Missing --bed")?;
    let out_dir = PathBuf::from(
        parse_arg("--out-dir").unwrap_or_else(|| "/tmp/wasp2_unified_profile".to_string()),
    );

    let threads = parse_usize("--threads", 8);
    let max_seqs = parse_usize("--max-seqs", 64);
    let channel_buffer = parse_usize("--channel-buffer", 50_000);
    let compression_threads = parse_usize("--compression-threads", 1);
    let compress_output = parse_arg("--compress-output")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    let parallel = parse_arg("--parallel")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(true);
    let indel_mode = parse_arg("--indel-mode")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false);

    std::fs::create_dir_all(&out_dir).context("Failed to create --out-dir")?;
    let r1 = out_dir.join("remap_r1.fq");
    let r2 = out_dir.join("remap_r2.fq");

    let config = UnifiedConfig {
        read_threads: threads,
        max_seqs,
        pair_buffer_reserve: 100_000,
        channel_buffer,
        compression_threads,
        compress_output,
        indel_mode,
        max_indel_size: 50,
        keep_no_flip_names_path: None,
        remap_names_path: None,
    };

    let run = || {
        if parallel {
            unified_make_reads_parallel(
                &bam,
                &bed,
                r1.to_string_lossy().as_ref(),
                r2.to_string_lossy().as_ref(),
                &config,
            )
        } else {
            unified_make_reads(
                &bam,
                &bed,
                r1.to_string_lossy().as_ref(),
                r2.to_string_lossy().as_ref(),
                &config,
            )
        }
    };

    // Match the Python binding behavior: use a per-run thread pool so we can control
    // Rayon worker threads precisely (e.g. for profiling).
    let stats = if parallel && threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("Failed to build Rayon thread pool")?;
        pool.install(run)?
    } else {
        run()?
    };

    eprintln!(
        "done: total_reads={} pairs={} haps={}",
        stats.total_reads, stats.pairs_processed, stats.haplotypes_written
    );
    Ok(())
}
