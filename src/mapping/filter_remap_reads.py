import os
import timeit
from typing import Optional

import pysam

try:
    from wasp2_rust import filter_bam_wasp
except ImportError:
    filter_bam_wasp = None

def filt_remapped_reads(
    to_remap_bam: str,
    remapped_bam: str,
    filt_out_bam: str,
    keep_read_file: Optional[str] = None,
    use_rust: bool = True,
    threads: int = 1,
) -> None:
    rust_allowed = (
        use_rust
        and filter_bam_wasp is not None
        and os.environ.get("WASP2_DISABLE_RUST") != "1"
    )

    if not rust_allowed:
        raise RuntimeError(
            "Rust WASP filter not available. Build the extension with "
            "`maturin develop --release` in the WASP2 env."
        )

    # Rust path only
    filter_bam_wasp(
        to_remap_bam,
        remapped_bam,
        filt_out_bam,
        keep_read_file=keep_read_file,
        threads=threads,
    )


def merge_filt_bam(keep_bam: str, remapped_filt_bam: str, out_bam: str) -> None:
    
    start_time = timeit.default_timer()
    
    # Merge for for complete filt bam
    pysam.merge("-f", "-o", out_bam, keep_bam, remapped_filt_bam, catch_stdout=False)
    print(f"Merged BAM in {timeit.default_timer() - start_time:.2f} seconds")
    
    start_sort = timeit.default_timer()
    pysam.sort(out_bam, "-o", out_bam, catch_stdout=False)
    pysam.index(out_bam, catch_stdout=False)
    
    print(f"Sorted and Indexed BAM in {timeit.default_timer() - start_sort:.2f} seconds")
    
    # print(f"\nWrote merged WASP filtered BAM to...\n{out_bam}")
