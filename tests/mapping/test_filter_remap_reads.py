from mapping import filter_remap_reads


def test_merge_filt_bam_uses_portable_samtools_syntax(monkeypatch):
    calls = []

    def record_run(args, *, check):
        calls.append((args, check))

    monkeypatch.setattr(filter_remap_reads.subprocess, "run", record_run)

    filter_remap_reads.merge_filt_bam("keep.bam", "remapped_filtered.bam", "output.bam", threads=3)

    assert calls == [
        (
            [
                "samtools",
                "merge",
                "-@",
                "3",
                "-f",
                "output.bam",
                "keep.bam",
                "remapped_filtered.bam",
            ],
            True,
        ),
        (["samtools", "index", "-@", "3", "output.bam"], True),
    ]
