"""Single-cell counting and command execution against independent Python inputs."""

import ast
import contextlib
import io
import subprocess
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pysam
import pytest
from scripts.compare_archive_analysis import compare_frames, load_archive, single_cell_counts
from typer.testing import CliRunner

from analysis.__main__ import app as analysis_app
from counting.__main__ import app as count_app
from counting.count_alleles_sc import CountStatsSC, _sparse_from_counts, count_bc_snp_alleles
from counting.run_counting_sc import WaspCountSC

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = "41354fc087404ba1fa99542803ee9cb9b7e0e6bf"


def archive_counter():
    namespace = {"bisect_left": bisect_left}
    for path, name in [
        ("src/counting/count_alleles.py", "find_read_aln_pos"),
        ("src/counting/count_alleles_sc.py", "count_bc_snp_alleles"),
    ]:
        source = subprocess.check_output(["git", "show", f"{ARCHIVE}:{path}"], cwd=ROOT)
        tree = ast.parse(source)
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        exec(compile(ast.Module(body=[function], type_ignores=[]), path, "exec"), namespace)
    return namespace["count_bc_snp_alleles"]


@pytest.mark.rust
@pytest.mark.parametrize(
    "scenario", ["overlap", "paired", "deletion", "skip", "barcodes", "missing_sequence"]
)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("threads", [1, 4])
def test_three_way_single_cell_counting(tmp_path, scenario, reverse, threads):
    import wasp2_rust

    bam_path = tmp_path / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": 1000} for c in ["chr1", "chr2"]],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for chrom in [0, 1]:
            for i in range(6):
                read = pysam.AlignedSegment(out.header)
                read.query_name = "pair" if scenario == "paired" and i < 2 else f"read{i}"
                read.reference_id, read.reference_start = chrom, 90 + i
                read.mapping_quality, read.flag = 60, 0
                if scenario == "paired" and i < 2:
                    read.flag = 65 if i == 0 else 129
                    read.next_reference_id = chrom
                    read.next_reference_start = 91 if i == 0 else 90
                read.cigarstring = {"deletion": "10M1D10M", "skip": "10M1N10M"}.get(scenario, "30M")
                length = 20 if scenario in ["deletion", "skip"] else 30
                read.query_sequence = ("A" if i % 2 == 0 else "G") * length
                if scenario == "missing_sequence" and i == 0:
                    read.query_sequence = None
                if scenario != "barcodes" or i != 0:
                    read.set_tag(
                        "CB", "unknown" if scenario == "barcodes" and i == 1 else ["A", "B"][i % 2]
                    )
                out.write(read)
    pysam.index(str(bam_path))
    snps = [(0, 101, "A", "G"), (1, 111, "A", "G")]
    if reverse:
        snps.reverse()
    barcodes = {"A": 0, "B": 1}
    expected = []
    for counter in [archive_counter(), count_bc_snp_alleles]:
        counts = CountStatsSC()
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for chrom_index, chrom in enumerate(["chr1", "chr2"]):
                chromosome_snps = [
                    (idx + 2 * chrom_index, pos, ref, alt) for idx, pos, ref, alt in snps
                ]
                counter(bam, barcodes, chrom, iter(chromosome_snps), counts)
        expected.append([dict(counts.ref_count), dict(counts.alt_count), dict(counts.other_count)])
    actual = [{}, {}, {}]
    regions = [
        (idx + 2 * chrom_index, chrom, pos, ref, alt)
        for chrom_index, chrom in enumerate(["chr1", "chr2"])
        for idx, pos, ref, alt in snps
    ]
    for snp, bc, ref, alt, other in wasp2_rust.BamCounterSC(str(bam_path)).count_bc_alleles(
        regions, barcodes, 0, threads
    ):
        for layer, value in zip(actual, [ref, alt, other], strict=True):
            if value:
                layer[snp, bc] = value
    assert expected[0] == expected[1] == actual


@pytest.mark.parametrize("value", [0, 65535, 65536, 70000, 2**32])
def test_sparse_counts_preserve_large_values_through_h5ad(tmp_path, value):
    counts = defaultdict(int, {(0, 0): value}) if value else defaultdict(int)
    matrix = _sparse_from_counts(counts, (1, 1))
    adata = ad.AnnData(X=matrix + matrix, layers={"ref": matrix, "alt": matrix})
    path = tmp_path / "counts.h5ad"
    adata.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert restored.layers["ref"][0, 0] == value
    assert restored.X[0, 0] == 2 * value


@pytest.mark.parametrize("form", ["string", "list", "mixed", "file", "file_list"])
def test_sample_arguments_normalize_consistently(tmp_path, form):
    path = tmp_path / "samples.txt"
    path.write_text("sample1\n sample2 \n\n")
    samples = {
        "string": "sample1, sample2",
        "list": ["sample1", "sample2"],
        "mixed": ["sample1,sample2"],
        "file": str(path),
        "file_list": [str(path)],
    }[form]
    options = WaspCountSC("reads.bam", "variants.vcf", "barcodes.tsv", None, samples=samples)
    assert options.samples == ["sample1", "sample2"]


@pytest.fixture
def synthetic_count_bundle(tmp_path):
    historical = single_cell_counts()
    positions = [100 + 100 * i for i in range(6)]
    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(
        str(bam),
        "wb",
        header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]},
    ) as out:
        for snp, pos in enumerate(positions):
            for cell in range(4):
                for allele, layer in [("A", "ref"), ("G", "alt")]:
                    for k in range(int(historical.layers[layer][snp, cell])):
                        read = pysam.AlignedSegment(out.header)
                        read.query_name = f"{snp}_{cell}_{allele}_{k}"
                        read.query_sequence, read.flag = allele, 0
                        read.reference_id, read.reference_start = 0, pos - 1
                        read.mapping_quality, read.cigarstring = 60, "1M"
                        read.set_tag("CB", f"cell{cell}")
                        out.write(read)
    pysam.index(str(bam))
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
        + "".join(
            f"chr1\t{pos}\t.\tA\tG\t60\tPASS\t.\tGT\t" + ["0|1", "1|0"][i % 2] + "\n"
            for i, pos in enumerate(positions)
        )
    )
    barcodes, groups, features = [
        tmp_path / name for name in ["barcodes.tsv", "groups.tsv", "features.bed"]
    ]
    barcodes.write_text("".join(f"cell{i}\n" for i in range(4)))
    groups.write_text("cell0\tA\ncell1\tA\ncell2\tB\ncell3\tB\n")
    features.write_text("chr1\t90\t210\tp0\nchr1\t290\t310\tp1\nchr1\t390\t610\tp2\n")
    return bam, vcf, barcodes, groups, features, historical


@pytest.mark.rust
@pytest.mark.parametrize("command", ["find-imbalance-sc", "compare-imbalance"])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
def test_bam_to_h5ad_to_analysis_matches_archive(
    tmp_path, synthetic_count_bundle, monkeypatch, command, phased, pc
):
    import shutil

    if not shutil.which("bedtools"):
        pytest.skip("bedtools required for feature intersection")
    bam, vcf, barcodes, groups, features, historical = synthetic_count_bundle
    # Metadata conversion must work without the undeclared optional PyArrow package.
    monkeypatch.setattr(
        pl.DataFrame, "to_pandas", lambda *a, **k: pytest.fail("PyArrow-dependent conversion")
    )
    counts = tmp_path / "counts.h5ad"
    result = CliRunner().invoke(
        count_app,
        [
            "count-variants-sc",
            str(bam),
            str(vcf),
            str(barcodes),
            "--samples",
            "sample1",
            "--feature",
            str(features),
            "-o",
            str(counts),
        ],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    measured = ad.read_h5ad(counts)
    assert measured.shape == (4, 6)
    assert measured.obs_names.tolist() == [f"cell{i}" for i in range(4)]
    for layer in ["ref", "alt"]:
        np.testing.assert_array_equal(
            measured.layers[layer].toarray(), historical.layers[layer].toarray().T
        )
    output = tmp_path / "results.tsv"
    result = CliRunner().invoke(
        analysis_app,
        [
            command,
            str(counts),
            str(groups),
            "--sample",
            "sample1",
            "--phased" if phased else "--unphased",
            "--pseudocount",
            str(pc),
            "-o",
            str(output),
        ],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    archive = load_archive()
    function = (
        archive["as_analysis_sc"].get_imbalance_sc
        if command == "find-imbalance-sc"
        else archive["compare_ai"].get_compared_imbalance
    )
    reference = measured.T.copy()
    reference.var["group"] = ["A", "A", "B", "B"]
    with contextlib.redirect_stdout(io.StringIO()):
        expected = function(
            reference,
            sample="sample1",
            groups=["A", "B"],
            phased=phased,
            pseudocount=pc,
            min_count=10,
        )
    assert expected
    for key, frame in expected.items():
        suffix = key if isinstance(key, str) else "_".join(key)
        observed = pd.read_csv(tmp_path / f"results_{suffix}.tsv", sep="\t")
        comparison = compare_frames(frame, observed)
        assert comparison["status"] == "PASS", comparison
