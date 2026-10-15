"""RNA annotation and Rust-backed grouping compatibility with Python 41354fc."""

import gzip
import subprocess
import types
from pathlib import Path

import pandas as pd
import polars as pl
import pysam
import pytest
from scripts.compare_archive_analysis import adapt_linear_reference, compare_frames, load_archive
from typer.testing import CliRunner

from analysis.__main__ import app as analysis_app
from counting.count_alleles import make_count_df
from counting.parse_gene_data import make_gene_data, parse_gene_file

ROOT = Path(__file__).resolve().parents[2]


def _archive_annotations():
    source = subprocess.check_output(
        ["git", "show", "41354fc:src/counting/parse_gene_data.py"], cwd=ROOT, text=True
    )
    module = types.ModuleType("archived_gene_data")
    exec(compile(source, "archived_gene_data", "exec"), module.__dict__)
    return module


def _annotation(path, rows):
    text = "".join(
        f"chr1\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}\n"
        for feature, start, end, attrs in rows
    )
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as stream:
            stream.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("suffix", ["gtf", "gtf.gz", "gff3"])
@pytest.mark.parametrize("feature", [None, "exon", "transcript", "gene"])
def test_annotation_selection_matches_archive(tmp_path, suffix, feature):
    rows = []
    for kind in ("gene", "transcript", "exon"):
        attrs = (
            f'gene_id "g1"; transcript_id "t1"; {kind}_id "{kind}1";'
            if suffix.startswith("gtf")
            else f"ID={kind}1;Parent=g1;"
        )
        rows.append((kind, 91, 220, attrs))
    path = _annotation(tmp_path / f"annotation.{suffix}", rows)
    expected = _archive_annotations().parse_gene_file(path, feature=feature)
    actual = parse_gene_file(path, feature=feature)
    assert actual[1:] == expected[1:]
    assert actual.gene_df.to_dicts() == expected.gene_df.to_dicts()
    if feature is None:
        assert actual.feature == "exon"


@pytest.mark.parametrize(
    "features,expected", [(["gene", "transcript"], "transcript"), (["gene"], "gene")]
)
def test_default_feature_priority_without_exons(tmp_path, features, expected):
    path = _annotation(
        tmp_path / "annotation.gff3", [(name, 1, 20, f"ID={name}1") for name in features]
    )
    assert parse_gene_file(path).feature == expected


@pytest.mark.parametrize(
    "attrs,options,identifier,parent",
    [
        ("ID=e1;Parent=t1", {}, "e1", "t1"),
        ("ID=e1;Parent=t1;", {}, "e1", "t1"),
        ('exon_id "e1"; gene_id "g1"', {}, "e1", "g1"),
        ("OtherID=wrong;ID=e1;OtherParent=wrong;Parent=t1", {}, "e1", "t1"),
        ("OtherID=wrong;Name=e1", {}, "e1", "e1"),
        ("ID=e1;Parent=t1,t2", {}, "e1", "t1,t2"),
        (
            "custom.id=e1;custom.parent=g1",
            {"attribute": "custom.id", "parent_attribute": "custom.parent"},
            "e1",
            "g1",
        ),
        ("ID=e1", {}, "e1", "e1"),
    ],
)
def test_final_attributes_and_exact_keys(tmp_path, attrs, options, identifier, parent):
    path = _annotation(tmp_path / "annotation.gff3", [("exon", 91, 220, attrs)])
    result = parse_gene_file(path, **options)
    assert result.gene_df.row(0) == ("chr1", 90, 220, identifier, parent)


@pytest.mark.parametrize("missing", ["exon_id", "gene_id"])
@pytest.mark.parametrize("value", [None, "", " "])
def test_missing_identifiers_fail_before_bed_write(tmp_path, missing, value):
    attrs = {"exon_id": "e2", "gene_id": "g1"}
    if value is None:
        del attrs[missing]
    else:
        attrs[missing] = value
    path = _annotation(
        tmp_path / "annotation.gtf",
        [
            ("exon", 1, 20, 'exon_id "e1"; gene_id "g1";'),
            ("exon", 31, 50, " ".join(f'{key} "{val}";' for key, val in attrs.items())),
        ],
    )
    bed = tmp_path / "features.bed"
    with pytest.raises(ValueError, match=rf"{missing}.*1.*chr1:31-50"):
        make_gene_data(path, bed, attribute="exon_id", parent_attribute="gene_id")
    assert not bed.exists()


@pytest.mark.parametrize(
    "attrs,options,message",
    [
        ("OtherID=e1", {}, "No feature identifier"),
        ("ID=e1", {"feature": "gene"}, "No 'gene' features"),
        ("ID=e1", {"attribute": "absent"}, "absent"),
    ],
)
def test_invalid_annotation_errors(tmp_path, attrs, options, message):
    path = _annotation(tmp_path / "annotation.gff3", [("exon", 1, 20, attrs)])
    with pytest.raises(ValueError, match=message):
        parse_gene_file(path, **options)


@pytest.fixture
def rna_counts(tmp_path):
    """Spliced reads and paired mates overlap shared exons in two genes."""
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as stream:

        def write(name, start, sequence, cigar="1M", flag=0, mate=-1):
            read = pysam.AlignedSegment(stream.header)
            read.query_name, read.query_sequence = name, sequence
            read.flag, read.reference_id, read.reference_start = flag, 0, start
            read.cigarstring, read.mapping_quality = cigar, 60
            if mate >= 0:
                read.next_reference_id, read.next_reference_start = 0, mate
            stream.write(read)

        for base in ("A", "G"):
            write(f"spliced_{base}", 90, base * 40, "20M90N20M")
        positions = [101, 151, 201, 251, 301, 401]
        for pos, (ref, alt) in zip(
            positions, [(16, 4), (5, 15), (12, 8), (8, 12), (18, 2), (7, 13)], strict=True
        ):
            if pos in (151, 301):
                write(
                    "pair",
                    pos - 1,
                    "A" if pos == 151 else "G",
                    flag=65 if pos == 151 else 129,
                    mate=300 if pos == 151 else 150,
                )
            for base, count in (("A", ref), ("G", alt)):
                for index in range(count):
                    write(f"{pos}_{base}_{index}", pos - 1, base)
    pysam.index(str(bam))
    # These are the SNP/exon intersection rows: positions 151 and 201 overlap two exons.
    frame = pl.DataFrame(
        {
            "chrom": ["chr1"] * 8,
            "pos": [101, 151, 151, 201, 201, 251, 301, 401],
            "ref": ["A"] * 8,
            "alt": ["G"] * 8,
            "GT": ["0|1", "1|0", "1|0", "0|1", "0|1", "1|0", "0|1", "1|0"],
            "exon_id": ["e1", "e1", "e2", "e1", "e2", "e2", "e3", "e4"],
            "gene_id": ["g1"] * 7 + ["g2"],
        }
    ).with_columns(pl.col("chrom").cast(pl.Categorical))
    return bam, frame


@pytest.mark.rust
@pytest.mark.parametrize("with_gt", [False, True])
@pytest.mark.parametrize("order", ["forward", "reverse", "gap_first"])
@pytest.mark.parametrize("threads", [1, 4])
def test_rna_counting_matches_archive(rna_counts, monkeypatch, with_gt, order, threads):
    from tests.test_rust_python_counting_parity import python_count_snp_alleles

    bam, frame = rna_counts
    if not with_gt:
        frame = frame.drop("GT")
    if order == "reverse":
        frame = frame.reverse()
    elif order == "gap_first":
        frame = pl.concat([frame.filter(pl.col("pos") == 151), frame.filter(pl.col("pos") != 151)])
    variants = list(frame.select("pos", "ref", "alt").unique(maintain_order=True).iter_rows())
    expected = {
        row[1]: tuple(row[2:]) for row in python_count_snp_alleles(str(bam), "chr1", variants)
    }
    monkeypatch.setenv("WASP2_RUST_THREADS", str(threads))
    actual = make_count_df(str(bam), frame)
    assert actual.height == frame.height
    for row in actual.iter_rows(named=True):
        assert (row["ref_count"], row["alt_count"], row["other_count"]) == expected[row["pos"]]


@pytest.mark.rust
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("pc", [0, 1])
@pytest.mark.parametrize("model", ["single", "linear"])
@pytest.mark.parametrize("grouping", ["exon_id", "gene_id"])
def test_rna_grouping_cli_matches_archive(tmp_path, rna_counts, phased, pc, model, grouping):
    bam, frame = rna_counts
    measured = make_count_df(str(bam), frame)
    counts = tmp_path / "counts.tsv"
    measured.write_csv(counts, separator="\t")
    modules = load_archive()
    if model == "linear":
        # Explicit reference-only repairs: array/scalar shapes and the archive's
        # ignored phased argument. Likelihoods and optimizers stay unchanged.
        adapt_linear_reference(modules)
        if phased:
            original = modules["as_analysis"].parse_opt

            def parse_phased(df, disp=None, phased=False):
                return original(df, disp=disp, phased=True)

            modules["as_analysis"].parse_opt = parse_phased
    expected = (
        modules["as_analysis"]
        .get_imbalance(
            pd.DataFrame(measured.to_dict(as_series=False)),
            region_col="exon_id",
            groupby="gene_id" if grouping == "gene_id" else None,
            method=model,
            phased=phased,
            pseudocount=pc,
        )
        .rename(columns={grouping: "region"})
    )
    aliases = (
        [None]
        if grouping == "exon_id"
        else [None, "--groupby", "--group", "--parent_col", "--parent"]
    )
    for alias in aliases:
        output = tmp_path / "results.tsv"
        options = (
            ["--region_col", grouping]
            if alias is None
            else ["--region_col", "exon_id", alias, "gene_id"]
        )
        result = CliRunner().invoke(
            analysis_app,
            [
                "find-imbalance",
                str(counts),
                *options,
                "--model",
                model,
                "--pseudocount",
                str(pc),
                *(["--phased"] if phased else []),
                "-o",
                str(output),
            ],
        )
        assert result.exit_code == 0, result.exception
        comparison = compare_frames(expected, pd.read_csv(output, sep="\t"))
        assert comparison["status"] == "PASS", comparison


@pytest.mark.rust
@pytest.mark.parametrize("suffix", ["gtf", "gtf.gz", "gff3"])
def test_rna_annotation_counting_cli(tmp_path, rna_counts, suffix):
    import shutil

    from counting.__main__ import app as count_app

    if shutil.which("bedtools") is None:
        pytest.skip("bedtools is required for annotation intersection")
    bam, frame = rna_counts
    features = [
        ("gene", 81, 330, "g1", "g1"),
        ("gene", 381, 430, "g2", "g2"),
        ("exon", 91, 210, "e1", "g1"),
        ("exon", 141, 260, "e2", "g1"),
        ("exon", 291, 320, "e3", "g1"),
        ("exon", 391, 420, "e4", "g2"),
    ]
    annotation = _annotation(
        tmp_path / f"annotation.{suffix}",
        [
            (
                kind,
                start,
                end,
                f'gene_id "{gene}"; exon_id "{identifier}";'
                if suffix.startswith("gtf")
                else f"ID={identifier};Parent={gene}",
            )
            for kind, start, end, identifier, gene in features
        ],
    )
    vcf = tmp_path / "variants.vcf"
    records = frame.select("pos", "GT").unique(maintain_order=True)
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        + "".join(
            f"chr1\t{pos}\t.\tA\tG\t60\tPASS\t.\tGT\t{gt}\n" for pos, gt in records.iter_rows()
        )
    )
    output = tmp_path / "counts.tsv"
    result = CliRunner().invoke(
        count_app,
        [
            "count-variants",
            str(bam),
            str(vcf),
            "--regions",
            str(annotation),
            "--sample",
            "sample1",
            "-o",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.exception
    actual = pl.read_csv(output, separator="\t")
    if suffix == "gff3":
        actual = actual.rename({"ID": "exon_id", "Parent": "gene_id"})
    expected = make_count_df(str(bam), frame)
    columns = ["pos", "exon_id", "gene_id", "ref_count", "alt_count", "other_count"]
    assert (
        actual.select(columns).sort(columns).to_dicts()
        == expected.select(columns).sort(columns).to_dicts()
    )
