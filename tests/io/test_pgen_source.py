"""Exercise the installed optional PGEN reader through the public BED interface."""

from pathlib import Path

import numpy as np
import pytest

from wasp2.io.compat import variants_to_bed

pgenlib = pytest.importorskip("pgenlib")


@pytest.fixture
def pgen(tmp_path: Path):
    path = tmp_path / "variants.pgen"
    path.with_suffix(".psam").write_text("#IID\ndonor.1\ndonor1\n")
    path.with_suffix(".pvar").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n1\t100\trs1\tA\tG\n1\t200\trs2\tC\tT\n1\t300\trs3\tA\tATT\n"
    )
    with pgenlib.PgenWriter(bytes(path), 2, variant_ct=3, hardcall_phase_present=True) as writer:
        writer.append_alleles(np.array([1, 0, 0, 0], dtype=np.int32), all_phased=True)
        writer.append_alleles(np.array([0, 1, 1, 1], dtype=np.int32), all_phased=False)
        writer.append_alleles(np.array([0, 1, 0, 0], dtype=np.int32), all_phased=True)
    return path


@pytest.mark.parametrize("samples", [["donor.1"], ["donor1", "donor.1"]])
def test_pgen_bed_keeps_phase_orientation_and_sample_columns(pgen, tmp_path, samples):
    output = variants_to_bed(pgen, tmp_path / "variants.bed", samples=samples)
    rows = [line.split("\t") for line in output.read_text().splitlines()]
    assert len(rows) == 2
    assert all(len(row) == 5 + len(samples) for row in rows)
    assert rows[0][5 + samples.index("donor.1")] == "1|0"
    assert rows[1][5 + samples.index("donor.1")] == "0/1"


def test_pgen_export_indel_flag(pgen, tmp_path):
    output = variants_to_bed(
        pgen, tmp_path / "variants.bed", samples=["donor.1"], include_indels=True
    )
    assert len(output.read_text().splitlines()) == 3
