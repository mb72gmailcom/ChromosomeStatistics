from pathlib import Path

import pytest

from sexxy.intervals import ExcludeIntervals, load_exclude_intervals
from sexxy.metadata import load_children_by_sex
from sexxy.vcf import compute_genotype_counts


@pytest.fixture
def repeat_intervals_file(tmp_path: Path) -> Path:
    path = tmp_path / "chrY.repeats.tsv"
    path.write_text(
        "chrY\t0\t2781790\n"
        "chrY\t2782979\t2783283\n"
        "chrY\t2783815\t2784149\n"
    )
    return path


def test_load_exclude_intervals(repeat_intervals_file: Path):
    intervals = load_exclude_intervals(repeat_intervals_file, "chrY")
    assert len(intervals) == 3
    assert intervals.contains(0)
    assert intervals.contains(2781789)
    assert not intervals.contains(2781790)
    assert intervals.contains(2782979)
    assert not intervals.contains(2783283)
    assert not intervals.contains(3000000)


def test_load_exclude_intervals_chrom_alias(repeat_intervals_file: Path):
    intervals = load_exclude_intervals(repeat_intervals_file, "Y")
    assert len(intervals) == 3


def test_load_exclude_intervals_requires_sorted(tmp_path: Path):
    path = tmp_path / "bad.tsv"
    path.write_text("chr1\t200\t300\nchr1\t100\t150\n")
    with pytest.raises(ValueError, match="sorted by start"):
        load_exclude_intervals(path, "chr1")


def test_load_exclude_intervals_no_matching_chrom(tmp_path: Path):
    path = tmp_path / "empty.tsv"
    path.write_text("chr2\t100\t200\n")
    with pytest.raises(ValueError, match="no intervals"):
        load_exclude_intervals(path, "chr1")


def test_compute_genotype_counts_exclude_repeats(tmp_path: Path, repeat_intervals_file: Path):
    meta = tmp_path / "metadata.tsv"
    meta.write_text("spid\tfather\tmother\tsex\nc1\tp1\tp2\tmale\n")
    vcf = tmp_path / "chrY.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n"
        "chrY\t100\trs_in\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\n"
        "chrY\t2781790\trs_out\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\n"
        "chrY\t2783000\trs_in2\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\n"
    )
    male, female, _ = load_children_by_sex(meta, sep="\t")
    result = compute_genotype_counts(
        vcf,
        male,
        female,
        chromosome="chrY",
        exclude_repeats=repeat_intervals_file,
    )
    assert result.excluded_repeat_rows == 2
    assert result.exclude_repeat_intervals == 3
    assert result.male_counts() == {"0/1": 1}


def test_exclude_intervals_preloaded():
    intervals = ExcludeIntervals([(100, 200)])
    assert intervals.contains(100)
    assert intervals.contains(199)
    assert not intervals.contains(200)
