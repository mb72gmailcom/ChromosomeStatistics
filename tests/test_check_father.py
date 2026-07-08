import json
from pathlib import Path

import pytest

from sexxy.cli import main
from sexxy.metadata import load_children_by_sex
from sexxy.vcf import compute_genotype_counts, is_chry


@pytest.fixture
def chry_metadata(tmp_path: Path) -> Path:
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "spid\tfather\tmother\tsex\n"
        "c1\tp1\tm1\tmale\n"
        "c2\tp1\tm1\tfemale\n"
    )
    return path


@pytest.fixture
def chry_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "chrY.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\tp1\n"
        "chrY\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/0:30\t0/0:30\n"
        "chrY\t200\trs2\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t0/0:30\t1/1:30\n"
        "chrY\t300\trs3\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t0/1:30\t1/1:30\n"
        "chrY\t400\trs4\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t0/0:30\t0/1:30\n"
    )
    return path


def _male_father_map(metadata: Path, male_children: list[str]) -> dict[str, str]:
    _, _, children_rows = load_children_by_sex(metadata, sep="\t")
    return {
        row["spid"]: row["father"]
        for row in children_rows
        if row["spid"] in male_children
    }


def test_is_chry():
    assert is_chry("chrY")
    assert is_chry("Y")
    assert not is_chry("chrX")


def test_check_father_chrY_counts(chry_vcf: Path, chry_metadata: Path):
    male, female, _ = load_children_by_sex(chry_metadata, sep="\t")
    result = compute_genotype_counts(
        chry_vcf,
        male,
        female,
        chromosome="chrY",
        check_father=True,
        male_father_by_child=_male_father_map(chry_metadata, male),
    )
    assert result.male_counts() == {
        "0/0": 1,
        "0/1": 2,
        "1/1": 1,
        "0/0_f": 1,
        "0/1_f": 1,
        "1/1_f": 1,
    }
    assert result.female_counts() == {"0/0": 2, "0/1": 1}


def test_check_father_off(chry_vcf: Path, chry_metadata: Path):
    male, female, _ = load_children_by_sex(chry_metadata, sep="\t")
    result = compute_genotype_counts(
        chry_vcf,
        male,
        female,
        chromosome="chrY",
    )
    assert result.male_counts() == {"0/0": 1, "0/1": 2, "1/1": 1}
    assert all(not key.endswith("_f") for key in result.male_counts())


def test_check_father_ignored_on_autosome(tmp_path: Path, chry_metadata: Path):
    vcf = tmp_path / "chr1.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tp1\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/0:30\n"
    )
    male, female, _ = load_children_by_sex(chry_metadata, sep="\t")
    result = compute_genotype_counts(
        vcf,
        male,
        female,
        chromosome="chr1",
        check_father=True,
        male_father_by_child=_male_father_map(chry_metadata, male),
    )
    assert result.male_counts() == {"0/0": 1}
    assert all(not key.endswith("_f") for key in result.male_counts())


def test_check_father_missing_father_in_vcf(tmp_path: Path, chry_metadata: Path):
    vcf = tmp_path / "chrY.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n"
        "chrY\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\n"
    )
    male, female, _ = load_children_by_sex(chry_metadata, sep="\t")
    result = compute_genotype_counts(
        vcf,
        male,
        female,
        chromosome="chrY",
        check_father=True,
        male_father_by_child=_male_father_map(chry_metadata, male),
    )
    assert result.male_counts() == {"0/0": 1}
    assert "0/0_f" not in result.male_counts()


def test_check_father_respects_child_filters(chry_vcf: Path, chry_metadata: Path):
    male, female, _ = load_children_by_sex(chry_metadata, sep="\t")
    result = compute_genotype_counts(
        chry_vcf,
        male,
        female,
        chromosome="chrY",
        min_dp=35,
        check_father=True,
        male_father_by_child=_male_father_map(chry_metadata, male),
    )
    assert result.male_counts() == {}


def test_cli_check_father_chrY(chry_vcf: Path, chry_metadata: Path, tmp_path: Path):
    out_dir = tmp_path / "results"
    rc = main(
        [
            str(chry_vcf),
            str(chry_metadata),
            "--chromosome",
            "chrY",
            "--check-father",
            "--output-dir",
            str(out_dir),
        ]
    )
    assert rc == 0
    male_data = json.loads((out_dir / "counts.chrY.male.json").read_text())
    assert male_data["gt_counts"]["0/0_f"] == 1
    assert male_data["gt_counts"]["0/1_f"] == 1
    assert male_data["gt_counts"]["1/1_f"] == 1
    params = json.loads((out_dir / "counts.chrY.params.json").read_text())
    assert params["filters"]["check_father"] is True


def test_cli_check_father_warns_on_non_chrY(chry_vcf: Path, chry_metadata: Path, tmp_path: Path, capsys):
    out_dir = tmp_path / "results"
    rc = main(
        [
            str(chry_vcf),
            str(chry_metadata),
            "--chromosome",
            "chr1",
            "--check-father",
            "--output-dir",
            str(out_dir),
        ]
    )
    assert rc == 0
    stderr = capsys.readouterr().err
    assert "--check-father applies only to chrY" in stderr
