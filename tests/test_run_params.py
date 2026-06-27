import json
from pathlib import Path

import pytest

from sexxy.cli import main
from sexxy.results import write_run_params


@pytest.fixture
def run_fixtures(tmp_path: Path):
    meta = tmp_path / "metadata.tsv"
    meta.write_text(
        "spid\tfather\tmother\tsex\n"
        "c1\tp1\tp2\tmale\n"
        "c2\tp1\tp2\tfemale\n"
    )
    vcf = tmp_path / "test.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
    )
    out_dir = tmp_path / "results"
    return vcf, meta, out_dir


def test_write_run_params(run_fixtures, tmp_path: Path):
    vcf, meta, out_dir = run_fixtures
    params = {
        "version": "0.1.0",
        "chromosome": "chr1",
        "inputs": {"vcf": str(vcf), "metadata": str(meta)},
        "cohort": {"male_children": 1, "female_children": 1},
        "output_files": [],
    }
    path = write_run_params(out_dir / "counts.chr1", "chr1", params)
    assert path.name == "counts.chr1.params.json"
    assert path.is_file()
    data = json.loads(path.read_text())
    assert data["chromosome"] == "chr1"
    assert data["params_file"] == str(path)
    assert str(path) in data["output_files"]


def test_cli_writes_params_file(run_fixtures, tmp_path: Path, capsys):
    vcf, meta, out_dir = run_fixtures
    rc = main(
        [
            str(vcf),
            str(meta),
            "--chromosome",
            "chr1",
            "--output-dir",
            str(out_dir),
        ]
    )
    assert rc == 0
    params_path = out_dir / "counts.chr1.params.json"
    assert params_path.is_file()
    data = json.loads(params_path.read_text())
    assert data["inputs"]["vcf"] == str(vcf)
    assert data["cohort"]["male_children"] == 1
    assert data["cohort"]["female_children"] == 1
    assert (out_dir / "counts.chr1.male.json").is_file()
    assert (out_dir / "counts.chr1.female.json").is_file()
    assert str(params_path) in data["output_files"]
    stderr = capsys.readouterr().err
    assert "counts.chr1.params.json" in stderr
