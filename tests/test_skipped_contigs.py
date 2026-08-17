import json
from pathlib import Path

from sexxy.cli import main
from sexxy.gnomad import GnomadAfStore
from sexxy.metadata import load_children_by_sex
from sexxy.results import write_skipped_contigs
from sexxy.vcf import _chrom_key, compute_genotype_counts


def test_compute_genotype_counts_skipped_contigs(tmp_path: Path):
    meta = tmp_path / "metadata.tsv"
    meta.write_text(
        "spid\tfather\tmother\tsex\n"
        "c1\tp1\tp2\tmale\n"
        "c2\tp1\tp2\tfemale\n"
    )
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
        "chr1_KI270706v1_random\t200\trs2\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/0:30\n"
        "chr1_KI270706v1_random\t300\trs3\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t0/0:30\n"
        "chr2\t400\trs4\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/0:30\n"
    )
    male, female, _ = load_children_by_sex(meta, sep="\t")
    result = compute_genotype_counts(vcf, male, female, chromosome="chr1")

    assert result.skipped_contigs == {
        "chr1_KI270706v1_random": 2,
        "chr2": 1,
    }
    assert result.male_counts() == {"0/0": 1}
    assert result.female_counts() == {"0/1": 1}


def test_skipped_contigs_not_treated_as_gnomad_hits(tmp_path: Path):
    meta = tmp_path / "metadata.tsv"
    meta.write_text(
        "spid\tfather\tmother\tsex\n"
        "c1\tp1\tp2\tmale\n"
        "c2\tp1\tp2\tfemale\n"
    )
    vcf = tmp_path / "mixed.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
        "chr1_KI270706v1_random\t100\trs2\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t1/1:30\n"
        "chr1\t300\trs3\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t0/0:30\n"
    )
    af_dir = tmp_path / "gnomad.v4" / "chr1"
    af_dir.mkdir(parents=True)
    (af_dir / "chr1-common-af.json").write_text(json.dumps({"chr1_100_A_G": 0.05}))
    male, female, _ = load_children_by_sex(meta, sep="\t")
    result = compute_genotype_counts(
        vcf,
        male,
        female,
        chromosome="chr1",
        gnomad_af=GnomadAfStore(tmp_path / "gnomad.v4"),
        common_freq_cutoff=0.01,
    )
    assert result.skipped_contigs == {"chr1_KI270706v1_random": 1}
    assert result.male_counts() == {"1/1": 1}
    assert result.female_counts() == {"0/0": 1}
    assert result.skipped_common_variants == 1


def test_write_skipped_contigs(tmp_path: Path):
    prefix = tmp_path / "counts.chr1"
    path = write_skipped_contigs(
        prefix,
        "chr1",
        vcf="/data/chr1.vcf.gz",
        target_chrom_key=_chrom_key("chr1"),
        skipped_contigs={"chr1_KI270706v1_random": 2},
    )
    assert path.name == "counts.chr1.skipped_contigs.json"
    data = json.loads(path.read_text())
    assert data["chromosome"] == "chr1"
    assert data["target_chrom_key"] == "1"
    assert data["total_skipped_rows"] == 2
    assert data["skipped_contigs"] == {"chr1_KI270706v1_random": 2}


def test_cli_writes_skipped_contigs_file(tmp_path: Path, capsys):
    meta = tmp_path / "metadata.tsv"
    meta.write_text("spid\tfather\tmother\tsex\nc1\tp1\tp2\tmale\n")
    vcf = tmp_path / "test.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\n"
        "chr1_GL000192v1_random\t200\trs2\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\n"
    )
    out_dir = tmp_path / "results"
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
    skipped_path = out_dir / "counts.chr1.skipped_contigs.json"
    assert skipped_path.is_file()
    data = json.loads(skipped_path.read_text())
    assert data["total_skipped_rows"] == 1
    assert data["skipped_contigs"] == {"chr1_GL000192v1_random": 1}
    params = json.loads((out_dir / "counts.chr1.params.json").read_text())
    assert params["total_skipped_rows"] == 1
    assert str(skipped_path) in params["output_files"]
    assert "non-matching contig" in capsys.readouterr().err
