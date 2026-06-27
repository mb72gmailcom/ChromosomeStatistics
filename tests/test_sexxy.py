import json
from pathlib import Path

import pytest

from sexxy.chrx import chrx_region, CHRX_REGION_ORDER
from sexxy.gnomad import GnomadAfStore, gnomad_af_path, load_gnomad_af_json
from sexxy.metadata import filter_children_to_vcf, load_children_by_sex
from sexxy.results import (
    CHRX_FEMALE_OUTPUT_KEYS,
    CHRX_MALE_OUTPUT_KEYS,
    FEMALE_OUTPUT_KEYS,
    MALE_OUTPUT_KEYS,
    resolve_output_target,
    write_genotype_count_results,
)
from sexxy.vcf import (
    _allele_balance,
    _allele_balance_from_parts,
    _field_indices,
    _passes_genotype_filters,
    _parse_sample_fields,
    chrom_matches,
    compute_genotype_counts,
    get_n_fields,
    is_snv,
    read_vcf_samples,
)


@pytest.fixture
def metadata_path(tmp_path: Path) -> Path:
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "spid\tsfid\tfather\tmother\tsex\n"
        "c1\tf1\tp1\tp2\tmale\n"
        "c2\tf1\tp1\tp2\tfemale\n"
        "c3\tf2\t\tp3\tfemale\n"
        "p1\tf1\t\t\tmale\n"
    )
    return path


@pytest.fixture
def vcf_path(tmp_path: Path) -> Path:
    path = tmp_path / "test.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
        "chr1\t200\t.\tA\tAT\t.\t.\t.\tGT:DP\t0/0:30\t0/0:30\n"
        "chr1\t300\trs3\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t0/0:30\n"
    )
    return path


def test_load_children_by_sex(metadata_path: Path):
    male, female, children = load_children_by_sex(metadata_path, sep="\t")
    assert male == ["c1"]
    assert female == ["c2"]
    assert len(children) == 2


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Male", "male"),
        ("Female", "female"),
        ("male", "male"),
        ("female", "female"),
        ("M", "male"),
        ("F", "female"),
        ("1", "male"),
        ("2", "female"),
        (1, "male"),
        (2, "female"),
        (1.0, "male"),
        (2.0, "female"),
        ("1.0", "male"),
        ("2.0", "female"),
    ],
)
def test_normalize_sex(raw, expected):
    assert _normalize_sex(raw) == expected


def test_load_children_by_sex_male_female_labels(tmp_path: Path):
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "spid\tsfid\tfather\tmother\tsex\n"
        "c1\tf1\tp1\tp2\tMale\n"
        "c2\tf1\tp1\tp2\tFemale\n"
    )
    male, female, _ = load_children_by_sex(path, sep="\t")
    assert male == ["c1"]
    assert female == ["c2"]


def test_load_children_by_sex_numeric_labels(tmp_path: Path):
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "spid\tsfid\tfather\tmother\tsex\n"
        "c1\tf1\tp1\tp2\t1\n"
        "c2\tf1\tp1\tp2\t2\n"
    )
    male, female, _ = load_children_by_sex(path, sep="\t")
    assert male == ["c1"]
    assert female == ["c2"]


def test_get_n_fields():
    line = "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/0\n"
    assert get_n_fields(line, 5) == ["chr1", "100", "rs1", "A", "G"]


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("A", "G", True),
        ("AT", "A", False),
        ("A", "AT", False),
        ("A", "G,T", False),
    ],
)
def test_is_snv(ref, alt, expected):
    assert is_snv(ref, alt) is expected


def test_compute_genotype_counts_skips_non_snvs(vcf_path: Path, metadata_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(vcf_path, male, female, chromosome="chr1")
    # vcf_path includes an indel at chr1:200 (ref=A, alt=AT) which must be skipped
    assert result.male_counts() == {"0/0": 1, "1/1": 1}
    assert result.female_counts() == {"0/1": 1, "0/0": 1}


def test_read_vcf_samples(vcf_path: Path):
    assert read_vcf_samples(vcf_path) == ["c1", "c2"]


def test_chromosome_filter_and_prefix(tmp_path: Path, metadata_path: Path):
    path = tmp_path / "multi_chrom.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
        "chr2\t100\trs9\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t1/1:30\n"
    )
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(path, male, female, chromosome="1")
    assert chrom_matches("chr1", "1")
    assert result.male_counts() == {"0/0": 1}
    assert result.female_counts() == {"0/1": 1}


def test_compute_genotype_counts_with_af_filter(vcf_path: Path, metadata_path: Path, tmp_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    af = {"rs1": 0.05, "rs3": 0.001}
    result = compute_genotype_counts(
        vcf_path,
        male,
        female,
        chromosome="chr1",
        allele_freqs=af,
        common_freq_cutoff=0.01,
    )
    # rs1 filtered (AF 0.05), only rs3 remains
    assert result.male_counts() == {"1/1": 1}
    assert result.female_counts() == {"0/0": 1}


@pytest.fixture
def gnomad_dir(tmp_path: Path) -> Path:
    chrm = "chr1"
    af_dir = tmp_path / "gnomad.v4"
    chrom_dir = af_dir / chrm
    chrom_dir.mkdir(parents=True)
    af_file = chrom_dir / f"{chrm}-common-af.json"
    af_file.write_text(json.dumps({"rs1": 0.05, "rs3": 0.001}))
    return af_dir


def test_gnomad_af_path(gnomad_dir: Path):
    assert gnomad_af_path(gnomad_dir, "chr1") == gnomad_dir / "chr1" / "chr1-common-af.json"


def test_load_gnomad_af_json(gnomad_dir: Path):
    path = gnomad_af_path(gnomad_dir, "chr1")
    assert load_gnomad_af_json(path) == {"rs1": 0.05, "rs3": 0.001}


def test_compute_genotype_counts_with_gnomad_af(
    vcf_path: Path, metadata_path: Path, gnomad_dir: Path
):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    store = GnomadAfStore(gnomad_dir)
    result = compute_genotype_counts(
        vcf_path, male, female, chromosome="chr1", gnomad_af=store, common_freq_cutoff=0.01
    )
    assert result.male_counts() == {"1/1": 1}
    assert result.female_counts() == {"0/0": 1}


@pytest.fixture
def filtered_vcf_path(tmp_path: Path) -> Path:
    fmt = "GT:DP:AD:SB:GQ:PL"
    path = tmp_path / "filtered.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        f"chr1\t100\trs1\tA\tG\t.\t.\t.\t{fmt}\t"
        "0/0:30:20,10:.:99:.\t0/1:25:12,13:.:50:.\n"
        f"chr1\t200\trs2\tA\tG\t.\t.\t.\t{fmt}\t"
        "1/1:30:2,28:.:99:.\t0/1:30:5,25:.:99:.\n"
        f"chr1\t300\trs3\tA\tG\t.\t.\t.\t{fmt}\t"
        "0/1:30:15,15:.:99:.\t0/0:30:30,0:.:99:.\n"
    )
    return path


def test_min_dp_filter(filtered_vcf_path: Path, metadata_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        filtered_vcf_path, male, female, chromosome="chr1", min_dp=30
    )
    # c1 0/1 at rs1 has DP 25 and is dropped; c2 0/1 at rs1 kept
    assert result.male_counts() == {"0/0": 1, "1/1": 1, "0/1": 1}
    assert result.female_counts() == {"0/1": 1, "0/0": 1}


def test_min_gq_filter(filtered_vcf_path: Path, metadata_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        filtered_vcf_path, male, female, chromosome="chr1", min_gq=60
    )
    # c2 0/1 at rs1 has GQ 50 and is dropped
    assert result.male_counts() == {"0/0": 1, "1/1": 1, "0/1": 1}
    assert result.female_counts() == {"0/1": 1, "0/0": 1}


def test_ab_threshold_filter(filtered_vcf_path: Path, metadata_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        filtered_vcf_path, male, female, chromosome="chr1", ab_threshold=0.2
    )
    # all 0/1 and 1/1 calls have AB > 0.2
    assert result.male_counts() == {"0/0": 1, "1/1": 1, "0/1": 1}
    assert result.female_counts() == {"0/1": 2, "0/0": 1}


def test_allele_balance_from_ad():
    fields = _parse_sample_fields("GT:DP:AD", "0/1:30:12,18")
    assert _allele_balance(fields) == pytest.approx(0.6)
    indices = _field_indices("GT:DP:AD")
    assert _allele_balance_from_parts("0/1:30:12,18".split(":"), indices) == pytest.approx(0.6)


def test_allele_balance_haploid_ad():
    fields = _parse_sample_fields("GT:DP:AD", "1:30:25")
    assert _allele_balance(fields) == 1.0
    indices = _field_indices("GT:DP:AD")
    assert _allele_balance_from_parts("1:30:25".split(":"), indices) == 1.0


def test_passes_genotype_filters_ab_only_on_het_hom_alt():
    fields = _parse_sample_fields("GT:AD", "0/0:30,0")
    assert _passes_genotype_filters(fields, min_gq=None, min_dp=None, ab_threshold=0.2)

    het_bad = _parse_sample_fields("GT:AD", "0/1:25,5")
    assert not _passes_genotype_filters(
        het_bad, min_gq=None, min_dp=None, ab_threshold=0.2
    )


@pytest.mark.parametrize(
    ("pos", "region"),
    [
        (10_001, "Par1"),
        (2_781_479, "Par1"),
        (2_781_489, "noPar"),
        (155_701_382, "noPar"),
        (155_701_383, "Par2"),
        (156_030_895, "Par2"),
        (2_781_480, None),
    ],
)
def test_chrx_region(pos, region):
    assert chrx_region(pos) == region


@pytest.fixture
def chrx_vcf_path(tmp_path: Path) -> Path:
    path = tmp_path / "chrx.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chrX\t10001\tp1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
        "chrX\t2781480\tgap\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t1/1:30\n"
        "chrX\t2781489\tnp\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t0/0:30\n"
        "chrX\t155701383\tp2\tA\tG\t.\t.\t.\tGT:DP\t1/1:30\t0/1:30\n"
    )
    return path


def test_compute_genotype_counts_chrx_regions(chrx_vcf_path: Path, metadata_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        chrx_vcf_path, male, female, chromosome="chrX"
    )
    assert result.regions == ("Par1", "noPar", "Par2")
    assert result.male_counts("Par1") == {"0/0": 1}
    assert result.female_counts("Par1") == {"0/1": 1}
    assert result.male_counts("noPar") == {"0/1": 1}
    assert result.female_counts("noPar") == {"0/0": 1}
    assert result.male_counts("Par2") == {"1/1": 1}
    assert result.female_counts("Par2") == {"0/1": 1}


def test_write_chrx_region_output_files(chrx_vcf_path: Path, metadata_path: Path, tmp_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        chrx_vcf_path, male, female, chromosome="chrX"
    )
    prefix = tmp_path / "counts.chrX"
    paths = write_genotype_count_results(
        result, prefix, male_children=result.male_cohort_size, female_children=result.female_cohort_size
    )
    assert len(paths) == 6
    names = {p.name for p in paths}
    assert names == {
        "counts.chrX.male.Par1.json",
        "counts.chrX.male.noPar.json",
        "counts.chrX.male.Par2.json",
        "counts.chrX.female.Par1.json",
        "counts.chrX.female.noPar.json",
        "counts.chrX.female.Par2.json",
    }
    male_par1 = json.loads((tmp_path / "counts.chrX.male.Par1.json").read_text())
    assert male_par1 == {
        "chromosome": "chrX",
        "gt_counts": {"0/0": 1},
        "male_children": 1,
        "region": "Par1",
        "sex": "male",
    }
    for region in CHRX_REGION_ORDER:
        for sex in ("male", "female"):
            data = json.loads((tmp_path / f"counts.chrX.{sex}.{region}.json").read_text())
            expected_keys = CHRX_MALE_OUTPUT_KEYS if sex == "male" else CHRX_FEMALE_OUTPUT_KEYS
            assert set(data.keys()) == set(expected_keys)
            assert data["region"] == region
            assert data["sex"] == sex
            assert "female_children" not in data if sex == "male" else "male_children" not in data


def test_chrx_nonpar_filter_overrides(tmp_path: Path, metadata_path: Path):
    fmt = "GT:DP:AD:SB:GQ:PL"
    path = tmp_path / "chrx_filters.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        f"chrX\t10001\tp1\tA\tG\t.\t.\t.\t{fmt}\t"
        "0/1:25:12,13:.:50:.\t0/0:30:20,10:.:99:.\n"
        f"chrX\t2781489\tnp\tA\tG\t.\t.\t.\t{fmt}\t"
        "0/1:25:12,13:.:50:.\t0/0:30:20,10:.:99:.\n"
    )
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")

    strict = compute_genotype_counts(
        path, male, female, chromosome="chrX", min_dp=30
    )
    assert strict.male_counts("Par1") == {}
    assert strict.male_counts("noPar") == {}
    assert strict.female_counts("Par1") == {"0/0": 1}
    assert strict.female_counts("noPar") == {"0/0": 1}

    relaxed_nonpar = compute_genotype_counts(
        path, male, female, chromosome="chrX", min_dp=30, min_dp_nonpar=20
    )
    assert relaxed_nonpar.male_counts("Par1") == {}
    assert relaxed_nonpar.male_counts("noPar") == {"0/1": 1}
    assert relaxed_nonpar.female_counts("Par1") == {"0/0": 1}
    assert relaxed_nonpar.female_counts("noPar") == {"0/0": 1}


def test_chrx_nonpar_overrides_do_not_apply_to_females(tmp_path: Path, metadata_path: Path):
    fmt = "GT:DP:AD:SB:GQ:PL"
    path = tmp_path / "chrx_female_noPar.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        f"chrX\t2781489\tnp\tA\tG\t.\t.\t.\t{fmt}\t"
        "0/1:25:12,13:.:50:.\t0/1:25:12,13:.:50:.\n"
    )
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        path, male, female, chromosome="chrX", min_dp=30, min_dp_nonpar=20
    )
    assert result.male_counts("noPar") == {"0/1": 1}
    assert result.female_counts("noPar") == {}


def test_write_output_dir(chrx_vcf_path: Path, metadata_path: Path, tmp_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(
        chrx_vcf_path, male, female, chromosome="chrX"
    )
    out_dir = tmp_path / "results"
    target = resolve_output_target(None, out_dir, "chrX")
    paths = write_genotype_count_results(
        result, target, male_children=result.male_cohort_size, female_children=result.female_cohort_size
    )
    assert out_dir.is_dir()
    assert all(p.parent == out_dir for p in paths)
    assert (out_dir / "counts.chrX.male.Par1.json") in paths


def test_write_autosome_sex_output_files(vcf_path: Path, metadata_path: Path, tmp_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(vcf_path, male, female, chromosome="chr1")
    prefix = tmp_path / "counts.chr1"
    paths = write_genotype_count_results(
        result, prefix, male_children=result.male_cohort_size, female_children=result.female_cohort_size
    )
    assert len(paths) == 2
    assert {p.name for p in paths} == {"counts.chr1.male.json", "counts.chr1.female.json"}
    male_data = json.loads((tmp_path / "counts.chr1.male.json").read_text())
    assert set(male_data.keys()) == set(MALE_OUTPUT_KEYS)
    assert male_data["sex"] == "male"
    assert "region" not in male_data
    female_data = json.loads((tmp_path / "counts.chr1.female.json").read_text())
    assert set(female_data.keys()) == set(FEMALE_OUTPUT_KEYS)
    assert female_data["sex"] == "female"


def test_write_creates_nested_output_dir(tmp_path: Path, metadata_path: Path, vcf_path: Path):
    male, female, _ = load_children_by_sex(metadata_path, sep="\t")
    result = compute_genotype_counts(vcf_path, male, female, chromosome="chr1")
    out_dir = tmp_path / "nested" / "results"
    target = resolve_output_target(None, out_dir, "chr1")
    paths = write_genotype_count_results(
        result, target, male_children=result.male_cohort_size, female_children=result.female_cohort_size
    )
    assert len(paths) == 2
    assert (out_dir / "counts.chr1.male.json").is_file()
    assert (out_dir / "counts.chr1.female.json").is_file()
