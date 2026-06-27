import pytest

from sexxy.metadata import filter_children_to_vcf, load_children_by_sex
from sexxy.vcf import compute_genotype_counts, read_vcf_samples


@pytest.fixture
def metadata_with_extra_child(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "spid\tsfid\tfather\tmother\tsex\n"
        "c1\tf1\tp1\tp2\tmale\n"
        "c2\tf1\tp1\tp2\tfemale\n"
        "c_missing\tf2\tp3\tp4\tfemale\n"
    )
    return path


@pytest.fixture
def vcf_without_c_missing(tmp_path):
    path = tmp_path / "test.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\t0/1:25\n"
    )
    return path


def test_filter_children_to_vcf_excludes_missing(
    metadata_with_extra_child, vcf_without_c_missing
):
    male, female, _ = load_children_by_sex(metadata_with_extra_child, sep="\t")
    samples = read_vcf_samples(vcf_without_c_missing)
    cohort = filter_children_to_vcf(samples, male, female)

    assert cohort.male_children == ["c1"]
    assert cohort.female_children == ["c2"]
    assert cohort.excluded_male == ()
    assert cohort.excluded_female == ("c_missing",)


def test_filter_children_to_vcf_strict_raises(
    metadata_with_extra_child, vcf_without_c_missing
):
    male, female, _ = load_children_by_sex(metadata_with_extra_child, sep="\t")
    samples = read_vcf_samples(vcf_without_c_missing)

    with pytest.raises(ValueError, match="not found in VCF header"):
        filter_children_to_vcf(samples, male, female, strict=True)


def test_compute_genotype_counts_excludes_missing_from_vcf(
    metadata_with_extra_child, vcf_without_c_missing
):
    male, female, _ = load_children_by_sex(metadata_with_extra_child, sep="\t")
    excluded: list[str] = []

    def on_excluded(excluded_male, excluded_female):
        excluded.extend(excluded_male)
        excluded.extend(excluded_female)

    result = compute_genotype_counts(
        vcf_without_c_missing,
        male,
        female,
        chromosome="chr1",
        on_excluded=on_excluded,
    )

    assert result.male_cohort_size == 1
    assert result.female_cohort_size == 1
    assert excluded == ["c_missing"]
    assert result.male_counts() == {"0/0": 1}
    assert result.female_counts() == {"0/1": 1}


def test_compute_genotype_counts_strict_raises_on_missing(
    metadata_with_extra_child, vcf_without_c_missing
):
    male, female, _ = load_children_by_sex(metadata_with_extra_child, sep="\t")

    with pytest.raises(ValueError, match="not found in VCF header"):
        compute_genotype_counts(
            vcf_without_c_missing,
            male,
            female,
            chromosome="chr1",
            strict=True,
        )


def test_compute_genotype_counts_raises_when_no_children_in_vcf(tmp_path):
    meta = tmp_path / "metadata.tsv"
    meta.write_text(
        "spid\tfather\tmother\tsex\n"
        "ghost\tp1\tp2\tfemale\n"
    )
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT:DP\t0/0:30\n"
    )
    _, female, _ = load_children_by_sex(meta, sep="\t")

    with pytest.raises(ValueError, match="no male or female children"):
        compute_genotype_counts(vcf, [], female, chromosome="chr1")
