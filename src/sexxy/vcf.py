"""Single-pass VCF scanning for SNV genotype counts by sex."""

from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, IO, Mapping

from sexxy.chrx import CHRX_REGION_ORDER, chrx_region, is_chrx
from sexxy.gnomad import GnomadAfStore
from sexxy.metadata import filter_children_to_vcf, sample_column_indices
from sexxy.results import GenotypeCountResult


def get_n_fields(line: str, n: int) -> list[str]:
    """Return the first *n* tab-separated fields from a VCF data line."""
    parts = line.rstrip("\n").split("\t")
    if len(parts) < n:
        raise ValueError(f"expected at least {n} fields, got {len(parts)}")
    return parts[:n]


def _open_vcf(path: str | Path) -> IO[str]:
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_vcf_samples(vcf_path: str | Path) -> list[str]:
    """Read sample IDs from the ``#CHROM`` header line."""
    with _open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError("VCF header line (#CHROM) not found")


def _parse_gt(sample_field: str) -> str:
    return sample_field.split(":", 1)[0]


def _unphased_gt(gt: str) -> str:
    return gt.replace("|", "/")


def _parse_sample_fields(format_str: str, sample_field: str) -> dict[str, str]:
    keys = format_str.split(":")
    vals = sample_field.split(":")
    if len(vals) < len(keys):
        vals = vals + ["."] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def _allele_balance(fields: Mapping[str, str]) -> float | None:
    ab = fields.get("AB", ".")
    if ab not in (".", ""):
        try:
            return float(ab)
        except ValueError:
            pass

    ad = fields.get("AD", ".")
    if ad in (".", ""):
        return None
    parts = ad.split(",")
    if len(parts) == 1:
        try:
            int(parts[0])
        except ValueError:
            return None
        return 1.0
    if len(parts) < 2:
        return None
    try:
        ref, alt = int(parts[0]), int(parts[1])
    except ValueError:
        return None
    total = ref + alt
    if total == 0:
        return None
    return alt / total


def _passes_genotype_filters(
    fields: Mapping[str, str],
    *,
    min_gq: float | None,
    min_dp: int | None,
    ab_threshold: float | None,
) -> bool:
    gt = _unphased_gt(fields.get("GT", "."))
    if gt in (".", "./."):
        return min_gq is None and min_dp is None and ab_threshold is None

    if min_gq is not None:
        gq = fields.get("GQ", ".")
        if gq in (".", ""):
            return False
        try:
            if float(gq) < min_gq:
                return False
        except ValueError:
            return False

    if min_dp is not None:
        dp = fields.get("DP", ".")
        if dp in (".", ""):
            return False
        try:
            if int(float(dp)) < min_dp:
                return False
        except ValueError:
            return False

    if ab_threshold is not None and gt in ("0/1", "1/1"):
        ab = _allele_balance(fields)
        if ab is None or ab <= ab_threshold:
            return False

    return True


def _count_if_passes(
    sample_field: str,
    format_str: str,
    counts: Counter[str],
    *,
    min_gq: float | None,
    min_dp: int | None,
    ab_threshold: float | None,
) -> None:
    fields = _parse_sample_fields(format_str, sample_field)
    if not _passes_genotype_filters(
        fields,
        min_gq=min_gq,
        min_dp=min_dp,
        ab_threshold=ab_threshold,
    ):
        return
    counts[fields.get("GT", ".").split(":")[0]] += 1


@dataclass(frozen=True)
class _FieldIndices:
    """Colon-field indices for FORMAT tags (GT is always index 0)."""

    gq: int | None
    dp: int | None
    ad: int | None
    ab: int | None


@dataclass(frozen=True)
class _FilterSettings:
    min_gq: float | None
    min_dp: int | None
    ab_threshold: float | None

    @property
    def count_mode(self) -> str:
        if self.min_gq is None and self.min_dp is None and self.ab_threshold is None:
            return "gt_only"
        if self.ab_threshold is not None:
            return "gq_dp_ab"
        return "gq_dp"


def _field_indices(format_str: str) -> _FieldIndices:
    keys = format_str.split(":")
    idx = {k: i for i, k in enumerate(keys)}
    return _FieldIndices(
        gq=idx.get("GQ"),
        dp=idx.get("DP"),
        ad=idx.get("AD"),
        ab=idx.get("AB"),
    )


def _allele_balance_from_parts(parts: list[str], indices: _FieldIndices) -> float | None:
    if indices.ab is not None and indices.ab < len(parts):
        ab = parts[indices.ab]
        if ab not in (".", ""):
            try:
                return float(ab)
            except ValueError:
                pass

    if indices.ad is None or indices.ad >= len(parts):
        return None
    ad = parts[indices.ad]
    if ad in (".", ""):
        return None
    ad_parts = ad.split(",")
    if len(ad_parts) == 1:
        try:
            int(ad_parts[0])
        except ValueError:
            return None
        return 1.0
    if len(ad_parts) < 2:
        return None
    try:
        ref, alt = int(ad_parts[0]), int(ad_parts[1])
    except ValueError:
        return None
    total = ref + alt
    if total == 0:
        return None
    return alt / total


def _passes_gq_dp(
    parts: list[str],
    indices: _FieldIndices,
    *,
    min_gq: float | None,
    min_dp: int | None,
) -> bool:
    gt = _unphased_gt(parts[0] if parts else ".")
    if gt in (".", "./."):
        return False

    if min_gq is not None:
        if indices.gq is None or indices.gq >= len(parts):
            return False
        gq = parts[indices.gq]
        if gq in (".", ""):
            return False
        try:
            if float(gq) < min_gq:
                return False
        except ValueError:
            return False

    if min_dp is not None:
        if indices.dp is None or indices.dp >= len(parts):
            return False
        dp = parts[indices.dp]
        if dp in (".", ""):
            return False
        try:
            if int(float(dp)) < min_dp:
                return False
        except ValueError:
            return False

    return True


def _passes_gq_dp_ab(
    parts: list[str],
    indices: _FieldIndices,
    settings: _FilterSettings,
) -> bool:
    gt = _unphased_gt(parts[0] if parts else ".")
    if gt in (".", "./."):
        return False

    if not _passes_gq_dp(parts, indices, min_gq=settings.min_gq, min_dp=settings.min_dp):
        return False

    if settings.ab_threshold is not None and gt in ("0/1", "1/1"):
        ab = _allele_balance_from_parts(parts, indices)
        if ab is None or ab <= settings.ab_threshold:
            return False

    return True


def _count_gt_only(dd: list[str], sample_indices: list[int], counts: Counter[str]) -> None:
    for ci in sample_indices:
        counts[dd[ci].split(":", 1)[0]] += 1


def _count_gq_dp(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    indices: _FieldIndices,
    settings: _FilterSettings,
) -> None:
    min_gq, min_dp = settings.min_gq, settings.min_dp
    for ci in sample_indices:
        parts = dd[ci].split(":")
        if not _passes_gq_dp(parts, indices, min_gq=min_gq, min_dp=min_dp):
            continue
        counts[parts[0]] += 1


def _count_gq_dp_ab(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    indices: _FieldIndices,
    settings: _FilterSettings,
) -> None:
    for ci in sample_indices:
        parts = dd[ci].split(":")
        if not _passes_gq_dp_ab(parts, indices, settings):
            continue
        counts[parts[0]] += 1


def _count_samples(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    settings: _FilterSettings,
    indices: _FieldIndices | None,
) -> None:
    mode = settings.count_mode
    if mode == "gt_only":
        _count_gt_only(dd, sample_indices, counts)
    elif mode == "gq_dp":
        assert indices is not None
        _count_gq_dp(dd, sample_indices, counts, indices, settings)
    else:
        assert indices is not None
        _count_gq_dp_ab(dd, sample_indices, counts, indices, settings)


def _region_filter_settings(
    region: str,
    sex: str,
    *,
    min_gq: float | None,
    min_dp: int | None,
    ab_threshold: float | None,
    min_gq_nonpar: float | None,
    min_dp_nonpar: int | None,
    ab_threshold_nonpar: float | None,
) -> _FilterSettings:
    if sex == "male" and region == "noPar":
        return _FilterSettings(
            min_gq=min_gq if min_gq_nonpar is None else min_gq_nonpar,
            min_dp=min_dp if min_dp_nonpar is None else min_dp_nonpar,
            ab_threshold=ab_threshold if ab_threshold_nonpar is None else ab_threshold_nonpar,
        )
    return _FilterSettings(
        min_gq=min_gq,
        min_dp=min_dp,
        ab_threshold=ab_threshold,
    )


def is_snv(ref: str, alt: str) -> bool:
    """Return True only for SNVs: ``len(ref) == 1`` and ``len(alt) == 1``.

    Multi-allelic ``ALT`` fields (comma-separated) are excluded.
    """
    if "," in alt:
        return False
    return len(ref) == 1 and len(alt) == 1


def _chrom_key(chrom: str) -> str:
    c = chrom.strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


def chrom_matches(vcf_chrom: str, target_chrom: str) -> bool:
    """Return whether a VCF ``#CHROM`` value matches the requested chromosome."""
    return _chrom_key(vcf_chrom) == _chrom_key(target_chrom)


def compute_genotype_counts(
    vcf_path: str | Path,
    male_children: list[str],
    female_children: list[str],
    *,
    chromosome: str,
    allele_freqs: Mapping[str, float] | None = None,
    gnomad_af: str | Path | GnomadAfStore | None = None,
    common_freq_cutoff: float = 0.01,
    af_key_col: str = "id",
    min_gq: float | None = None,
    min_dp: int | None = None,
    ab_threshold: float | None = None,
    min_gq_nonpar: float | None = None,
    min_dp_nonpar: int | None = None,
    ab_threshold_nonpar: float | None = None,
    strict: bool = False,
    on_excluded: Callable[[list[str], list[str]], None] | None = None,
) -> GenotypeCountResult:
    """Scan *vcf_path* once and count genotypes for male and female children.

    Analysis is per chromosome: only rows whose ``#CHROM`` matches
    *chromosome* are included (``chr1`` and ``1`` are treated as equivalent).
    All input files are expected to be chromosome-specific.

    For chrX, counts are accumulated separately in three regions: ``Par1``,
    ``noPar``, and ``Par2``. Other chromosomes use a single ``all`` region.

    Only SNVs are included. Variants with frequency above *common_freq_cutoff*
    are skipped when *allele_freqs* or *gnomad_af* is provided.

    Sample FORMAT fields are assumed to list ``GT`` first. When quality filters
    are enabled, ``GQ``, ``DP``, and ``AD``/``AB`` are read by index from the
    FORMAT column.

    *allele_freqs* is a static variant-id -> AF map for this chromosome.
    *gnomad_af* is a base directory (or :class:`~sexxy.gnomad.GnomadAfStore`);
    the file ``{chromosome}/{chromosome}-common-af.json`` is loaded once.

    Per-genotype quality filters (each optional; filtering is enabled only
    when the parameter is set):

    *min_gq*
        Skip calls with genotype quality (``GQ``) below this value.
    *min_dp*
        Skip calls with read depth (``DP``) below this value.
    *ab_threshold*
        Allele-balance filter applied only to ``0/1`` and ``1/1`` genotypes.
        Require ``AB > ab_threshold``. ``AB`` is read from the sample field
        when present, otherwise computed from ``AD``.
    *min_gq_nonpar*, *min_dp_nonpar*, *ab_threshold_nonpar*
        Optional overrides used only for **male** calls in the chrX ``noPar``
        region. Each defaults to the corresponding global filter when unset.
        Female calls always use the global filters in all chrX regions.
    *strict*
        When ``True``, require every metadata child to appear in the VCF
        header; otherwise exclude children not present in the VCF.
    *on_excluded*
        Optional callback ``(excluded_male, excluded_female)`` invoked when
        children are dropped because they are absent from the VCF header.

    Returns a :class:`~sexxy.results.GenotypeCountResult`.
    """
    if allele_freqs is not None and gnomad_af is not None:
        raise ValueError("pass only one of allele_freqs or gnomad_af")

    gnomad_store: GnomadAfStore | None
    if gnomad_af is None:
        gnomad_store = None
    elif isinstance(gnomad_af, GnomadAfStore):
        gnomad_store = gnomad_af
    else:
        gnomad_store = GnomadAfStore(gnomad_af)

    samples = read_vcf_samples(vcf_path)
    cohort = filter_children_to_vcf(
        samples, male_children, female_children, strict=strict
    )
    if not cohort.male_children and not cohort.female_children:
        raise ValueError(
            "no male or female children from metadata are present in the VCF header"
        )
    if on_excluded is not None and (cohort.excluded_male or cohort.excluded_female):
        on_excluded(list(cohort.excluded_male), list(cohort.excluded_female))
    male_children = cohort.male_children
    female_children = cohort.female_children

    male_ind = sample_column_indices(samples, male_children) if male_children else []
    female_ind = sample_column_indices(samples, female_children) if female_children else []

    if is_chrx(chromosome):
        regions = CHRX_REGION_ORDER
        dgt_m = {r: Counter() for r in regions}
        dgt_f = {r: Counter() for r in regions}
    else:
        regions = ("all",)
        dgt_m = {"all": Counter()}
        dgt_f = {"all": Counter()}

    filter_kw = {
        "min_gq": min_gq,
        "min_dp": min_dp,
        "ab_threshold": ab_threshold,
        "min_gq_nonpar": min_gq_nonpar,
        "min_dp_nonpar": min_dp_nonpar,
        "ab_threshold_nonpar": ab_threshold_nonpar,
    }
    male_settings = {
        r: _region_filter_settings(r, "male", **filter_kw) for r in regions
    }
    female_settings = {
        r: _region_filter_settings(r, "female", **filter_kw) for r in regions
    }
    need_field_indices = any(
        male_settings[r].count_mode != "gt_only"
        or female_settings[r].count_mode != "gt_only"
        for r in regions
    )

    chrom_af: Mapping[str, float] | None = None
    if gnomad_store is not None:
        chrom_af = gnomad_store.for_chromosome(chromosome)

    target_chrom_key = _chrom_key(chromosome)
    scan_chrx = is_chrx(chromosome)
    use_variant_af_key = allele_freqs is not None and af_key_col == "variant"
    field_indices: _FieldIndices | None = None

    with _open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if _chrom_key(chrom) != target_chrom_key:
                continue
            if not is_snv(ref, alt):
                continue

            if scan_chrx:
                region = chrx_region(int(pos))
                if region is None:
                    continue
            else:
                region = "all"

            if chrom_af is not None:
                if float(chrom_af.get(vid, 0)) > common_freq_cutoff:
                    continue
            elif allele_freqs is not None:
                key = f"{chrom}:{pos}:{ref}:{alt}" if use_variant_af_key else vid
                if float(allele_freqs.get(key, 0)) > common_freq_cutoff:
                    continue

            if need_field_indices and field_indices is None:
                field_indices = _field_indices(parts[8])

            dd = parts[9:]
            if male_ind:
                _count_samples(dd, male_ind, dgt_m[region], male_settings[region], field_indices)
            if female_ind:
                _count_samples(dd, female_ind, dgt_f[region], female_settings[region], field_indices)

    return GenotypeCountResult(
        chromosome=chromosome,
        regions=regions,
        male={r: dict(dgt_m[r]) for r in regions},
        female={r: dict(dgt_f[r]) for r in regions},
        male_cohort_size=len(male_children),
        female_cohort_size=len(female_children),
        excluded_male=cohort.excluded_male,
        excluded_female=cohort.excluded_female,
    )
