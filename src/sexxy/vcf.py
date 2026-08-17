"""Single-pass VCF scanning for SNV genotype counts by sex."""

from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, IO, Mapping

from sexxy.chrx import CHRX_REGION_ORDER, chrx_region, is_chrx
from sexxy.gnomad import GnomadAfStore, gnomad_af_key
from sexxy.intervals import ExcludeIntervals, load_exclude_intervals
from sexxy.metadata import filter_children_to_vcf, sample_column_indices
from sexxy.results import GenotypeCountResult


def get_n_fields(line: str, n: int) -> list[str]:
    """Return the first *n* tab-separated fields from a VCF data line."""
    parts = line.rstrip("\n").split("\t", n)
    if len(parts) < n:
        raise ValueError(f"expected at least {n} fields, got {len(parts)}")
    return parts[:n]


def _split_vcf_prefix(line: str) -> tuple[str, str, str, str, str, str] | None:
    """Return ``CHROM, POS, ID, REF, ALT, remainder``, or ``None`` if too few fields.

    Uses a capped split so sample columns are not parsed until a site is kept.
    """
    parts = line.split("\t", 5)
    if len(parts) < 6:
        return None
    return parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]


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


def af_variant_key(chrom: str, pos: str | int, ref: str, alt: str) -> str:
    """Return ``chrom:pos:ref:alt`` AF lookup key for one allele."""
    return f"{chrom}:{pos}:{ref}:{alt}"


# Common allele-1 genotypes: identity after unphasing. Avoids the remap scan
# on ordinary biallelic SNVs (``allele_number == 1``).
_ALLELE1_GT: dict[str, str] = {
    "0/0": "0/0",
    "0/1": "0/1",
    "1/0": "1/0",
    "1/1": "1/1",
    "0|0": "0/0",
    "0|1": "0/1",
    "1|0": "1/0",
    "1|1": "1/1",
    "0": "0",
    "1": "1",
    ".": ".",
    "./.": "./.",
    ".|.": "./.",
    "": ".",
    "0/.": "0/.",
    "./0": "./0",
    "1/.": "1/.",
    "./1": "./1",
    "0|.": "0/.",
    ".|0": "./0",
    "1|.": "1/.",
    ".|1": "./1",
}


def _remap_gt_scan(gt: str, target: str) -> str | None:
    """Remap *gt* so VCF allele *target* becomes ``1``; other alts → ``None``."""
    gt_u = _unphased_gt(gt)
    if gt_u in (".", "./.", ""):
        return gt_u if gt_u else "."
    remapped: list[str] = []
    for allele in gt_u.split("/"):
        if allele == ".":
            remapped.append(".")
        elif allele == "0":
            remapped.append("0")
        elif allele == target:
            remapped.append("1")
        else:
            return None
    return "/".join(remapped)


def _remap_gt(gt: str, allele_number: int) -> str | None:
    """Remap *gt* so VCF allele *allele_number* becomes ``1``.

    Other non-ref alleles make the call unusable for this allele (``None``).
    Phased genotypes are converted to unphased ``/`` form.
    """
    if allele_number == 1:
        cached = _ALLELE1_GT.get(gt)
        if cached is not None:
            return cached
        return _remap_gt_scan(gt, "1")
    return _remap_gt_scan(gt, str(allele_number))


def _parse_sample_fields(format_str: str, sample_field: str) -> dict[str, str]:
    keys = format_str.split(":")
    vals = sample_field.split(":")
    if len(vals) < len(keys):
        vals = vals + ["."] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def _allele_balance_from_ad(ad: str, allele_index: int) -> float | None:
    """Return ``AD[allele_index] / sum(AD)``, or ``None`` if unusable."""
    if ad in (".", ""):
        return None
    ad_parts = ad.split(",")
    if len(ad_parts) == 1:
        try:
            int(ad_parts[0])
        except ValueError:
            return None
        return 1.0
    if allele_index < 0 or allele_index >= len(ad_parts):
        return None
    try:
        depths = [int(x) for x in ad_parts]
    except ValueError:
        return None
    total = sum(depths)
    if total == 0:
        return None
    return depths[allele_index] / total


def _allele_balance(fields: Mapping[str, str], *, allele_index: int = 1) -> float | None:
    """Allele balance from ``AD`` as ``AD[allele_index] / sum(AD)``."""
    return _allele_balance_from_ad(fields.get("AD", "."), allele_index)


def _passes_genotype_filters(
    fields: Mapping[str, str],
    *,
    min_gq: float | None,
    min_dp: int | None,
    ab_threshold: float | None,
    allele_index: int = 1,
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
        ab = _allele_balance(fields, allele_index=allele_index)
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
    allele_number: int = 1,
) -> None:
    fields = _parse_sample_fields(format_str, sample_field)
    remapped = _remap_gt(fields.get("GT", "."), allele_number)
    if remapped is None:
        return
    fields = dict(fields)
    fields["GT"] = remapped
    if not _passes_genotype_filters(
        fields,
        min_gq=min_gq,
        min_dp=min_dp,
        ab_threshold=ab_threshold,
        allele_index=allele_number,
    ):
        return
    counts[remapped] += 1


@dataclass(frozen=True)
class _FieldIndices:
    """Colon-field indices for FORMAT tags (GT is always index 0)."""

    gq: int | None
    dp: int | None
    ad: int | None


def _count_mode(
    min_gq: float | None, min_dp: int | None, ab_threshold: float | None
) -> str:
    if min_gq is None and min_dp is None and ab_threshold is None:
        return "gt_only"
    if ab_threshold is not None:
        return "gq_dp_ab"
    return "gq_dp"


@dataclass(frozen=True)
class _FilterSettings:
    min_gq: float | None
    min_dp: int | None
    ab_threshold: float | None
    count_mode: str


def _field_indices(format_str: str) -> _FieldIndices:
    keys = format_str.split(":")
    idx = {k: i for i, k in enumerate(keys)}
    return _FieldIndices(
        gq=idx.get("GQ"),
        dp=idx.get("DP"),
        ad=idx.get("AD"),
    )


def _allele_balance_from_parts(
    parts: list[str],
    indices: _FieldIndices,
    *,
    allele_index: int = 1,
) -> float | None:
    if indices.ad is None or indices.ad >= len(parts):
        return None
    return _allele_balance_from_ad(parts[indices.ad], allele_index)


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
    *,
    allele_index: int = 1,
) -> bool:
    gt = _unphased_gt(parts[0] if parts else ".")
    if gt in (".", "./."):
        return False

    if not _passes_gq_dp(parts, indices, min_gq=settings.min_gq, min_dp=settings.min_dp):
        return False

    if settings.ab_threshold is not None and gt in ("0/1", "1/1"):
        ab = _allele_balance_from_parts(parts, indices, allele_index=allele_index)
        if ab is None or ab <= settings.ab_threshold:
            return False

    return True


def _with_remapped_gt(parts: list[str], allele_number: int) -> list[str] | None:
    if not parts:
        return None
    remapped = _remap_gt(parts[0], allele_number)
    if remapped is None:
        return None
    if remapped == parts[0]:
        return parts
    out = list(parts)
    out[0] = remapped
    return out


def _count_gt_only(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    *,
    allele_number: int = 1,
) -> None:
    if allele_number == 1:
        lookup = _ALLELE1_GT
        for ci in sample_indices:
            gt = dd[ci].split(":", 1)[0]
            remapped = lookup.get(gt)
            if remapped is None:
                remapped = _remap_gt_scan(gt, "1")
                if remapped is None:
                    continue
            counts[remapped] += 1
        return
    target = str(allele_number)
    for ci in sample_indices:
        remapped = _remap_gt_scan(dd[ci].split(":", 1)[0], target)
        if remapped is None:
            continue
        counts[remapped] += 1


def _count_gq_dp(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    indices: _FieldIndices | None,
    settings: _FilterSettings,
    *,
    allele_number: int = 1,
) -> None:
    if indices is None:
        return
    min_gq, min_dp = settings.min_gq, settings.min_dp
    for ci in sample_indices:
        parts = _with_remapped_gt(dd[ci].split(":"), allele_number)
        if parts is None:
            continue
        if not _passes_gq_dp(parts, indices, min_gq=min_gq, min_dp=min_dp):
            continue
        counts[parts[0]] += 1


def _count_gq_dp_ab(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    indices: _FieldIndices | None,
    settings: _FilterSettings,
    *,
    allele_number: int = 1,
) -> None:
    if indices is None:
        return
    for ci in sample_indices:
        parts = _with_remapped_gt(dd[ci].split(":"), allele_number)
        if parts is None:
            continue
        if not _passes_gq_dp_ab(parts, indices, settings, allele_index=allele_number):
            continue
        counts[parts[0]] += 1


def _count_samples(
    dd: list[str],
    sample_indices: list[int],
    counts: Counter[str],
    settings: _FilterSettings,
    indices: _FieldIndices | None,
    *,
    allele_number: int = 1,
) -> None:
    mode = settings.count_mode
    if mode == "gt_only":
        _count_gt_only(dd, sample_indices, counts, allele_number=allele_number)
    elif mode == "gq_dp":
        assert indices is not None
        _count_gq_dp(dd, sample_indices, counts, indices, settings, allele_number=allele_number)
    else:
        assert indices is not None
        _count_gq_dp_ab(dd, sample_indices, counts, indices, settings, allele_number=allele_number)


class _IndicesRef:
    """Mutable FORMAT-index holder so bound counters see the first parsed FORMAT."""

    __slots__ = ("value",)

    def __init__(self) -> None:
        self.value: _FieldIndices | None = None


def _noop_count(
    dd: list[str], counts: Counter[str], allele_number: int = 1
) -> None:
    return


def _bind_count(
    sample_indices: list[int],
    settings: _FilterSettings,
    indices_ref: _IndicesRef,
):
    """Return ``(dd, counts, allele_number) -> None`` with mode/settings closed over."""
    if not sample_indices:
        return _noop_count
    mode = settings.count_mode
    if mode == "gt_only":
        def _count(dd: list[str], counts: Counter[str], allele_number: int = 1) -> None:
            _count_gt_only(dd, sample_indices, counts, allele_number=allele_number)

        return _count
    if mode == "gq_dp":
        def _count(dd: list[str], counts: Counter[str], allele_number: int = 1) -> None:
            _count_gq_dp(
                dd,
                sample_indices,
                counts,
                indices_ref.value,
                settings,
                allele_number=allele_number,
            )

        return _count

    def _count(dd: list[str], counts: Counter[str], allele_number: int = 1) -> None:
        _count_gq_dp_ab(
            dd,
            sample_indices,
            counts,
            indices_ref.value,
            settings,
            allele_number=allele_number,
        )

    return _count


def _bind_father_count(
    pairs: list[tuple[int, int]],
    settings: _FilterSettings,
    indices_ref: _IndicesRef,
):
    def _count(dd: list[str], counts: Counter[str], allele_number: int = 1) -> None:
        _count_samples_with_father_match(
            dd,
            pairs,
            counts,
            settings,
            indices_ref.value,
            allele_number=allele_number,
        )

    return _count


_FATHER_MATCH_GTS = frozenset({"0/0", "0/1", "1/1"})


def _child_gt_if_passes(
    dd: list[str],
    ci: int,
    settings: _FilterSettings,
    indices: _FieldIndices | None,
    *,
    allele_number: int = 1,
) -> str | None:
    parts = _with_remapped_gt(dd[ci].split(":"), allele_number)
    if parts is None:
        return None
    mode = settings.count_mode
    if mode == "gt_only":
        return parts[0]
    assert indices is not None
    if mode == "gq_dp":
        if not _passes_gq_dp(parts, indices, min_gq=settings.min_gq, min_dp=settings.min_dp):
            return None
    elif not _passes_gq_dp_ab(parts, indices, settings, allele_index=allele_number):
        return None
    return parts[0]


def _count_samples_with_father_match(
    dd: list[str],
    child_father_pairs: list[tuple[int, int]],
    counts: Counter[str],
    settings: _FilterSettings,
    indices: _FieldIndices | None,
    *,
    allele_number: int = 1,
) -> None:
    for child_ci, father_ci in child_father_pairs:
        child_gt = _child_gt_if_passes(
            dd, child_ci, settings, indices, allele_number=allele_number
        )
        if child_gt is None:
            continue
        counts[child_gt] += 1
        child_gt_u = _unphased_gt(child_gt)
        if child_gt_u in _FATHER_MATCH_GTS:
            father_gt = _remap_gt(dd[father_ci].split(":", 1)[0], allele_number)
            if father_gt is not None and _unphased_gt(father_gt) == child_gt_u:
                counts[f"{child_gt_u}_f"] += 1


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
        min_gq = min_gq if min_gq_nonpar is None else min_gq_nonpar
        min_dp = min_dp if min_dp_nonpar is None else min_dp_nonpar
        ab_threshold = ab_threshold if ab_threshold_nonpar is None else ab_threshold_nonpar
    return _FilterSettings(
        min_gq=min_gq,
        min_dp=min_dp,
        ab_threshold=ab_threshold,
        count_mode=_count_mode(min_gq, min_dp, ab_threshold),
    )


def is_snv(ref: str, alt: str) -> bool:
    """Return True only for single-allele SNVs: ``len(ref) == 1`` and ``len(alt) == 1``.

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


def is_chry(chromosome: str) -> bool:
    return _chrom_key(chromosome).upper() == "Y"


def _male_child_father_pairs(
    samples: list[str],
    male_children: list[str],
    male_father_by_child: Mapping[str, str],
) -> list[tuple[int, int]]:
    sample_index = {s: i for i, s in enumerate(samples)}
    pairs: list[tuple[int, int]] = []
    for child_id in male_children:
        father_id = male_father_by_child.get(child_id)
        if not father_id or father_id not in sample_index:
            continue
        pairs.append((sample_index[child_id], sample_index[father_id]))
    return pairs


def compute_genotype_counts(
    vcf_path: str | Path,
    male_children: list[str],
    female_children: list[str],
    *,
    chromosome: str,
    allele_freqs: Mapping[str, float] | None = None,
    gnomad_af: str | Path | GnomadAfStore | None = None,
    common_freq_cutoff: float = 0.01,
    min_gq: float | None = None,
    min_dp: int | None = None,
    ab_threshold: float | None = None,
    min_gq_nonpar: float | None = None,
    min_dp_nonpar: int | None = None,
    ab_threshold_nonpar: float | None = None,
    exclude_repeats: str | Path | ExcludeIntervals | None = None,
    strict: bool = False,
    on_excluded: Callable[[list[str], list[str]], None] | None = None,
    check_father: bool = False,
    male_father_by_child: Mapping[str, str] | None = None,
    include_multiallelic: bool = False,
) -> GenotypeCountResult:
    """Scan *vcf_path* once and count genotypes for male and female children.

    Analysis is per chromosome: only rows whose ``#CHROM`` matches
    *chromosome* are included (``chr1`` and ``1`` are treated as equivalent).
    All input files are expected to be chromosome-specific.

    For chrX, counts are accumulated separately in three regions: ``Par1``,
    ``noPar``, and ``Par2``. Other chromosomes use a single ``all`` region.

    Only SNV alleles are included (``len(REF) == 1`` and ``len(ALT allele) == 1``).
    By default multi-allelic rows (comma-separated ``ALT``) are skipped. When
    *include_multiallelic* is ``True``, each ALT allele is processed separately:
    genotypes are remapped so that allele becomes ``1``, and each passing allele
    contributes its own genotype counts.

    Alleles with frequency above *common_freq_cutoff* are skipped using only
    the site fields (``POS``, ``REF``, ``ALT``), before sample columns are
    parsed. *gnomad_af* maps are keyed by ``chrN_pos_ref_alt``; *allele_freqs* maps
    are keyed by ``chrom:pos:ref:alt``. Missing keys are treated as AF ``0``
    (kept). Input VCFs are chromosome-specific; non-matching contigs are
    counted in *skipped_contigs* and not AF-filtered. The number of SNV alleles
    skipped as common is returned in *skipped_common_variants*.

    Sample FORMAT fields are assumed to list ``GT`` first. When quality filters
    are enabled, ``GQ``, ``DP``, and ``AD`` are read by index from the FORMAT
    column.

    *allele_freqs* is a static ``chrom:pos:ref:alt`` -> AF map for this chromosome.
    *gnomad_af* is a base directory (or :class:`~sexxy.gnomad.GnomadAfStore`);
    the file ``{chromosome}/{chromosome}-common-af.json`` is loaded once
    (keys ``chrN_pos_ref_alt``).

    Per-genotype quality filters (each optional; filtering is enabled only
    when the parameter is set):

    *min_gq*
        Skip calls with genotype quality (``GQ``) below this value.
    *min_dp*
        Skip calls with read depth (``DP``) below this value.
    *ab_threshold*
        Allele-balance filter applied only to remapped ``0/1`` and ``1/1``
        genotypes. Require ``AD[i] / sum(AD) > ab_threshold`` for the allele
        being processed (``i`` is the VCF allele index).
    *min_gq_nonpar*, *min_dp_nonpar*, *ab_threshold_nonpar*
        Optional overrides used only for **male** calls in the chrX ``noPar``
        region. Each defaults to the corresponding global filter when unset.
        Female calls always use the global filters in all chrX regions.
    *exclude_repeats*
        Optional path to a tab-separated file of ``chrom start end`` intervals,
        or a pre-loaded :class:`~sexxy.intervals.ExcludeIntervals`. Variants
        with ``start <= POS < end`` are skipped.
    *strict*
        When ``True``, require every metadata child to appear in the VCF
        header; otherwise exclude children not present in the VCF.
    *on_excluded*
        Optional callback ``(excluded_male, excluded_female)`` invoked when
        children are dropped because they are absent from the VCF header.
    *check_father*
        When ``True`` on chrY, for each male child also count how often
        ``0/0``, ``0/1``, and ``1/1`` match the father's genotype under
        keys ``0/0_f``, ``0/1_f``, and ``1/1_f``. Child and father genotypes
        are remapped to the allele being processed before comparison.
    *male_father_by_child*
        Child sample ID to father sample ID map (required when
        *check_father* is ``True``).
    *include_multiallelic*
        When ``True``, expand comma-separated ``ALT`` fields and process each
        SNV allele independently.

    Returns a :class:`~sexxy.results.GenotypeCountResult`.
    """
    if allele_freqs is not None and gnomad_af is not None:
        raise ValueError("pass only one of allele_freqs or gnomad_af")

    repeat_filter: ExcludeIntervals | None
    if exclude_repeats is None:
        repeat_filter = None
    elif isinstance(exclude_repeats, ExcludeIntervals):
        repeat_filter = exclude_repeats
    else:
        repeat_filter = load_exclude_intervals(exclude_repeats, chromosome)
    exclude_repeat_intervals = len(repeat_filter) if repeat_filter is not None else 0

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

    use_father_check = (
        check_father
        and is_chry(chromosome)
        and male_father_by_child
        and male_children
    )
    male_father_pairs = (
        _male_child_father_pairs(samples, male_children, male_father_by_child)
        if use_father_check
        else []
    )

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

    def allele_passes_af(chrom: str, pos: str, ref: str, allele: str) -> bool:
        if chrom_af is not None:
            return not (
                float(chrom_af.get(gnomad_af_key(chrom, pos, ref, allele), 0))
                > common_freq_cutoff
            )
        if allele_freqs is not None:
            return not (
                float(allele_freqs.get(af_variant_key(chrom, pos, ref, allele), 0))
                > common_freq_cutoff
            )
        return True

    target_chrom_key = _chrom_key(chromosome)
    scan_chrx = is_chrx(chromosome)
    indices_ref = _IndicesRef()
    skipped_contigs: Counter[str] = Counter()
    excluded_repeat_rows = 0
    skipped_common_variants = 0

    if male_father_pairs:
        count_male = {
            r: _bind_father_count(male_father_pairs, male_settings[r], indices_ref)
            for r in regions
        }
    else:
        count_male = {
            r: _bind_count(male_ind, male_settings[r], indices_ref) for r in regions
        }
    count_female = {
        r: _bind_count(female_ind, female_settings[r], indices_ref) for r in regions
    }

    with _open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            prefix = _split_vcf_prefix(line)
            if prefix is None:
                continue
            chrom, pos, _vid, ref, alt, rest = prefix
            if chrom != chromosome and _chrom_key(chrom) != target_chrom_key:
                skipped_contigs[chrom] += 1
                continue
            if len(ref) != 1:
                continue

            if "," not in alt:
                if len(alt) != 1:
                    continue
                if not allele_passes_af(chrom, pos, ref, alt):
                    skipped_common_variants += 1
                    continue
                keep_alleles: tuple[int, ...] | list[int] = (1,)
            else:
                if not include_multiallelic:
                    continue
                keep_list: list[int] = []
                for allele_idx, allele in enumerate(alt.split(",")):
                    if len(allele) != 1:
                        continue
                    if not allele_passes_af(chrom, pos, ref, allele):
                        skipped_common_variants += 1
                        continue
                    keep_list.append(allele_idx + 1)
                if not keep_list:
                    continue
                keep_alleles = keep_list

            pos_i = int(pos)
            if repeat_filter is not None and repeat_filter.contains(pos_i):
                excluded_repeat_rows += 1
                continue

            if scan_chrx:
                region = chrx_region(pos_i)
                if region is None:
                    continue
            else:
                region = "all"

            tail = rest.rstrip("\n").split("\t")
            if len(tail) < 5:
                continue
            if need_field_indices and indices_ref.value is None:
                indices_ref.value = _field_indices(tail[3])
            dd = tail[4:]

            count_m = count_male[region]
            count_f = count_female[region]
            dgt_mr = dgt_m[region]
            dgt_fr = dgt_f[region]
            for allele_number in keep_alleles:
                count_m(dd, dgt_mr, allele_number)
                count_f(dd, dgt_fr, allele_number)

    return GenotypeCountResult(
        chromosome=chromosome,
        regions=regions,
        male={r: dict(dgt_m[r]) for r in regions},
        female={r: dict(dgt_f[r]) for r in regions},
        male_cohort_size=len(male_children),
        female_cohort_size=len(female_children),
        excluded_male=cohort.excluded_male,
        excluded_female=cohort.excluded_female,
        skipped_contigs=dict(skipped_contigs),
        excluded_repeat_rows=excluded_repeat_rows,
        exclude_repeat_intervals=exclude_repeat_intervals,
        skipped_common_variants=skipped_common_variants,
    )
