"""Result containers and output helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from sexxy.chrx import CHRX_REGION_ORDER, is_chrx

MALE_OUTPUT_KEYS = ("chromosome", "sex", "male_children", "gt_counts")
FEMALE_OUTPUT_KEYS = ("chromosome", "sex", "female_children", "gt_counts")
CHRX_MALE_OUTPUT_KEYS = ("chromosome", "region", "sex", "male_children", "gt_counts")
CHRX_FEMALE_OUTPUT_KEYS = ("chromosome", "region", "sex", "female_children", "gt_counts")


@dataclass(frozen=True)
class GenotypeCountResult:
    """Genotype counts by sex and, for chrX, by pseudoautosomal region."""

    chromosome: str
    regions: tuple[str, ...]
    male: dict[str, dict[str, int]]
    female: dict[str, dict[str, int]]
    male_cohort_size: int
    female_cohort_size: int
    excluded_male: tuple[str, ...] = ()
    excluded_female: tuple[str, ...] = ()

    def male_counts(self, region: str | None = None) -> dict[str, int]:
        return self._counts(self.male, region)

    def female_counts(self, region: str | None = None) -> dict[str, int]:
        return self._counts(self.female, region)

    def _counts(self, by_region: dict[str, dict[str, int]], region: str | None) -> dict[str, int]:
        if region is not None:
            return by_region[region]
        if len(self.regions) == 1:
            return by_region[self.regions[0]]
        raise ValueError(f"region required for chrX results; choose from {self.regions}")


def _ensure_parent_dir(path: Path) -> None:
    parent = path.parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)


def resolve_output_target(
    output: str | Path | None,
    output_dir: str | Path | None,
    chromosome: str,
) -> str | None:
    """Combine *output* basename/prefix with *output_dir*, creating the directory."""
    if output_dir is None:
        return str(output) if output is not None else None

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)

    stem = f"counts.{chromosome}"
    if output is not None:
        p = Path(output)
        stem = p.stem if p.suffix == ".json" else p.name

    return str(directory / stem)


def output_prefix(path: str | Path | None, chromosome: str) -> str:
    if path is None:
        return f"counts.{chromosome}"
    p = Path(path)
    if p.suffix == ".json":
        return str(p.with_suffix(""))
    return str(p)


def _sex_payload(
    result: GenotypeCountResult,
    sex: str,
    *,
    male_children: int,
    female_children: int,
    region: str | None = None,
) -> dict:
    if region is not None:
        counts = result.male_counts(region) if sex == "male" else result.female_counts(region)
    else:
        counts = result.male_counts() if sex == "male" else result.female_counts()

    payload: dict = {
        "chromosome": result.chromosome,
        "sex": sex,
        "gt_counts": counts,
    }
    if region is not None:
        payload["region"] = region
    if sex == "male":
        payload["male_children"] = male_children
    else:
        payload["female_children"] = female_children
    return payload


def _write_sex_file(
    written: list[Path],
    prefix: str,
    sex: str,
    payload: dict,
    *,
    suffix: str = ".json",
) -> None:
    path = Path(f"{prefix}.{sex}{suffix}")
    _ensure_parent_dir(path)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    written.append(path)


def write_genotype_count_results(
    result: GenotypeCountResult,
    output: str | Path | None,
    *,
    male_children: int,
    female_children: int,
) -> list[Path]:
    """Write result JSON file(s). Returns paths written.

    Autosomes and chrY: ``{prefix}.male.json`` and ``{prefix}.female.json``.

    chrX: six files ``{prefix}.{sex}.{region}.json`` for ``Par1``, ``noPar``,
    and ``Par2``. Each file lists only the cohort count for that sex; there is
    no ``region`` field on autosomes/chrY.
    """
    written: list[Path] = []
    prefix = output_prefix(output, result.chromosome)

    if is_chrx(result.chromosome):
        for region in CHRX_REGION_ORDER:
            for sex in ("male", "female"):
                path = Path(f"{prefix}.{sex}.{region}.json")
                _ensure_parent_dir(path)
                payload = _sex_payload(
                    result,
                    sex,
                    male_children=male_children,
                    female_children=female_children,
                    region=region,
                )
                path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
                written.append(path)
        return written

    for sex in ("male", "female"):
        payload = _sex_payload(
            result,
            sex,
            male_children=male_children,
            female_children=female_children,
        )
        _write_sex_file(written, prefix, sex, payload)

    return written


def result_to_json(result: GenotypeCountResult, *, male_children: int, female_children: int) -> str:
    """Serialize male and female payloads as one JSON object (for stdout)."""
    if is_chrx(result.chromosome):
        payload = {
            "chromosome": result.chromosome,
            "regions": {
                region: {
                    "male": _sex_payload(
                        result, "male",
                        male_children=male_children,
                        female_children=female_children,
                        region=region,
                    ),
                    "female": _sex_payload(
                        result, "female",
                        male_children=male_children,
                        female_children=female_children,
                        region=region,
                    ),
                }
                for region in CHRX_REGION_ORDER
            },
        }
    else:
        payload = {
            "chromosome": result.chromosome,
            "male": _sex_payload(
                result, "male",
                male_children=male_children,
                female_children=female_children,
            ),
            "female": _sex_payload(
                result, "female",
                male_children=male_children,
                female_children=female_children,
            ),
        }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def write_run_params(
    output: str | Path | None,
    chromosome: str,
    params: dict,
) -> Path:
    """Write run parameters to ``{prefix}.params.json`` and return the path."""
    prefix = output_prefix(output, chromosome)
    path = Path(f"{prefix}.params.json")
    _ensure_parent_dir(path)
    payload = dict(params)
    payload["params_file"] = str(path)
    output_files = list(payload.get("output_files", []))
    output_files.append(str(path))
    payload["output_files"] = output_files
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path
