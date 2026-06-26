"""Result containers and output helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from sexxy.chrx import CHRX_REGION_ORDER, is_chrx


@dataclass(frozen=True)
class GenotypeCountResult:
    """Genotype counts by sex and, for chrX, by pseudoautosomal region."""

    chromosome: str
    regions: tuple[str, ...]
    male: dict[str, dict[str, int]]
    female: dict[str, dict[str, int]]

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


def single_output_path(output: str | Path | None, chromosome: str) -> Path:
    if output is None:
        return Path(f"counts.{chromosome}.json")
    p = Path(output)
    if p.suffix == ".json":
        return p
    return p.with_suffix(".json")


def write_genotype_count_results(
    result: GenotypeCountResult,
    output: str | Path | None,
    *,
    male_children: int,
    female_children: int,
) -> list[Path]:
    """Write result JSON file(s). Returns paths written."""
    written: list[Path] = []

    if is_chrx(result.chromosome):
        prefix = output_prefix(output, result.chromosome)
        for region in CHRX_REGION_ORDER:
            for sex, counts, n_children in (
                ("male", result.male[region], male_children),
                ("female", result.female[region], female_children),
            ):
                path = Path(f"{prefix}.{sex}.{region}.json")
                _ensure_parent_dir(path)
                payload = {
                    "chromosome": result.chromosome,
                    "region": region,
                    "sex": sex,
                    "children": n_children,
                    "gt_counts": counts,
                }
                path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
                written.append(path)
        return written

    path = single_output_path(output, result.chromosome)
    _ensure_parent_dir(path)
    payload = {
        "chromosome": result.chromosome,
        "male_children": male_children,
        "female_children": female_children,
        "male_gt_counts": result.male_counts(),
        "female_gt_counts": result.female_counts(),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    written.append(path)
    return written


def result_to_json(result: GenotypeCountResult, *, male_children: int, female_children: int) -> str:
    """Serialize a non-chrX result, or full chrX result, as JSON text."""
    if is_chrx(result.chromosome):
        payload = {
            "chromosome": result.chromosome,
            "male_children": male_children,
            "female_children": female_children,
            "male_gt_counts": result.male,
            "female_gt_counts": result.female,
        }
    else:
        payload = {
            "chromosome": result.chromosome,
            "male_children": male_children,
            "female_children": female_children,
            "male_gt_counts": result.male_counts(),
            "female_gt_counts": result.female_counts(),
        }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
