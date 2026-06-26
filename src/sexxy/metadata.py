"""Load metadata and identify male/female children with valid parent IDs."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

from sexxy.table import read_table

# Values treated as missing parent IDs.
_MISSING = frozenset({None, "", ".", "NA", "na", "NaN", "nan", "0", 0})


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def _is_valid_id(value) -> bool:
    if _is_missing(value):
        return False
    return str(value).strip() not in _MISSING


def _normalize_sex(value) -> str | None:
    """Map metadata sex labels to ``'male'`` or ``'female'``.

    Supported labels:
    - Text: ``Male``, ``Female``, ``male``, ``female``, ``M``, ``F``
    - Numeric codes: ``1`` (male), ``2`` (female), including ``1.0`` / ``2.0``
    """
    if _is_missing(value):
        return None

    if isinstance(value, bool):
        return None

    if isinstance(value, int):
        if value == 1:
            return "male"
        if value == 2:
            return "female"

    if isinstance(value, float):
        if value == 1.0:
            return "male"
        if value == 2.0:
            return "female"

    s = str(value).strip().lower()
    if s in {"1", "1.0"}:
        return "male"
    if s in {"2", "2.0"}:
        return "female"
    if s in {"m", "male", "xy"}:
        return "male"
    if s in {"f", "female", "xx"}:
        return "female"
    return None


def load_children_by_sex(
    metadata_path: str,
    *,
    patient_col: str = "patient_id",
    father_col: str = "father_id",
    mother_col: str = "mother_id",
    sex_col: str = "sex",
    sep: str | None = None,
) -> tuple[list[str], list[str], list[dict[str, str]]]:
    """Load metadata and return (male_children, female_children, children_rows).

    Children are rows with valid father and mother IDs. Sex must be recognized
    as male or female; rows with unknown sex are skipped.

    Sex labels may be text (``Male``/``Female``) or numeric codes (``1``/``2``).
    """
    _fieldnames, rows = read_table(metadata_path, sep=sep)
    for col in (patient_col, father_col, mother_col, sex_col):
        if col not in _fieldnames:
            raise ValueError(f"metadata missing required column: {col!r}")

    male_children: list[str] = []
    female_children: list[str] = []
    children_rows: list[dict[str, str]] = []

    for row in rows:
        if not (_is_valid_id(row[father_col]) and _is_valid_id(row[mother_col])):
            continue
        sex_norm = _normalize_sex(row[sex_col])
        children_rows.append(row)
        if sex_norm == "male":
            male_children.append(row[patient_col])
        elif sex_norm == "female":
            female_children.append(row[patient_col])

    return male_children, female_children, children_rows


def sample_column_indices(
    vcf_samples: Sequence[str],
    sample_ids: Sequence[str],
) -> list[int]:
    """Map sample IDs to 0-based indices into the VCF genotype columns."""
    index_by_id: Mapping[str, int] = {s: i for i, s in enumerate(vcf_samples)}
    missing = [sid for sid in sample_ids if sid not in index_by_id]
    if missing:
        raise ValueError(
            f"{len(missing)} sample(s) not found in VCF header, e.g. {missing[:5]}"
        )
    return [index_by_id[sid] for sid in sample_ids]
