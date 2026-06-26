"""Load metadata and identify male/female children with valid parent IDs."""

from __future__ import annotations

from typing import Mapping, Sequence

import pandas as pd

# Values treated as missing parent IDs.
_MISSING = frozenset({None, "", ".", "NA", "na", "NaN", "nan", "0", 0})


def _is_valid_id(value) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    return str(value).strip() not in _MISSING


def _normalize_sex(value) -> str | None:
    """Map metadata sex labels to ``'male'`` or ``'female'``.

    Supported labels:
    - Text: ``Male``, ``Female``, ``male``, ``female``, ``M``, ``F``
    - Numeric codes: ``1`` (male), ``2`` (female), including ``1.0`` / ``2.0``
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
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
) -> tuple[list[str], list[str], pd.DataFrame]:
    """Load metadata and return (male_children, female_children, children_df).

    Children are rows with valid father and mother IDs. Sex must be recognized
    as male or female; rows with unknown sex are skipped.

    Sex labels may be text (``Male``/``Female``) or numeric codes (``1``/``2``).
    """
    df = pd.read_csv(metadata_path, sep=sep, dtype=str, keep_default_na=False)
    for col in (patient_col, father_col, mother_col, sex_col):
        if col not in df.columns:
            raise ValueError(f"metadata missing required column: {col!r}")

    has_parents = df[father_col].map(_is_valid_id) & df[mother_col].map(_is_valid_id)
    children = df.loc[has_parents].copy()
    children["_sex_norm"] = children[sex_col].map(_normalize_sex)

    male_children = children.loc[children["_sex_norm"] == "male", patient_col].tolist()
    female_children = children.loc[children["_sex_norm"] == "female", patient_col].tolist()
    children = children.drop(columns=["_sex_norm"])

    return male_children, female_children, children


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
