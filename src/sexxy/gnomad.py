"""Load gnomAD v4 per-chromosome common allele-frequency JSON files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

DEFAULT_GNOMAD_AF_DIR = "/mnt/home/mbershadsky/ceph/gnomad.v4"


def gnomad_af_path(base_dir: str | Path, chrm: str) -> Path:
    """Return path to ``{base_dir}/{chrm}/{chrm}-common-af.json``."""
    base = Path(base_dir)
    return base / chrm / f"{chrm}-common-af.json"


def load_gnomad_af_json(path: str | Path) -> dict[str, float]:
    """Load a gnomAD common-AF JSON file (variant id -> frequency)."""
    with open(path, "r") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}, got {type(data).__name__}")
    return data


class GnomadAfStore:
    """Lazy per-chromosome gnomAD AF lookup backed by JSON files.

    Files are expected at::

        {base_dir}/{chrm}/{chrm}-common-af.json

    Each file is a JSON object mapping variant IDs (VCF ID column) to AF.
    """

    def __init__(self, base_dir: str | Path = DEFAULT_GNOMAD_AF_DIR):
        self.base_dir = Path(base_dir)
        self._by_chrom: dict[str, Mapping[str, float]] = {}

    def path_for(self, chrm: str) -> Path:
        return gnomad_af_path(self.base_dir, chrm)

    def for_chromosome(self, chrm: str) -> Mapping[str, float]:
        """Load and cache the AF map for *chrm*."""
        if chrm not in self._by_chrom:
            gfile = self.path_for(chrm)
            self._by_chrom[chrm] = load_gnomad_af_json(gfile)
        return self._by_chrom[chrm]

    def get(self, chrm: str, variant_id: str) -> float:
        return float(self.for_chromosome(chrm).get(variant_id, 0))
