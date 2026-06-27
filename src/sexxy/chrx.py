"""chrX pseudoautosomal and non-PAR region boundaries."""

from __future__ import annotations

CHRX_REGIONS: dict[str, tuple[int, int]] = {
    "Par1": (10_001, 2_781_479),
    "noPar": (2_781_489, 155_701_382),
    "Par2": (155_701_383, 156_030_895),
}

CHRX_REGION_ORDER: tuple[str, ...] = ("Par1", "noPar", "Par2")


def is_chrx(chromosome: str) -> bool:
    chrom = chromosome.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper() == "X"


def chrx_region(pos: int) -> str | None:
    """Return the chrX region name for *pos*, or None if outside all regions."""
    for name in CHRX_REGION_ORDER:
        start, end = CHRX_REGIONS[name]
        if start <= pos <= end:
            return name
    return None
