"""Semi-open genomic intervals for excluding repeat regions from counting."""

from __future__ import annotations

from pathlib import Path


def _chrom_key(chrom: str) -> str:
    c = chrom.strip()
    if c.lower().startswith("chr"):
        return c[3:]
    return c


class ExcludeIntervals:
    """Sorted semi-open intervals ``[start, end)`` for fast position lookup.

    VCF rows are scanned in increasing ``POS`` order, so a moving index is used.
    """

    def __init__(self, intervals: list[tuple[int, int]]):
        self.intervals = intervals
        self._i = 0

    def __len__(self) -> int:
        return len(self.intervals)

    def contains(self, pos: int) -> bool:
        """Return whether *pos* lies in any interval (``start <= pos < end``)."""
        while self._i < len(self.intervals) and self.intervals[self._i][1] <= pos:
            self._i += 1
        j = self._i
        while j < len(self.intervals) and self.intervals[j][0] <= pos:
            start, end = self.intervals[j]
            if start <= pos < end:
                return True
            j += 1
        return False


def load_exclude_intervals(
    path: str | Path,
    chromosome: str,
) -> ExcludeIntervals:
    """Load tab-separated ``chrom start end`` intervals for *chromosome*.

    Each line defines a semi-open interval ``[start, end)``. Lines whose
    chromosome does not match *chromosome* (``chr1`` / ``1`` equivalent) are
    ignored. Intervals must be sorted by ``start``.
    """
    target = _chrom_key(chromosome)
    intervals: list[tuple[int, int]] = []

    with open(path) as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"{path}:{line_no}: expected chrom, start, end tab-separated fields"
                )
            if _chrom_key(parts[0]) != target:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_no}: interval start/end must be integers"
                ) from exc
            if start >= end:
                raise ValueError(
                    f"{path}:{line_no}: invalid interval [{start}, {end}); require start < end"
                )
            intervals.append((start, end))

    if not intervals:
        raise ValueError(
            f"no intervals for chromosome {chromosome!r} found in {path}"
        )

    for i in range(1, len(intervals)):
        if intervals[i][0] < intervals[i - 1][0]:
            raise ValueError(
                f"intervals in {path} must be sorted by start position"
            )

    return ExcludeIntervals(intervals)
