"""Read delimiter-separated tables without third-party dependencies."""

from __future__ import annotations

import csv
from pathlib import Path


def read_table(path: str | Path, sep: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    """Read a headered CSV/TSV file into a list of row dicts."""
    with open(path, newline="") as f:
        if sep is None:
            sample = f.read(8192)
            f.seek(0)
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
            except csv.Error:
                delimiter = "\t"
        else:
            delimiter = sep

        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"no header row found in {path}")
        fieldnames = list(reader.fieldnames)
        rows = [
            {name: (row.get(name) or "") for name in fieldnames}
            for row in reader
        ]
    return fieldnames, rows
