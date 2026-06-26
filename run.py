#!/usr/bin/env python3
"""Convenience entry point for the sexxy CLI.

Usage:
    python run.py cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 -o counts.chr1.json
    python run.py cohort.chrX.vcf.gz metadata.tsv --chromosome chrX --output-dir results/
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running from a checkout without ``pip install -e .``
_SRC = Path(__file__).resolve().parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from sexxy.cli import main

if __name__ == "__main__":
    sys.exit(main())
