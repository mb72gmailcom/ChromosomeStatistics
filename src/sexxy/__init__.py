"""Compute genotype statistics from VCF files for male and female children."""

from sexxy.chrx import CHRX_REGIONS, chrx_region, is_chrx
from sexxy.gnomad import GnomadAfStore, DEFAULT_GNOMAD_AF_DIR, gnomad_af_path, load_gnomad_af_json
from sexxy.metadata import load_children_by_sex
from sexxy.results import GenotypeCountResult, resolve_output_target, write_genotype_count_results
from sexxy.vcf import chrom_matches, compute_genotype_counts, is_snv

__all__ = [
    "CHRX_REGIONS",
    "DEFAULT_GNOMAD_AF_DIR",
    "GenotypeCountResult",
    "GnomadAfStore",
    "chrx_region",
    "chrom_matches",
    "compute_genotype_counts",
    "gnomad_af_path",
    "is_chrx",
    "is_snv",
    "load_children_by_sex",
    "load_gnomad_af_json",
    "resolve_output_target",
    "write_genotype_count_results",
]
__version__ = "0.1.0"
