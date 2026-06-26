"""Command-line interface for sexxy."""

from __future__ import annotations

import argparse
import sys

from sexxy.chrx import is_chrx
from sexxy.gnomad import DEFAULT_GNOMAD_AF_DIR
from sexxy.metadata import load_children_by_sex
from sexxy.results import resolve_output_target, result_to_json, write_genotype_count_results
from sexxy.table import read_table
from sexxy.vcf import compute_genotype_counts


def _load_allele_freqs(path: str, key_col: str) -> dict[str, float]:
    fieldnames, rows = read_table(path)
    if key_col not in fieldnames:
        raise ValueError(f"allele frequency file missing column: {key_col!r}")
    af_col = "af" if "af" in fieldnames else fieldnames[-1]
    out: dict[str, float] = {}
    for row in rows:
        key = row[key_col]
        raw = row.get(af_col, "")
        try:
            out[key] = float(raw) if raw != "" else 0.0
        except ValueError:
            out[key] = 0.0
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute SNV genotype counts from a VCF for male/female children.",
    )
    parser.add_argument("vcf", help="Chromosome-specific VCF (.vcf or .vcf.gz)")
    parser.add_argument("metadata", help="Sample metadata (CSV/TSV)")
    parser.add_argument(
        "-c",
        "--chromosome",
        required=True,
        help="Chromosome name (e.g. chr1, chrX); must match VCF and gnomAD paths",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output filename or prefix (within --output-dir when set). "
            "Autosomes: one JSON file. chrX: six region files."
        ),
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        help="Directory for output file(s); created if missing",
    )
    parser.add_argument("--patient-col", default="patient_id")
    parser.add_argument("--father-col", default="father_id")
    parser.add_argument("--mother-col", default="mother_id")
    parser.add_argument("--sex-col", default="sex")
    parser.add_argument("--metadata-sep", default=None, help="Metadata delimiter")
    parser.add_argument(
        "--allele-freqs",
        help="Optional CSV/TSV with variant IDs and allele frequencies",
    )
    parser.add_argument(
        "--gnomad-af-dir",
        default=None,
        help=(
            "gnomAD v4 base directory with per-chromosome JSON files "
            f"(default when set: {DEFAULT_GNOMAD_AF_DIR})"
        ),
    )
    parser.add_argument("--af-key-col", default="id", help="Key column in AF file")
    parser.add_argument(
        "--common-freq-cutoff",
        type=float,
        default=0.01,
        help="Skip variants with AF above this cutoff (default: 0.01)",
    )
    parser.add_argument(
        "--min-gq",
        type=float,
        default=None,
        help="Skip genotype calls with GQ below this value",
    )
    parser.add_argument(
        "--min-dp",
        type=int,
        default=None,
        help="Skip genotype calls with DP below this value",
    )
    parser.add_argument(
        "--ab-threshold",
        type=float,
        default=None,
        help="For 0/1 and 1/1 calls, require AB > threshold",
    )
    parser.add_argument(
        "--min-gq-nonpar",
        type=float,
        default=None,
        help="chrX noPar region GQ cutoff (default: --min-gq)",
    )
    parser.add_argument(
        "--min-dp-nonpar",
        type=int,
        default=None,
        help="chrX noPar region DP cutoff (default: --min-dp)",
    )
    parser.add_argument(
        "--ab-threshold-nonpar",
        type=float,
        default=None,
        help="chrX noPar region AB cutoff (default: --ab-threshold)",
    )
    args = parser.parse_args(argv)

    male_children, female_children, _ = load_children_by_sex(
        args.metadata,
        patient_col=args.patient_col,
        father_col=args.father_col,
        mother_col=args.mother_col,
        sex_col=args.sex_col,
        sep=args.metadata_sep,
    )

    allele_freqs = None
    gnomad_af = None
    if args.allele_freqs and args.gnomad_af_dir:
        parser.error("use only one of --allele-freqs or --gnomad-af-dir")
    if args.allele_freqs:
        allele_freqs = _load_allele_freqs(args.allele_freqs, args.af_key_col)
    elif args.gnomad_af_dir is not None:
        gnomad_af = args.gnomad_af_dir

    result = compute_genotype_counts(
        args.vcf,
        male_children,
        female_children,
        chromosome=args.chromosome,
        allele_freqs=allele_freqs,
        gnomad_af=gnomad_af,
        common_freq_cutoff=args.common_freq_cutoff,
        af_key_col=args.af_key_col,
        min_gq=args.min_gq,
        min_dp=args.min_dp,
        ab_threshold=args.ab_threshold,
        min_gq_nonpar=args.min_gq_nonpar,
        min_dp_nonpar=args.min_dp_nonpar,
        ab_threshold_nonpar=args.ab_threshold_nonpar,
    )

    n_male = len(male_children)
    n_female = len(female_children)
    output_target = resolve_output_target(
        args.output, args.output_dir, args.chromosome
    )

    if is_chrx(args.chromosome):
        paths = write_genotype_count_results(
            result,
            output_target,
            male_children=n_male,
            female_children=n_female,
        )
        for path in paths:
            print(path, file=sys.stderr)
        return 0

    if output_target is not None:
        write_genotype_count_results(
            result,
            output_target,
            male_children=n_male,
            female_children=n_female,
        )
    else:
        print(
            result_to_json(
                result,
                male_children=n_male,
                female_children=n_female,
            ),
            end="",
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
