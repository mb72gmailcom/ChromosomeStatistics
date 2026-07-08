"""Command-line interface for sexxy."""

from __future__ import annotations

import argparse
import sys

from sexxy import __version__
from sexxy.gnomad import DEFAULT_GNOMAD_AF_DIR
from sexxy.metadata import load_children_by_sex
from sexxy.results import (
    resolve_output_target,
    write_genotype_count_results,
    write_run_params,
    write_skipped_contigs,
)
from sexxy.table import read_table
from sexxy.vcf import _chrom_key, compute_genotype_counts


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
            "Autosomes/chrY: two files ({prefix}.male.json, .female.json). "
            "chrX: six files (male/female x Par1/noPar/Par2)."
        ),
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        help="Directory for output file(s); created if missing",
    )
    parser.add_argument("--patient-col", default="spid")
    parser.add_argument("--father-col", default="father")
    parser.add_argument("--mother-col", default="mother")
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
        help="chrX noPar GQ cutoff for males only (default: --min-gq)",
    )
    parser.add_argument(
        "--min-dp-nonpar",
        type=int,
        default=None,
        help="chrX noPar DP cutoff for males only (default: --min-dp)",
    )
    parser.add_argument(
        "--ab-threshold-nonpar",
        type=float,
        default=None,
        help="chrX noPar AB cutoff for males only (default: --ab-threshold)",
    )
    parser.add_argument(
        "--exclude-repeats",
        default=None,
        help=(
            "Tab-separated file of chrom/start/end intervals; skip variants "
            "with start <= POS < end (semi-open)"
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Require every metadata child to appear in the VCF header; "
            "default is to exclude children not present in the VCF"
        ),
    )
    parser.add_argument(
        "--check-father",
        action="store_true",
        help=(
            "chrY only: for male children, count how often 0/0, 0/1, and 1/1 "
            "match the father's genotype (keys 0/0_f, 0/1_f, 1/1_f)"
        ),
    )
    args = parser.parse_args(argv)

    male_children, female_children, children_rows = load_children_by_sex(
        args.metadata,
        patient_col=args.patient_col,
        father_col=args.father_col,
        mother_col=args.mother_col,
        sex_col=args.sex_col,
        sep=args.metadata_sep,
    )

    male_father_by_child = None
    if args.check_father:
        if _chrom_key(args.chromosome).upper() != "Y":
            print("--check-father applies only to chrY; ignoring", file=sys.stderr)
        else:
            male_father_by_child = {
                row[args.patient_col]: row[args.father_col]
                for row in children_rows
                if row[args.patient_col] in male_children
            }

    allele_freqs = None
    gnomad_af = None
    if args.allele_freqs and args.gnomad_af_dir:
        parser.error("use only one of --allele-freqs or --gnomad-af-dir")
    if args.allele_freqs:
        allele_freqs = _load_allele_freqs(args.allele_freqs, args.af_key_col)
    elif args.gnomad_af_dir is not None:
        gnomad_af = args.gnomad_af_dir

    def _report_excluded(excluded_male: list[str], excluded_female: list[str]) -> None:
        nm, nf = len(excluded_male), len(excluded_female)
        print(
            f"Excluded {nm} male and {nf} female children not in VCF header",
            file=sys.stderr,
        )
        examples = excluded_male[:3] + excluded_female[:3]
        if examples:
            print(f"  e.g. {examples}", file=sys.stderr)

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
        exclude_repeats=args.exclude_repeats,
        strict=args.strict,
        on_excluded=_report_excluded,
        check_father=args.check_father,
        male_father_by_child=male_father_by_child,
    )

    n_male = result.male_cohort_size
    n_female = result.female_cohort_size
    if n_male or n_female:
        print(f"Cohort: {n_male} male, {n_female} female", file=sys.stderr)
    if result.excluded_repeat_rows:
        print(
            f"Excluded {result.excluded_repeat_rows} variant rows in "
            f"{result.exclude_repeat_intervals} repeat interval(s)",
            file=sys.stderr,
        )
    output_target = resolve_output_target(
        args.output, args.output_dir, args.chromosome
    )

    paths = write_genotype_count_results(
        result,
        output_target,
        male_children=n_male,
        female_children=n_female,
    )

    skipped_path = write_skipped_contigs(
        output_target,
        args.chromosome,
        vcf=args.vcf,
        target_chrom_key=_chrom_key(args.chromosome),
        skipped_contigs=result.skipped_contigs,
    )
    paths.append(skipped_path)
    total_skipped = sum(result.skipped_contigs.values())
    if total_skipped:
        n_contigs = len(result.skipped_contigs)
        print(
            f"Skipped {total_skipped} variant rows on {n_contigs} non-matching contig(s)",
            file=sys.stderr,
        )

    params = {
        "version": __version__,
        "command": sys.argv,
        "inputs": {
            "vcf": args.vcf,
            "metadata": args.metadata,
        },
        "chromosome": args.chromosome,
        "metadata_columns": {
            "patient_col": args.patient_col,
            "father_col": args.father_col,
            "mother_col": args.mother_col,
            "sex_col": args.sex_col,
            "metadata_sep": args.metadata_sep,
        },
        "allele_frequency": {
            "gnomad_af_dir": args.gnomad_af_dir,
            "allele_freqs": args.allele_freqs,
            "af_key_col": args.af_key_col,
            "common_freq_cutoff": args.common_freq_cutoff,
        },
        "filters": {
            "min_gq": args.min_gq,
            "min_dp": args.min_dp,
            "ab_threshold": args.ab_threshold,
            "min_gq_nonpar": args.min_gq_nonpar,
            "min_dp_nonpar": args.min_dp_nonpar,
            "ab_threshold_nonpar": args.ab_threshold_nonpar,
            "strict": args.strict,
            "check_father": args.check_father,
            "exclude_repeats": args.exclude_repeats,
        },
        "cohort": {
            "metadata_male": len(male_children),
            "metadata_female": len(female_children),
            "male_children": n_male,
            "female_children": n_female,
            "excluded_male": list(result.excluded_male),
            "excluded_female": list(result.excluded_female),
        },
        "skipped_contigs_file": str(skipped_path),
        "total_skipped_rows": total_skipped,
        "exclude_repeat_intervals": result.exclude_repeat_intervals,
        "excluded_repeat_rows": result.excluded_repeat_rows,
        "output_files": [str(p) for p in paths],
    }
    params_path = write_run_params(output_target, args.chromosome, params)
    paths.append(params_path)

    for path in paths:
        print(path, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
