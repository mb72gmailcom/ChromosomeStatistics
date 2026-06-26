# sexXY

Compute SNV genotype counts from VCF files, separately for male and female children.

## Install

```bash
pip install -e ".[dev]"
```

## Metadata

Tab- or comma-separated file with at least:

| Column       | Description                          |
|-------------|--------------------------------------|
| `patient_id` | Sample / patient identifier          |
| `father_id`  | Father ID (empty if not a child)     |
| `mother_id`  | Mother ID (empty if not a child)     |
| `sex`        | `Male`/`Female`, `1`/`2`, `M`/`F`, etc. |

Children are rows with **both** a valid father and mother ID. Two ID lists are built: `male_children` and `female_children`.

Column names are configurable via function arguments or CLI flags.

## Usage

All inputs are **chromosome-specific**. Pass the chromosome name with
``chromosome=`` (API) or ``--chromosome`` / ``-c`` (CLI). ``chr1`` and ``1``
are treated as the same chromosome when matching VCF rows.

### Python API

```python
from sexxy import load_children_by_sex, compute_genotype_counts

chrm = "chr1"
male_children, female_children, children_df = load_children_by_sex(
    "metadata.tsv",
    sep="\t",
)

result = compute_genotype_counts(
    f"cohort.{chrm}.vcf.gz",
    male_children,
    female_children,
    chromosome=chrm,
    gnomad_af="/mnt/home/mbershadsky/ceph/gnomad.v4",
    common_freq_cutoff=0.01,
)

print(result.male_counts())   # e.g. {'0/0': 12345, '0/1': 678, '1/1': 90}
print(result.female_counts())
```

The VCF is scanned **once**; male and female genotype counts are accumulated in the same pass.

### gnomAD v4 allele frequencies

For a given chromosome, the AF file is loaded once:

```
{gnomad_af_dir}/{chrm}/{chrm}-common-af.json
```

Example (matches the reference workflow):

```python
chrm = "chr1"
gfile = f"/mnt/home/mbershadsky/ceph/gnomad.v4/{chrm}/{chrm}-common-af.json"
with open(gfile) as f:
    df = json.load(f)
# variants with df.get(variant_id, 0) > 0.01 are skipped
```

```python
from sexxy import GnomadAfStore, compute_genotype_counts

chrm = "chr1"
store = GnomadAfStore("/mnt/home/mbershadsky/ceph/gnomad.v4")
result = compute_genotype_counts(
    vcf, male_children, female_children, chromosome=chrm, gnomad_af=store,
)
print(result.male_counts())
```

### chrX regions

For chrX, counts are split into three regions (sites outside these intervals are
skipped):

| Region | Start | End |
|--------|------:|----:|
| `par1` | 10,001 | 2,781,479 |
| `noPar` | 2,781,489 | 155,701,382 |
| `par2` | 155,701,383 | 156,030,895 |

```python
from sexxy import compute_genotype_counts, write_genotype_count_results

result = compute_genotype_counts(vcf, male_children, female_children, chromosome="chrX")
print(result.male_counts("par1"))
print(result.female_counts("noPar"))

write_genotype_count_results(result, "counts.chrX", male_children=len(male_children), female_children=len(female_children))
# writes counts.chrX.male.par1.json, counts.chrX.male.noPar.json, counts.chrX.male.par2.json,
#         counts.chrX.female.par1.json, counts.chrX.female.noPar.json, counts.chrX.female.par2.json
```

Use `--min-gq-nonpar`, `--min-dp-nonpar`, and `--ab-threshold-nonpar` to apply different
quality cutoffs in the `noPar` region (defaults match the global flags).

```bash
sexxy cohort.chrX.vcf.gz metadata.tsv --chromosome chrX -o counts.chrX
```

### CLI

From an installed package:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 -o counts.chr1.json
```

From a checkout (no install required):

```bash
python run.py cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 -o counts.chr1.json
```

Write outputs to a directory with ``--output-dir`` / ``-d``. The directory is
created automatically if it does not exist.

```bash
python run.py cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 --output-dir results/
# -> results/counts.chr1.json

python run.py cohort.chrX.vcf.gz metadata.tsv --chromosome chrX --output-dir results/
# -> results/counts.chrX.male.par1.json, ... (six files)
```

Combine with ``-o`` to set the filename/prefix inside the directory:

```bash
python run.py cohort.chrX.vcf.gz metadata.tsv --chromosome chrX \
  --output-dir results/ -o mycohort.chrX
```

Filter common variants with gnomAD v4 JSON files:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv \
  --chromosome chr1 \
  --gnomad-af-dir /mnt/home/mbershadsky/ceph/gnomad.v4 \
  -o counts.chr1.json
```

Or a flat chromosome-specific allele-frequency file:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv \
  --chromosome chr1 \
  --allele-freqs chr1-common-af.tsv \
  --af-key-col id \
  -o counts.chr1.json
```

## Behavior

- **SNVs only**: rows where `len(REF) == 1` and `len(ALT) == 1`
- **Genotype field**: first sub-field of `FORMAT` (e.g. `0/1` from `0/1:25:15,10`)
- **Missing samples**: raises if a child ID is absent from the VCF header
- Supports plain `.vcf` and `.vcf.gz`

### Optional genotype quality filters

Set any of these parameters to enable per-call filtering (unset = no filter):

| Parameter | VCF field | Rule |
|-----------|-----------|------|
| `min_gq` | `GQ` | skip call if `GQ < min_gq` |
| `min_dp` | `DP` | skip call if `DP < min_dp` |
| `ab_threshold` | `AB` or `AD` | for `0/1` and `1/1`: require `AB > ab_threshold` |

`0/0` calls are not subject to the AB filter. When `AB` is absent, it is computed as `alt / (ref + alt)` from `AD`.

For chrX `noPar`, you can override the cutoffs with `min_gq_nonpar`, `min_dp_nonpar`, and
`ab_threshold_nonpar` (each defaults to the global value when unset).

```python
result = compute_genotype_counts(
    vcf, male_children, female_children,
    chromosome="chrX",
    min_gq=20,
    min_dp=10,
    ab_threshold=0.2,
    min_dp_nonpar=5,
)
```

CLI:

```bash
sexxy cohort.chrX.vcf.gz metadata.tsv --chromosome chrX \
  --min-gq 20 --min-dp 10 --ab-threshold 0.2 \
  --min-dp-nonpar 5 -o counts.chrX
```

Autosomes:

```python
result = compute_genotype_counts(
    vcf, male_children, female_children,
    chromosome="chr1",
    min_gq=20,
    min_dp=10,
    ab_threshold=0.2,
)
```

CLI:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 \
  --min-gq 20 --min-dp 10 --ab-threshold 0.2 -o counts.chr1.json
```
