# sexXY

Compute SNV genotype counts from VCF files, separately for male and female children.

## Install

```bash
pip install -e ".[dev]"
```

## Metadata

Tab- or comma-separated file with at least:

| Column | Description |
|--------|-------------|
| `spid` | Patient / sample identifier |
| `sfid` | Family identifier (optional; not used in counting) |
| `father` | Father ID (empty if not a child) |
| `mother` | Mother ID (empty if not a child) |
| `sex` | `Male`/`Female`, `1`/`2`, `M`/`F`, etc. |

Children are rows with **both** a valid father and mother ID. Two ID lists are built: `male_children` and `female_children`.

Only children whose **sample ID appears in the VCF header** are included in counting. Children present in metadata but absent from the VCF are excluded automatically; the CLI prints how many were dropped. Use ``--strict`` to fail instead if any child is missing from the VCF.

Column names are configurable via function arguments or CLI flags
(``--patient-col``, ``--father-col``, ``--mother-col``, ``--sex-col``) if your
file uses different headers.

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

### Output files

**Autosomes and chrY** — two files (no ``region`` field):

```
counts.chr19.male.json
counts.chr19.female.json
```

Example ``counts.chr19.female.json``:

```json
{
  "chromosome": "chr19",
  "sex": "female",
  "female_children": 120,
  "gt_counts": {"0/0": 520, "0/1": 28}
}
```

Male files use ``"sex": "male"`` and ``male_children`` only.

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

For chrX, **six** output files are written (male and female for each region):

```
counts.chrX.male.Par1.json
counts.chrX.male.noPar.json
counts.chrX.male.Par2.json
counts.chrX.female.Par1.json
counts.chrX.female.noPar.json
counts.chrX.female.Par2.json
```

Example ``counts.chrX.male.Par1.json``:

```json
{
  "chromosome": "chrX",
  "region": "Par1",
  "sex": "male",
  "male_children": 115,
  "gt_counts": {"0/0": 500, "0/1": 30}
}
```

Female files have ``"sex": "female"``, ``female_children``, and female ``gt_counts`` only.

| Region | Start | End |
|--------|------:|----:|
| `Par1` | 10,001 | 2,781,479 |
| `noPar` | 2,781,489 | 155,701,382 |
| `Par2` | 155,701,383 | 156,030,895 |

```python
from sexxy import compute_genotype_counts, write_genotype_count_results

result = compute_genotype_counts(vcf, male_children, female_children, chromosome="chrX")
print(result.male_counts("Par1"))
print(result.female_counts("noPar"))

write_genotype_count_results(result, "counts.chrX", male_children=len(male_children), female_children=len(female_children))
# writes counts.chrX.male.Par1.json, ... (six files total)
```

Use `--min-gq-nonpar`, `--min-dp-nonpar`, and `--ab-threshold-nonpar` to apply different
quality cutoffs for **male** calls in the `noPar` region only (defaults match the global
flags). Female calls always use the global cutoffs in all chrX regions.

```bash
sexxy cohort.chrX.vcf.gz metadata.tsv --chromosome chrX -o counts.chrX
```

### CLI

From an installed package:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 --output-dir results/
# -> results/counts.chr1.male.json, results/counts.chr1.female.json,
#    results/counts.chr1.params.json
```

The **params file** records inputs, filters, cohort sizes (including any children
excluded because they are absent from the VCF), and paths to all output files.
Use it to reproduce or audit a run.

From a checkout (no install required):

```bash
python run.py cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 --output-dir results/
```

Write outputs to a directory with ``--output-dir`` / ``-d``. The directory is
created automatically if it does not exist.

```bash
python run.py cohort.chr1.vcf.gz metadata.tsv --chromosome chr1 --output-dir results/
# -> results/counts.chr1.male.json, results/counts.chr1.female.json

python run.py cohort.chrX.vcf.gz metadata.tsv --chromosome chrX --output-dir results/
# -> results/counts.chrX.male.Par1.json, ... (six files)
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
  --output-dir results/
```

Or a flat chromosome-specific allele-frequency file:

```bash
sexxy cohort.chr1.vcf.gz metadata.tsv \
  --chromosome chr1 \
  --allele-freqs chr1-common-af.tsv \
  --af-key-col id \
  --output-dir results/
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

For chrX `noPar` **male** calls, you can override the cutoffs with `min_gq_nonpar`,
`min_dp_nonpar`, and `ab_threshold_nonpar` (each defaults to the global value when unset).
Female calls always use the global filters in all chrX regions.

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
  --min-gq 20 --min-dp 10 --ab-threshold 0.2 --output-dir results/
```
