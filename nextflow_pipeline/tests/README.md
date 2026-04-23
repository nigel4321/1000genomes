# variant-scan tests

Tests are written against Python's stdlib `unittest` — no pip dependencies.

## Layout

| File | Scope |
|---|---|
| `synthetic_vcf.py` | VCF 4.1-compliant synthetic data generator (1000G Phase 3 conventions: `##source=1000GenomesPhase3Pipeline`, `hs37d5` reference, per-super-pop AF INFO fields, phased diploid genotypes, `ALL.chrN.phase3_shapeit2_mvncall_integrated_v5b.<date>.genotypes.vcf.gz` filenames). |
| `fixtures.py` | Shared cohort + variant fixtures. The default cohort is 20 real 1000G sample IDs (4 per super-pop: EUR/AFR/EAS/SAS/AMR). |
| `test_inspect_vcf.py` | `inspect_vcf.py` metadata extraction — contigs, sample count, sample-hash stability, reference, pipeline/date tags. |
| `test_scan_variant.py` | `scan_variant.py` — one test per pipeline status (`not_applicable`, `position_empty`, `absent`, `present_in_range`, `present_below_threshold`, `present_above_threshold`, `present_af_unknown`) plus carrier extraction and `chr`-prefix resolution. |
| `test_reports.py` | The three markdown builders (`build_metadata_report.py`, `build_variant_report.py`, `build_carrier_report.py`) — aggregation, cohort-mismatch flagging, het+2·hom integrity check. |
| `test_pipeline_e2e.py` | End-to-end `nextflow run main.nf` against two synthetic VCFs, asserts the four published outputs look right. |

## Requirements

- Python 3.8+ (stdlib only)
- `bcftools`, `tabix`, `bgzip` (htslib) on `PATH` — the synthetic data generator uses them to bgzip and index.
- `nextflow` on `PATH` (only needed for `test_pipeline_e2e.py`; the test auto-skips if missing).

## Running

```bash
# from the pipeline root
python3 -m unittest discover -s tests -v

# or a specific module
python3 -m unittest tests.test_scan_variant -v

# or a single test
python3 -m unittest tests.test_scan_variant.ScanVariantTest.test_status_present_in_range_and_carrier_counts -v
```

The e2e test takes ~15 s (Nextflow startup). Set `NF_KEEP=1` to preserve the Nextflow `work/` dir from that test for debugging.

## Synthetic data overview

The generator produces bgzipped + tabix-indexed VCFs that share all the structural features the pipeline reads:

- VCFv4.1 header with `##reference`, `##source=1000GenomesPhase3Pipeline`, and all GRCh37 `##contig` lines.
- `INFO` fields: `AC`, `AF`, `AN`, `NS`, `EAS_AF`, `EUR_AF`, `AFR_AF`, `AMR_AF`, `SAS_AF`, `VT`.
- `GT`-only `FORMAT`, phased genotypes (`a|b`).
- Per-super-pop allele frequencies are *computed* from the drawn genotypes (or directly from supplied genotypes), so carrier counts and INFO AF are always internally consistent.
- Filename matches the real 1000G convention so `inspect_vcf.py`'s regex-based pipeline-tag and date-stamp extraction is exercised.

See `synthetic_vcf.Variant` for ways to construct a record (direct genotype list, AF-per-super-pop draw, or edge cases like `strip_af_info=True` / `all_missing=True`).

## Generating synthetic VCFs from the command line

`synthetic_vcf.py` is also runnable as a CLI for quick one-off data generation outside the test suite — useful when you want real input to explore the pipeline with, without downloading the full 1000G chromosome files.

```bash
python3 tests/synthetic_vcf.py --help
```

### Minimal example

```bash
python3 tests/synthetic_vcf.py --chrom 22 --outdir /tmp/synth --n-variants 100
```

Produces `/tmp/synth/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` (plus `.tbi`) with 100 random SNPs drawn from a skewed-low MAF distribution and the default 20-sample 1000G-style cohort.

### Including a known target variant

Splice in one specific variant at a fixed position and per-super-pop AF — e.g. rs12913832:

```bash
python3 tests/synthetic_vcf.py \
    --chrom 15 --outdir /tmp/synth --n-variants 200 \
    --target-name rs12913832 \
        --target-pos 28365618 --target-ref A --target-alt G \
        --target-af-eur 0.75 --target-af-afr 0.02 --target-af-eas 0.05 \
        --target-af-sas 0.10 --target-af-amr 0.30 \
    --seed 42
```

### Feeding the generated data into the pipeline

```bash
python3 tests/synthetic_vcf.py --chrom 15 --outdir /tmp/synth --n-variants 50 \
    --target-name rs12913832 --target-pos 28365618 --target-ref A --target-alt G \
    --target-af-eur 0.75

nextflow run main.nf \
    --input '/tmp/synth/*.vcf.gz' \
    -params-file variants/rs12913832.yaml
```

### Key flags

| Flag | Purpose |
|---|---|
| `--chrom` (required) | GRCh37 contig name, no `chr` prefix (`1`–`22`, `X`, `Y`). |
| `--output` / `--outdir` | Either override the full output path, or just the directory (the 1000G-style filename is generated). |
| `--n-variants` | How many random SNPs to emit (default 100). |
| `--pos-start` / `--pos-stop` | Bound the position range. Defaults to 1 Mb → chromosome length. |
| `--seed` | Deterministic — drives both positions and genotypes. |
| `--target-name`, `--target-pos`, `--target-ref`, `--target-alt` | Splice in one named variant at fixed coordinates. |
| `--target-af-{eur,afr,eas,sas,amr}` | Per-super-pop AF for the target (used to draw realistic carrier genotypes). |
| `--date-stamp` | Date component of the default filename (default `20130502`). |

For more complex scenarios (multiple named variants, custom cohorts, stripped-AF or all-missing edge cases) import the module — the CLI deliberately covers just the common "generate one plausible per-chromosome VCF" case.
