# variant-scan — Nextflow pipeline

Scan a directory of VCF files for a target variant, producing two markdown reports:

1. **`metadata_report.md`** — cohort-level and per-file overview for researchers orienting themselves to an unfamiliar dataset.
2. **`variant_report.md`** — which files contain the target variant within an acceptable allele-frequency range, with per-population breakdowns.

Designed for a standalone VM (local executor). Scales to hundreds of VCFs via per-file parallelism.

## Requirements

- Nextflow ≥ 22.10 — `curl -s https://get.nextflow.io | bash`
- `bcftools` ≥ 1.9 on PATH
- Python 3.8+ (stdlib only — no pip install required)
- Each input VCF must be bgzipped (`.vcf.gz`) with a matching tabix index (`.vcf.gz.tbi`)

## Usage

```bash
nextflow run main.nf \
    --input '/data/vcfs/*.vcf.gz' \
    --outdir results \
    -params-file variants/rs12913832.yaml
```

All variant fields can also be passed as CLI flags (`--variant_name`, `--variant_chrom`, etc.) instead of a params file.

### Parameters

| Parameter | Required | Description |
|---|---|---|
| `--input` | yes | Glob for input VCFs (quote it) |
| `--variant_name` | yes | Label used in report headings |
| `--variant_chrom` | yes | Chromosome; `chr`-prefix optional, pipeline resolves |
| `--variant_pos` | yes | 1-based position |
| `--variant_ref` | yes | Reference allele (strict match) |
| `--variant_alt` | yes | Alternate allele (strict match) |
| `--variant_min_af` | no | Lower AF bound (default 0.0) |
| `--variant_max_af` | no | Upper AF bound (default 1.0) |
| `--outdir` | no | Output directory (default `results`) |

### Adding new variants

Copy `variants/rs12913832.yaml` and edit. One file per variant keeps the parameter history auditable.

## Why match by position + allele rather than rsID

rsIDs are not guaranteed to be present. The 1000 Genomes Phase 3 release, for example, was frozen against a 2013 dbSNP snapshot — many records carry `.` in the ID column even for variants that today have a well-known rsID. Matching by coordinate + REF + ALT is portable across releases and correctly excludes overlapping structural variants (e.g. `<CN0>`, `<CN2>`) at the same position.

## Pipeline layout

```
nextflow_pipeline/
├── main.nf                          # entry workflow
├── nextflow.config                  # executor, resources, params
├── modules/
│   ├── inspect_vcf.nf               # per-VCF metadata collection
│   ├── scan_variant.nf              # per-VCF variant lookup
│   └── reports.nf                   # markdown aggregation
├── bin/                             # scripts auto-added to PATH by Nextflow
│   ├── inspect_vcf.py
│   ├── scan_variant.py
│   ├── build_metadata_report.py
│   └── build_variant_report.py
└── variants/
    └── rs12913832.yaml              # example variant spec
```

## Execution model

Per input VCF, Nextflow runs two independent parallel tasks:

- `INSPECT_VCF` — contigs, sample count, variant count, reference build heuristic, pipeline/date tags, sample-list hash (to detect cohort mismatches across files).
- `SCAN_VARIANT` — tabix region query, strict REF/ALT match, AF extraction from INFO (or recomputed via `bcftools +fill-tags` if missing), classification into one of seven statuses.

Aggregation runs once both per-file fan-outs finish:

- `METADATA_REPORT` → `metadata_report.md`
- `VARIANT_REPORT` → `variant_report.md` (highlights files where the variant is in range)

## Variant statuses

Each file's scan emits one of:

| Status | Meaning |
|---|---|
| `not_applicable` | Target chromosome not in this VCF |
| `position_empty` | Chromosome present, but no variant call at the position |
| `absent` | Variant(s) at position but none with the specified REF/ALT |
| `present_in_range` | Variant found and AF ∈ [min_af, max_af] |
| `present_below_threshold` | Variant found but AF < min_af |
| `present_above_threshold` | Variant found but AF > max_af |
| `present_af_unknown` | Variant found but AF could not be determined |

## Profiles

Tune parallelism via `-profile`:

- `standard` — default, `queueSize=16`
- `small_vm` — `queueSize=4`
- `big_vm` — `queueSize=32`

## Output directory

```
results/
├── metadata_report.md
└── variant_report.md
```

Per-file intermediate JSONs stay in Nextflow's `work/` directory and can be inspected for debugging.
