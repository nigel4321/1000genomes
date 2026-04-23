# synthetic_people

Generate synthetic single-person VCFs: each file represents one person with **one clinically-highlighted variant** drawn from ClinVar, on top of a **realistic background of common variants** sampled from local 1000 Genomes VCFs with genotypes drawn under Hardy-Weinberg equilibrium.

All output is VCFv4.2, bgzipped, and tabix-indexed.

## Requirements

- Python 3.8+ (stdlib only — no pip dependencies)
- `bcftools`, `tabix`, `bgzip` (htslib) on `PATH`
- Network access on first run only (to download ClinVar into `cache/`)
- Optional: local 1000 Genomes Phase 3 VCFs alongside this directory, for the common-variant background. Without them the output is highlighted-variant-only.

## Quick start

```bash
# 10 person VCFs with reproducible output
python3 synthetic_people/generate_people.py --n 10 --seed 42

# 100 non-deterministic people — each run produces different output
python3 synthetic_people/generate_people.py --n 100

# Fresh output dir + only benign ClinVar variants + more background per person
python3 synthetic_people/generate_people.py --n 5 \
    --output-dir /tmp/synth_people \
    --clinvar-sig Benign,Likely_benign \
    --n-background 1500
```

Output lands in `synthetic_people/out/person_0001.vcf.gz`, `person_0002.vcf.gz`, etc.

## Data sources

### ClinVar (highlighted variants)

Downloaded from NCBI once and cached in `./cache/`:

- GRCh37: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz` (~85 MB)
- GRCh38: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`

The cache persists across runs. To force a refresh, delete `./cache/clinvar_<build>.vcf.gz*` and re-run.

Candidates are filtered on the `CLNSIG` (clinical significance) INFO field. Default: `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`. Override with `--clinvar-sig`.

### 1000 Genomes (background)

By default the script globs `../ALL.chr*.phase3_*.genotypes.vcf.gz` relative to itself — i.e. the 1000G VCFs in the parent `1000genomes/` directory. Override with `--background-glob` (repeat the flag to combine multiple sources).

For each source, 5,000 variants with `INFO/AF >= 0.05` (default) are reservoir-sampled into a global pool. For each person, `--n-background` (default 500) are drawn from the pool, then their genotypes are rolled under Hardy-Weinberg equilibrium from each site's allele frequency. Hom-ref draws are discarded (no `0|0` lines in the output), so a person typically ends up with ~200–300 emitted background records.

## Randomness

- `--seed N` — deterministic: same inputs + same seed → byte-identical output (same sample IDs, same variants picked, same genotypes drawn).
- Omit `--seed` — each invocation produces different people: different sample IDs (HG/NA-prefixed, 1000G-style), different highlighted variants, different background genotypes.

Every output is still VCF-compliant, internally consistent (AC/AN/AF match the emitted genotype), and uses human chromosome names + a declared GRCh37/GRCh38 reference.

## Output layout

```
synthetic_people/
├── cache/                         # ClinVar download (persists across runs)
│   ├── clinvar_GRCh37.vcf.gz
│   └── clinvar_GRCh37.vcf.gz.tbi
└── out/
    ├── person_0001.vcf.gz
    ├── person_0001.vcf.gz.tbi
    ├── person_0002.vcf.gz
    └── ...
```

Each `person_<N>.vcf.gz` has:

- Header: `##fileformat=VCFv4.2`, `##reference=...`, full GRCh37/GRCh38 `##contig` lines, `INFO` declarations for `AC`/`AN`/`AF` + `HIGHLIGHT`/`CLNSIG`/`CLNDN`, `FORMAT=GT`.
- One sample column with a random HG/NA-prefixed ID.
- One row per emitted variant, sorted by chromosome order (matching the reference build) and position.
- The highlighted row carries `HIGHLIGHT` flag plus `CLNSIG` and `CLNDN` values from ClinVar.

## Compatibility with the nextflow pipeline

Output files pass `bin/qc_validate.py` cleanly (GRCh37/GRCh38 reference recognised, human contigs, `GT` declared, `AF`/`AC`+`AN` present). You can feed them straight into the pipeline:

```bash
nextflow run nextflow_pipeline/main.nf \
    --input 'synthetic_people/out/*.vcf.gz' \
    --outdir /tmp/synth_scan_results \
    -params-file nextflow_pipeline/variants/rs12913832.yaml
```

Because each file has one sample, cohort-level AF and carrier counts in the variant/carriers reports will be per-person rather than population-level.

## CLI reference

| Flag | Purpose | Default |
|---|---|---|
| `--n` | Number of person VCFs to generate | `10` |
| `--seed` | RNG seed; omit for fresh randomness each run | `None` |
| `--build` | `GRCh37` or `GRCh38` (must match background source) | `GRCh37` |
| `--output-dir` | Where to write `person_<N>.vcf.gz` | `./out` |
| `--cache-dir` | Where ClinVar is downloaded | `./cache` |
| `--background-glob` | Source glob(s) for common variants (repeatable) | 1000G files in parent dir |
| `--n-background` | Max background variants sampled per person (pre-filter) | `500` |
| `--af-min` | Minimum AF for background variants | `0.05` |
| `--clinvar-sig` | Comma-separated CLNSIG values to include | `Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic` |
