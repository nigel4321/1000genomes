# synthetic_people

Generate cohorts of synthetic whole-genome VCFs with realistic variant
content: an **LD-correct coalescent-simulated background** under a
configurable demographic model, **ClinVar-pathogenic variants** highlighted
per person, and full VCF 4.2 **per-sample quality metrics** (GT / DP /
GQ / AD).

The spec is `SYHTHETIC_PROJECT.md`; incremental build plan and
per-milestone status is in `IMPLEMENTATION_PLAN.md`.

As of **M5** the default path is an msprime + stdpopsim coalescent; the
M4 1000G-pool + power-law SFS sampler is retained behind
`--legacy-background`.

---

## Install

System binaries (not pip-installable):

```
bcftools tabix bgzip      # htslib — header parse, bgzip, tabix-index
```

Python deps go in a project venv (repo root):

```bash
sudo apt install python3.12-venv          # only if ensurepip missing
python3 -m venv .venv
.venv/bin/pip install -r synthetic_people/requirements.txt
```

`requirements.txt` pins:

| Package | Used by | Purpose |
|---|---|---|
| `numpy>=1.24` | M4+ | sampling, histogram arithmetic |
| `msprime>=1.3`, `tskit>=0.5` | M5 | coalescent simulation + tree sequences |
| `stdpopsim>=0.2` | M5 | human demographic catalogue (`OutOfAfrica_3G09`, etc.) |
| `matplotlib>=3.7`, `scikit-allel>=1.3` | M10 | validation plots (reserved) |

Probe the environment:

```bash
.venv/bin/python synthetic_people/generate_people.py --check-deps
```

First run only: ClinVar is downloaded (~50 MB) into
`synthetic_people/cache/` and re-used thereafter.

---

## Usage

### Default: coalescent-simulated cohort (M5)

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 42 \
    --chromosomes 22 --chr-length-mb 5.0 \
    --demo-model OutOfAfrica_3G09 --population CEU
```

- `--n` — cohort size (number of person VCFs).
- `--chromosomes` — comma-separated list (e.g. `19,20,21,22`).
- `--chr-length-mb` — simulated prefix per chromosome, 0 = full length.
- `--demo-model` — stdpopsim model id; `none` falls back to constant-size
  Ne = 10 000 msprime draw.
- `--population` — sampling population for the demo model (CEU, YRI, CHB
  for `OutOfAfrica_3G09`).
- `--rec-rate`, `--mu` — only consulted when `--demo-model=none`.

### Legacy: 1000G-pool + power-law SFS (M4)

```bash
.venv/bin/python synthetic_people/generate_people.py \
    --n 50 --seed 42 --legacy-background --build GRCh37 \
    --background-glob "ALL.chr22.phase3_*.vcf.gz" \
    --n-background 500 --sfs-alpha 2.0
```

Retained for comparison, offline-only use, and any workflow that wants
to draw coordinates directly from the local 1000G Phase 3 data.
Legacy-only flags (`--background-glob`, `--n-background`, `--af-min`,
`--sfs-alpha`) are marked `[legacy]` in `--help`.

### Reproducibility

- `--seed N` — same inputs + same seed → byte-identical output.
- Omit `--seed` — each invocation produces different people (different
  sample IDs, different highlighted variants, different genotypes).

---

## Data sources

### ClinVar (highlighted variants)

Downloaded from NCBI once, cached in `./cache/`:

- GRCh37: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`
- GRCh38: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`

Filtered on `CLNSIG`. Default `Pathogenic`, `Likely_pathogenic`,
`Pathogenic/Likely_pathogenic`; override via `--clinvar-sig`. Delete
`cache/clinvar_<build>.vcf.gz*` to force a refresh.

### 1000 Genomes (legacy background only)

Default glob: `../ALL.chr*.phase3_*.genotypes.vcf.gz`. For each source,
5 000 variants with `MAX(INFO/AF) >= --af-min` (default 0.05) are
reservoir-sampled. Only coordinates + allele strings carry through —
source AFs are ignored downstream because M4's SFS sampler redraws
them per site.

### stdpopsim HomSap catalogue (coalescent default)

`OutOfAfrica_3G09` (CEU / YRI / CHB) drives the demographic history by
default. `Africa_1T12`, `OutOfAfrica_2T12`, etc. are also available.
First coalescent run may download model metadata; subsequent runs are
offline.

---

## Output layout

```
out/
├── person_0001.vcf.gz
├── person_0001.vcf.gz.tbi
├── person_0002.vcf.gz
├── ...
└── summary/
    └── sfs.tsv          # cohort AC histogram (columns: ac, n_sites)
```

Per-person VCF:

- **Header**: `VCFv4.2`, `##reference=<accession URL>`, `##contig` lines
  for every standard chromosome with `assembly=GRCh37` or `GRCh38`, full
  `INFO`/`FORMAT`/`ALT` declarations (AC, AN, AF, SVTYPE, SVLEN, END,
  CIPOS, HIGHLIGHT, CLNSIG, CLNDN; GT/DP/GQ/AD; `<DEL>`/`<DUP>`/
  `<INV>`/`<INS>`).
- **Records**: one ClinVar-highlighted variant flagged `HIGHLIGHT`
  (plus `CLNSIG` / `CLNDN` when present) + the cohort background
  projected to that person (hom-ref calls dropped).
- **Per-record FORMAT** (`GT:DP:GQ:AD`): DP ~ Poisson(λ ≈ 30, per-sample
  jitter σ = 3), AD ~ Binomial(DP, genotype_p) with ~5% ref-bias on
  hets, support-weighted Phred GQ capped [0, 99].
- **Multi-allelic** (legacy path): comma-separated ALT / AC / AF / AD;
  samples can carry `1|2`-style hets.

Compatibility: output passes `nextflow_pipeline/bin/qc_validate.py
--strict` cleanly (GRCh37/GRCh38 recognised, human contigs, GT/DP/GQ/AD
declared, AF + AC/AN present).

---

## Architecture

```
synthetic_people/
├── SYHTHETIC_PROJECT.md
├── IMPLEMENTATION_PLAN.md
├── README.md
├── requirements.txt
├── generate_people.py        # 15-line shim → syntheticgen.cli:main
├── syntheticgen/
│   ├── builds.py             # GRCh37 + GRCh38 contig tables, ClinVar URLs
│   ├── clinvar.py            # ClinVar fetch + pathogenic candidate load
│   ├── background.py         # 1000G coordinate pool loader (reservoir)
│   ├── cohort.py             # M4 shared-site cohort + haplotype slotting
│   ├── sfs.py                # M4 P(k) ∝ 1/k^α sampler + histogram
│   ├── coalescent.py         # M5 msprime + stdpopsim driver
│   ├── titv.py               # M3+ Ti/Tv calibrator for de-novo SNVs
│   ├── quality.py            # M2 DP / GQ / AD simulation
│   ├── header.py             # VCF header assembly
│   ├── writer.py             # bgzip + tabix single-sample write
│   └── cli.py                # argparse + orchestration
├── tests/
│   ├── test_quality.py       # M2 + M3 (generalised to N alleles)
│   ├── test_multiallelic.py  # M3
│   ├── test_titv.py          # M3+
│   ├── test_sfs.py           # M4
│   ├── test_cohort.py        # M4
│   └── test_coalescent.py    # M5 (skips cleanly without msprime/stdpopsim)
└── out/                      # generated VCFs + summary/
```

---

## Functionality by milestone

### M1 — package scaffolding, GRCh38 default, rich header

Split the original single-file `generate_people.py` into the
`syntheticgen/` package (`generate_people.py` is now a thin shim).
Default `--build` flipped to GRCh38; GRCh37 retained. `##contig` lines
carry `assembly=GRCh38` (or `GRCh37`). FORMAT declares the full
GT/DP/GQ/AD tag set from day one; INFO declares SVTYPE / SVLEN / END /
CIPOS plus symbolic `<DEL>`/`<DUP>`/`<INV>`/`<INS>` ALTs — populated in
later milestones, declared early so the header stops changing shape.
`--check-deps` audits htslib binaries and Python deps.

### M2 — per-variant quality metrics (DP / GQ / AD)

`syntheticgen/quality.py` simulates:

- **DP** ~ Poisson(λ = 30, per-sample jitter σ = 3).
- **AD** ~ Binomial(DP, p) with p = 0.0 / 0.475 / 1.0 by genotype (0.475
  on hets models empirical ref-bias in short-read WGS).
- **GQ** — support-weighted Phred, capped at 99, with depth-dependent
  ceiling `10 · log10(DP) · 6`.

Writer emits `GT:DP:GQ:AD`. 5-person batch verified: DP mean 29.4,
AD `sum == DP` on 100% of rows, het alt-fraction 0.475.

### M3 — indels, multi-allelic, Ti/Tv in target range

Indels flow through unchanged from the 1000G source (already
left-aligned, parsimonious — verified no indel has a common prefix
> 1 char across a 50-person batch).

Multi-allelic opened end-to-end: loader keeps per-allele AFs, draw
samples two haplotypes from a categorical over
`{REF, alt_1, …, alt_k}` so `1|2` hets can occur, per-allele AC / AF /
AD, `Number=R` for AD. 50-person chr22 batch Ti/Tv = 2.11 naturally
from the 1000G source.

### M3+ — Ti/Tv calibrator

`syntheticgen/titv.py` — `choose_alt(ref, rng, target=2.1)` draws a
non-REF base weighted so long-run Ti/Tv converges on `target`
(transition partner weight `target`; each transversion weight `0.5`).
Landed defensively ahead of M5, where de-novo SNV generation drops the
unbiased ratio to ~0.5. `is_transition` / `titv_ratio` helpers for
downstream validation.

### M4 — cohort-level generation + power-law SFS

Pivot from per-person independent HWE draws to one-pass cohort
simulation:

- `syntheticgen/sfs.py` — `draw_minor_count(n_hap, α)` samples k ∈
  {1, …, 2N-1} with `P(k) ∝ 1/k^α` (default α = 2.0). Steeper than
  Watterson (α = 1.0) to match gnomAD-like singleton-dominated spectra.
  `draw_allele_counts` handles multi-allelic via rejection so total
  AC ≤ 2N-1.
- `syntheticgen/cohort.py` — `assign_haplotypes` places alt alleles
  into specific 2N haplotype slots without replacement. Diploid GTs per
  person come from pairing consecutive slots, so realised AC matches
  drawn AC exactly (no HWE-resampling smoothing at the site level).

SFS histogram persisted to `out/summary/sfs.tsv`; singleton count +
fraction printed to the run log. 50-person chr22 legacy batch: 317
singletons / 511 alt observations = **62% singleton fraction** (clears
the >50% exit threshold).

### M5 — coalescent backbone (msprime + stdpopsim)

`syntheticgen/coalescent.py` drives msprime through stdpopsim's engine;
the chosen demographic model supplies population-size history and
per-chromosome metadata. REF/ALT bases for binary tree-sequence
mutations are synthesised via the M3+ Ti/Tv calibrator. Output matches
the M4 cohort-site dict shape, so the writer is unchanged.

New flags: `--chromosomes`, `--chr-length-mb`, `--demo-model`
(default `OutOfAfrica_3G09`), `--population` (default CEU), `--rec-rate`,
`--mu`, `--legacy-background`.

200-sample × 10 Mb chr22 exit check: 28 054 variable sites in 16 s;
8 805 common (MAF ≥ 5%); monotonic LD decay:

| distance bin | mean r² |
|---|---|
| 100–500 bp | 0.55 |
| 0.5–1 kb | 0.46 |
| 1–5 kb | 0.35 |
| 5–20 kb | 0.20 |
| 20–100 kb | 0.05 |
| 100–500 kb | 0.006 |
| ≥ 500 kb | <0.003 |

r² < 0.1 reached by ~20 kb — well inside the "<0.1 by 1 Mb" plan
threshold. Short-range anchor sits at ~0.5 rather than the plan's 0.9
because recombination is uniform; wiring in `HapMapII_GRCh38` hit a
stdpopsim "missing data" error on sub-chromosome regions and is deferred
to M6/M10.

---

## Test suite

71 tests across six files; all passing with deps installed.

```bash
cd synthetic_people && ../.venv/bin/python -m unittest discover -s tests -v
```

Without msprime/stdpopsim installed, `test_coalescent.py` skips cleanly
and **61/61** remaining still pass.

| File | Count | Coverage |
|---|---|---|
| `test_quality.py` | 12 | Poisson DP distribution, bi/multi-allelic AD (`sum==DP`, ref-bias on `0\|1` / `0\|2`, 50/50 split on `1\|2`), GQ range, depth-dependent cap, genotype consistency |
| `test_multiallelic.py` | 5 | Bi-allelic HWE reduction, multi-allelic categorical, `1\|2`-style het rate, per-alt dosage vectors |
| `test_titv.py` | 14 | Transition-partner table, transversion enumeration, case-insensitivity, `titv_ratio` corner cases (empty / no-Tv / indel skip), `choose_alt` uniformity, convergence at targets 0.5 / 1.0 / 2.1 / 3.0 (±5%), parameter validation |
| `test_sfs.py` | 16 | `draw_minor_count` range + near-uniform at α→0, singleton fraction >55% at default α = 2.0, `draw_allele_counts` total bound, histogram aggregation, TSV round-trip, parameter validation |
| `test_cohort.py` | 14 | `assign_haplotypes` exact-count preservation, random-slot placement, overflow rejection, cohort reproducibility under seed, every-site-variable invariant, coord-sharing across people, hom-ref drop-out |
| `test_coalescent.py` | 10 | Output shape, monotone positions, realised AC = declared AC, no fixed sites, seed reproducibility, Ti/Tv ∈ [1.7, 2.6], multi-chromosome, error on unknown chromosome, stdpopsim end-to-end (`skipUnless` on msprime/stdpopsim import) |

Per-milestone exit check: `nextflow_pipeline/bin/qc_validate.py --vcf
<person.vcf.gz> --name <id> --out <out.json> --strict` (exit 1 on any
hard failure).

---

## Known gaps

Tracked in `IMPLEMENTATION_PLAN.md`:

- **M6** — UK-cohort admixture (EUR + SAS + AFR pulses, per-person local
  ancestry BED).
- **M7** — ClinVar / dbSNP / COSMIC overlays on the coalescent output.
- **M8** — structural variants (symbolic ALTs, SVTYPE / SVLEN / END).
- **M9** — sequencer-noise / genotyping-error injection.
- **M10** — validation suite (PCA vs 1000G, LD decay curves, stats).
- **M11** — delivery packaging (manifest, truth sets, smoke script).

## CLI reference

| Flag | Purpose | Default |
|---|---|---|
| `--n` | Cohort size | `10` |
| `--seed` | RNG seed; omit for fresh randomness each run | `None` |
| `--build` | `GRCh37` or `GRCh38` | `GRCh38` |
| `--output-dir` | Per-person VCF output | `./out` |
| `--cache-dir` | ClinVar download cache | `./cache` |
| `--check-deps` | Print dependency status and exit | `False` |
| `--clinvar-sig` | Comma-separated CLNSIG values | `Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic` |
| `--chromosomes` | [coalescent] Comma-separated chroms | `22` |
| `--chr-length-mb` | [coalescent] Simulated prefix per chrom | `5.0` |
| `--demo-model` | [coalescent] stdpopsim model id; `none` for uniform | `OutOfAfrica_3G09` |
| `--population` | [coalescent] Sampling population | `CEU` |
| `--rec-rate` | [coalescent, `--demo-model=none`] Uniform recombination rate | `1e-8` |
| `--mu` | [coalescent, `--demo-model=none`] Mutation rate | `1.29e-8` |
| `--legacy-background` | Use M4 1000G-pool + power-law SFS sampler | `False` |
| `--background-glob` | [legacy] Source glob(s) for common variants | 1000G files in parent dir |
| `--n-background` | [legacy] Shared background site count | `500` |
| `--af-min` | [legacy] Minimum AF when loading the pool | `0.05` |
| `--sfs-alpha` | [legacy] Power-law exponent for the SFS | `2.0` |
