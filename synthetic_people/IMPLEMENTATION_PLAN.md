# High-Fidelity Synthetic Generator — Implementation Plan

Plan for upgrading `synthetic_people/generate_people.py` from its current
stdlib-only, independent-site HWE sampler into the high-fidelity,
admixture-aware, LD-correct generator specified in `SYHTHETIC_PROJECT.md`.

Progress is tracked by checkboxes on each milestone. Tick off items as work
lands.

---

## Gap analysis — what the current generator does vs. what the spec needs

| Area | Current (`generate_people.py`) | Spec requirement |
|---|---|---|
| Output scope | One VCF per person, one "highlighted" ClinVar variant + independent HWE background | Cohort-level, LD-correct, admixed genomes |
| FORMAT tags | `GT` only | `GT`, `DP`, `GQ`, `AD` |
| Variant types | SNV-like + simple small indels, biallelic only | SNPs, left-aligned indels, **multi-allelic**, **SVs** with `SVTYPE` |
| Header | Basic contigs (no `assembly=`), minimal `INFO` | Full contigs with `assembly=GRCh38`, all `INFO`/`FORMAT` declared |
| Reference build | Defaults to GRCh37 | **GRCh38** required |
| Linkage | None (sites independent) | LD blocks via coalescent (`msprime`) or HMM |
| Mutation metrics | Not enforced | Ti/Tv ≈ 2.1, power-law SFS (gnomAD-like), dbSNP/ClinVar/COSMIC grounding |
| Population | Single anonymous individual | Mixed UK cohort: EUR + SAS + AFR with admixture pulses |
| Ancestry truth | None | Per-individual local ancestry BED |
| Error model | None | Sequencer noise (ART/SimNGS or equivalent) |
| Validation | None | PCA vs 1000G, LD decay curve, Ti/Tv + Het/Hom + AF report |
| Dependencies | Python stdlib + htslib binaries | + numpy, msprime, tskit, stdpopsim, matplotlib |

**Verdict:** incremental rewrite is feasible. The single-file script becomes
a small package (`synthetic_people/`) with cohort-level simulation at its
core. The existing ClinVar overlay logic is re-used in the database-grounding
milestone.

---

## Dependency strategy

The spec explicitly recommends `msprime` and `stdpopsim`, which are the
industry-standard way to get realistic LD and admixture. There is no
reasonable pure-stdlib path to those goals. Therefore:

* `requirements.txt` gains: `numpy`, `msprime`, `tskit`, `stdpopsim`,
  `matplotlib`, `scikit-allel` (for PCA / LD decay), optionally `pandas`.
* htslib binaries (`bcftools`, `tabix`, `bgzip`) remain required.
* A `pip install -r requirements.txt` instruction lands in `README.md`.
* Each milestone after M1 is gated on this dependency bump landing.

---

## Target package layout (by the end of M1)

```
synthetic_people/
├── SYHTHETIC_PROJECT.md        # spec (existing, read-only)
├── IMPLEMENTATION_PLAN.md      # this file
├── README.md                   # updated per milestone
├── requirements.txt            # expanded per M1
├── generate_people.py          # thin CLI entry point (re-exports main)
├── syntheticgen/               # new package
│   ├── __init__.py
│   ├── cli.py                  # arg parsing, top-level orchestration
│   ├── builds.py               # contig tables, reference URLs, assembly tags
│   ├── header.py               # VCF header assembly (contigs, INFO, FORMAT)
│   ├── writer.py               # bgzip+tabix write, record formatting
│   ├── clinvar.py              # ClinVar fetch + pathogenic candidate load
│   ├── dbsnp.py                # dbSNP overlay (M7)
│   ├── cosmic.py               # COSMIC overlay (M7, optional)
│   ├── background.py           # legacy 1000G background sampler (kept for M1-M3)
│   ├── coalescent.py           # msprime / stdpopsim driver (M5+)
│   ├── admixture.py            # UK-cohort admixture + local ancestry BED (M6)
│   ├── quality.py              # DP / GQ / AD simulation (M2)
│   ├── sv.py                   # structural variant generator (M8)
│   ├── errors.py               # sequencing-error injection (M9)
│   └── validate.py             # PCA, LD decay, Ti/Tv, Het/Hom plots (M10)
└── out/                        # VCFs + BED truth + ancestry BED + plots
```

---

## Milestones

Each milestone is an independently shippable increment. The script remains
runnable at every step. Milestones are ordered so that validation (M10) can
sensibly be run against the output of earlier milestones — e.g. Ti/Tv and
SFS checks become meaningful at M4, LD decay at M5, PCA at M6.

### M1 — Scaffolding, GRCh38 default, rich header

- [x] Split `generate_people.py` into the package layout above. `generate_people.py`
      becomes a 3-line shim that calls `syntheticgen.cli.main`.
- [x] Default `--build` flips to **GRCh38**. GRCh37 stays selectable.
- [x] Contig lines include `assembly=GRCh38` (or `GRCh37`) per spec §2.1.
- [x] Header declares `FORMAT` for `GT`, `DP`, `GQ`, `AD` (values still
      placeholder until M2) and `INFO` for `SVTYPE`, `SVLEN`, `END` (unused
      until M8).
- [x] `pip install -r requirements.txt` documented in README; CI-friendly
      `--check-deps` flag added.
- [x] Existing one-highlighted-variant-plus-HWE-background behaviour
      preserved byte-for-byte under a seeded run — verified via
      `diff` on record lines between the pre-refactor baseline and
      `--build GRCh37 --seed 42` output.

**Exit check:** `python3 generate_people.py --n 3 --seed 42 --build GRCh38`
produces 3 valid GRCh38 VCFs that pass `bcftools view -h` and
`bin/qc_validate.py`. ✅ Verified 2026-04-24: `qc_validate.py --strict`
exits 0 with 0 errors / 0 warnings; reference build recognised as
`GRCh38`; header declares full FORMAT + SV INFO tag set.

---

### M2 — Realistic per-variant quality metrics (DP / GQ / AD)

- [x] Add `syntheticgen/quality.py` implementing:
      * `DP` ~ Poisson(λ=30) with per-sample λ jitter.
      * `AD` from a Binomial(DP, p) where p reflects genotype
        (0.0 / 0.5 / 1.0) plus realistic sampling bias (~5% ref-bias
        for hets, so not 50/50 exactly).
      * `GQ` decreasing with DP variance; floor at 0, cap at 99.
- [x] Per-record FORMAT becomes `GT:DP:GQ:AD`.
- [x] Write a small unit-test file (`tests/test_quality.py`) that drives
      distributions and checks means/variances. 12 tests, all passing.

**Exit check:** `bcftools stats` on a fresh batch reports DP mean ≈ 30,
GQ distribution non-degenerate, AD sum == DP on every row. ✅ Verified
2026-04-24 on a 5-person / 1,289-record batch: DP mean 29.41 (13–51),
GQ span 1–99 (184 high-confidence ≥90), AD mismatches 0/1,289, het
alt-fraction 0.475 (confirming ref-bias). `qc_validate.py --strict`
remains clean.

---

### M3 — Variant diversity: indels, multi-allelic, Ti/Tv calibration

- [x] Indels (1–50 bp) flow through cleanly from the 1000G source; they
      arrive already left-aligned and parsimonious per VCF 4.2, and the
      pool filter preserves that. Verified no indel has a common prefix
      > 1 char across a 50-person batch.
- [x] Multi-allelic records enabled: loader now passes `MAX(INFO/AF)` and
      keeps per-allele AFs; draw path samples two haplotypes from a
      categorical over {REF, alt_1, ..., alt_k}; writer emits
      `ALT=A,G` with comma-separated per-allele `AC` / `AF` / `AD`.
- [x] Ti/Tv "calibrator" not needed: natural sampling from 1000G Phase 3
      already lands at **2.11** on a 50-person batch, dead-centre of the
      [1.9, 2.3] target. Skipped the calibrator per the working-agreement
      to avoid premature abstraction.

**Exit check:** ✅ 2026-04-24 on 50-person / 15,587-record batch:
Ti/Tv = 2.11; 2,148 indels (13.8%), all left-aligned; 313 multi-allelic
sites (2.0%); 23/23 unit tests pass; `qc_validate.py --strict` clean.

---

### M4 — Site Frequency Spectrum + cohort-level generation

- [ ] Switch from per-person independent sampling to **cohort generation**:
      simulate an N-sample cohort once, write one VCF per person with the
      same variant coordinates, different genotypes.
- [ ] Draw site frequencies from a power-law / neutral-coalescent SFS
      (~1/f) rather than uniform, so singletons dominate.
- [ ] Report SFS in the run log; persist a histogram to
      `out/summary/sfs.tsv`.

**Exit check:** fraction of singletons > 50% of variable sites; AF spectrum
plot shows the expected 1/f shape.

---

### M5 — Coalescent backbone (single population, with LD)

- [ ] Add `syntheticgen/coalescent.py` wrapping `msprime.sim_ancestry` +
      `msprime.sim_mutations` on a `stdpopsim` human genetic map
      (HapMapII_GRCh38 or deCODE).
- [ ] For each chromosome simulated, emit tree sequence → tskit → per-sample
      genotypes written through the existing VCF writer.
- [ ] Provide `--demo-model` knob (default: OutOfAfrica single-pop draw, for
      now).
- [ ] Keep background-from-1000G path behind a `--legacy-background` flag so
      old behaviour remains testable.

**Exit check:** LD decay plot (M10 preview) produced from a 200-sample batch
shows r² falling from ~0.9 at 1 kb to <0.1 by 1 Mb — matches real human
expectations.

---

### M6 — UK-cohort admixture + local ancestry truth BED

- [ ] Use `stdpopsim`'s 3-population human models to simulate EUR + SAS +
      AFR separately, then apply an "admixture pulse" using
      `msprime.demes` / `demes` to mix at user-configurable proportions
      (default 60/25/15 reflecting London-ish demographics — make
      configurable, do not hard-code).
- [ ] Record **local ancestry** per haplotype segment; write one
      `out/ancestry/person_<N>.bed` per individual:
      `chrom  start  end  ancestry_haplotype1  ancestry_haplotype2`.
- [ ] Sample IDs still HG/NA-prefixed; add an ancestry-summary column in
      the batch manifest (`EUR=0.61,SAS=0.24,AFR=0.15`).

**Exit check:** PCA plot (M10) shows synthetic samples occupying the
EUR–SAS–AFR triangle in 1000G reference space, with admixed individuals
lying between clusters proportionally to their ancestry fractions.

---

### M7 — ClinVar / dbSNP / COSMIC grounding

- [ ] Refactor existing ClinVar overlay into `syntheticgen/clinvar.py`. For
      any coalescent-produced variant whose (chrom,pos,ref,alt) tuple
      collides with a ClinVar record, copy `CLNSIG` / `CLNDN` onto the
      output record (don't force the overlay — keep positions consistent
      with the simulation).
- [ ] Add `syntheticgen/dbsnp.py`: download the dbSNP GRCh38 VCF (cached
      like ClinVar) and **inject** a small set of known variants into the
      cohort (annotated with rsIDs) so `ID` columns are non-empty and
      realistic rsID densities are hit.
- [ ] Optional `syntheticgen/cosmic.py` (gated behind `--somatic` flag;
      COSMIC registration is required, so it stays optional).

**Exit check:** a significant fraction of output records carry rsIDs from
dbSNP; ClinVar-known positions are annotated; `bcftools view -i 'CLNSIG!="."'`
returns the injected pathogenic variants.

---

### M8 — Structural variants

- [ ] Add `syntheticgen/sv.py`: generate a handful of SVs per individual
      (deletions 50 bp – 10 kb, tandem duplications, simple inversions).
      Emit with proper `SVTYPE`, `SVLEN`, `END`, `CIPOS`, and symbolic
      ALTs (`<DEL>`, `<DUP>`, `<INV>`).
- [ ] Update header with `ALT=<ID=DEL,Description=...>` lines etc.

**Exit check:** `bcftools view -i 'INFO/SVTYPE!="."'` returns structured SV
records; `bcftools stats` counts them in its SV section.

---

### M9 — Sequencing error modelling

- [ ] Lightweight path (default): inject per-call genotyping noise
      probabilistically — occasional GT flips, lowered GQ, depth dropouts
      — parameterised by a target FDR (~0.1%).
- [ ] Heavy path (optional, gated by `--art` flag): shell out to `ART` to
      simulate reads from the simulated genome and call back with
      `bcftools call`, closing the loop. Disabled by default because ART
      isn't on the typical bioinformatics box.

**Exit check:** noisy batch shows expected FDR when cross-checked against
the coalescent truth (M5 provides the "truth" comparator).

---

### M10 — Validation suite (PCA, LD decay, stats)

- [ ] `syntheticgen/validate.py` + a `validate_batch.py` CLI producing:
      * **PCA plot** using `scikit-allel` — project synthetic samples into
        a 1000G PCA space (pre-computed from the local chr19–22 VCFs).
      * **LD decay curve**: r² vs. physical distance in kb, split by
        population, on log-x.
      * **Variant statistics report** (HTML or Markdown):
        Ti/Tv, Het/Hom per sample, AF histogram, singleton fraction, indel
        length distribution, SV count summary.
- [ ] Plots land under `out/validation/`; the Markdown report links them.

**Exit check:** `python3 validate_batch.py out/` succeeds end-to-end and
produces all three artefacts; acceptance criteria from spec §6 are visibly
met on the plots.

---

### M11 — Delivery packaging

- [ ] Per-batch manifest `out/manifest.json` listing each person with:
      ancestry fractions, highlighted variants, file paths.
- [ ] Truth set BED: per person, a BED listing "golden" variant positions
      (injected ClinVar + simulated true variants) vs. noise-introduced
      calls from M9.
- [ ] `README.md` final rewrite describing the new pipeline end-to-end,
      including how to regenerate the validation report.
- [ ] Smoke-test script `scripts/smoke.sh` running `--n 20` end-to-end in
      CI-friendly time (<5 min).

**Exit check:** a fresh clone can `pip install -r requirements.txt`,
run `python3 generate_people.py --n 20 --seed 42`, and have all spec
deliverables (VCF.gz, tbi, truth BED, local ancestry BED, validation plots,
stats report) land under `out/` without manual intervention.

---

## Out of scope (for now)

* Somatic VCFs with tumour-normal pairs — feasible extension but not
  explicitly required beyond the optional COSMIC hook.
* Multi-sample joint VCFs — current design is one-file-per-person; joint
  output can be added at the end as a post-processing step (`bcftools
  merge`) if needed.
* Real BAM/FASTQ generation — `--art` mode produces reads only as an
  internal step to simulate sequencer errors; full FASTQ delivery is not a
  goal of this spec.

---

## Risk register

| Risk | Mitigation |
|---|---|
| `msprime`/`stdpopsim` installs fail on the user's machine | Pin known-good versions in `requirements.txt`; document conda as fallback |
| Full-genome coalescent simulation is slow / RAM-hungry | Default to chr19–22 (matching the local 1000G data already on disk); expose `--chromosomes` knob |
| dbSNP download is huge (~25 GB full, ~1 GB common-only) | Default to dbSNP "common" subset; make size trade-off visible in the CLI |
| Validation plots require 1000G reference already on disk | Re-use the chr19–22 VCFs already present in the repo root |

---

## Working agreement

* Every milestone lands as its own commit (or small PR-sized stack) so the
  history reads as an incremental evolution, not a single mega-rewrite.
* Each milestone's **exit check** must pass before the next starts.
* The existing, simpler generator continues to work under a
  `--legacy-background` flag through at least M7, so downstream users have
  a bail-out while the new pipeline is maturing.
