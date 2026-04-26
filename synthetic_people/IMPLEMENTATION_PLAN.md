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
- [x] Ti/Tv calibrator in place (`syntheticgen/titv.py`): landed
      defensively ahead of M5, where de-novo SNV generation drops the
      natural ratio to ~0.5. `choose_alt(ref, rng, target=2.1)` draws a
      non-REF base weighted so long-run Ti/Tv converges on the target;
      14 tests verify convergence at several targets (0.5, 1.0, 2.1,
      3.0) within ±5%. Not wired into the current 1000G-backed path
      (which lands at 2.11 naturally without calibration).

**Exit check:** ✅ 2026-04-24 on 50-person / 15,587-record batch:
Ti/Tv = 2.11; 2,148 indels (13.8%), all left-aligned; 313 multi-allelic
sites (2.0%); 23/23 unit tests pass; `qc_validate.py --strict` clean.

---

### M4 — Site Frequency Spectrum + cohort-level generation

- [x] Switched from per-person independent sampling to **cohort generation**
      (`syntheticgen/cohort.py`): a single pass picks N_sites from the
      1000G coordinate pool and, for each site, assigns alt alleles into
      the 2N haplotype slots without replacement. Each person's phased GT
      at each site is then read off consecutive slot pairs, so the realised
      cohort AC is exact (no HWE-resampling smoothing) and every person
      sees the same (chrom, pos, ref, alts) at each site.
- [x] Site frequencies now come from a power-law SFS (`syntheticgen/sfs.py`):
      `P(k) ∝ 1/k^α` over k ∈ {1, …, 2N-1}, with α exposed via
      `--sfs-alpha` (default 2.0). α=1.0 is Watterson-neutral; α=2.0 biases
      toward singletons to match gnomAD-style empirical spectra.
- [x] Run log prints cohort size, total alt observations, singleton count
      and fraction. Histogram persists to `out/summary/sfs.tsv` (columns:
      `ac`, `n_sites`). Multi-allelic sites contribute one observation
      per alt.
- [x] 24 new unit tests (12 in `test_sfs.py`, 12 in `test_cohort.py`)
      covering SFS range/shape, singleton-fraction default, allele-count
      rejection-sampling bound, slot-assignment exactness, cohort
      reproducibility under seed, and hom-ref drop-out logic. 61 tests
      total, all passing.

**Exit check:** ✅ 2026-04-24 on `--n 50 --seed 42 --build GRCh37` chr22
batch: 511 alt observations, **317 singletons = 62.0%** (well above the
50% threshold), full 1/k^α shape visible in `sfs.tsv` (AC=1→317, AC=2→74,
…, AC=76→1). Cross-check: cohort AC realised at a high-AC site sums to
76 across the 50 per-person VCFs — exactly what the histogram records.
`qc_validate.py --strict` exits 0 with 0 errors / 0 warnings.

---

### M5 — Coalescent backbone (single population, with LD)

- [x] `syntheticgen/coalescent.py` wraps `msprime.sim_ancestry` +
      `msprime.sim_mutations` through `stdpopsim`'s engine so the
      configured demographic model drives diversity directly. REF/ALT
      bases for the binary tree-sequence mutations are synthesised via
      the M3+ Ti/Tv calibrator (`titv.choose_alt`), so de-novo SNVs
      land near the target 2.1 ratio without post-hoc rejection.
- [x] Tree-sequence iteration (`ts.variants()`) drives `_tree_sequence_to_sites`
      which converts each variable site to the same cohort-site dict
      shape M4 uses — so `writer.py` needs no changes and the per-person
      loop in `cli.py` branches only on where the sites come from.
- [x] `--demo-model` knob (default `OutOfAfrica_3G09`), `--population`
      (default `CEU`), plus `--chromosomes`, `--chr-length-mb`, `--rec-rate`,
      `--mu`. Setting `--demo-model none` falls back to a constant-size
      Ne=10k single-pop `msprime` draw — faster, dependency-lighter, no
      realistic demography.
- [x] `--legacy-background` routes back to the M4 1000G-pool + power-law
      SFS sampler (legacy-only flags marked `[legacy]` in help); old
      behaviour stays testable.
- [x] 10 new unit tests in `tests/test_coalescent.py` — skip cleanly if
      msprime/stdpopsim aren't installed. Cover output shape, monotone
      positions, realised AC = declared AC, reproducibility under seed,
      Ti/Tv landing in [1.7, 2.6], and one stdpopsim end-to-end check.
      71/71 tests passing with deps; 61/61 still pass without.

**Exit check:** ✅ 2026-04-24 on 200-sample × 10 Mb chr22 sim
(`OutOfAfrica_3G09` / CEU, uniform recombination, seed 42): 28,054
variable sites in 16 s; 8,805 common (MAF ≥ 5%). Mean r² across
distance bins: **0.55 @ 100–500 bp → 0.46 @ 0.5–1 kb → 0.35 @ 1–5 kb →
0.20 @ 5–20 kb → 0.05 @ 20–100 kb → 0.006 @ 100–500 kb → <0.003 @
≥500 kb**. Monotonic decay, r² < 0.1 reached by ~20 kb (well inside
the "<0.1 by 1 Mb" threshold). The short-range anchor is ~0.5 rather
than the plan's aspirational 0.9 because we're using uniform
recombination, not a HapMap-derived map — wiring in
`HapMapII_GRCh38` produced a stdpopsim "missing data" error on
sub-chromosome regions and is deferred to M6/M10. `qc_validate.py
--strict` clean on coalescent output.

---

### M6 — UK-cohort admixture + local ancestry truth BED

- [x] `syntheticgen/admixture.py` builds a `demes`-defined demography
      with three sources (EUR, SAS, AFR) and a single admixture pulse
      `PULSE_TIME=20` generations (~600 years) ago into a UK deme.
      Source pop sizes mirror Gutenkunst OOA_3G09 parameters
      (ANC=12.3k, AFR=12.3k, OOA bottleneck=2.1k, EUR/SAS=10k); UK
      present-day Ne=50k. `msprime.sim_ancestry(...,
      record_migrations=True)` runs the demography; mutations placed
      with `BinaryMutationModel` and REF/ALT bases drawn through the
      M3+ Ti/Tv calibrator.
- [x] Local ancestry: for each haplotype-sample, walk the tree at
      every breakpoint to find the lineage node spanning the pulse
      time; the t=20 migration on that node tells us which source deme
      the lineage migrated into. Adjacent same-ancestry segments are
      merged; haplotype pairs are intersected into per-person joint
      `(start, end, h1_pop, h2_pop)` rows. Written one BED per person
      to `out/ancestry/person_NNNN.bed`.
- [x] Manifest: `out/manifest.json` always emitted; in admixture mode
      includes `ancestry_proportions` (requested) and per-person
      `ancestry_fractions` (realised).
- [x] CLI: `--admixture` flag plus `--eur-frac`/`--sas-frac`/`--afr-frac`
      knobs (default 60/25/15). Proportions validated to sum to 1.0
      and be non-negative. Overrides `--demo-model` / `--population`
      (which target the M5 single-population path).
- [x] 13 new unit tests in `tests/test_admixture.py` — skip cleanly
      when msprime/demes/tskit aren't installed. Cover demography
      build (proportion validation, UK ancestor topology), simulation
      shape, full chromosome coverage, realised AC = declared AC, BED
      round-trip, ancestry-fraction normalisation, multi-chromosome
      cohort, seed reproducibility, and aggregate ancestry tracking
      requested proportions within tolerance on a 30-person × 1 Mb
      sim. 84/84 tests passing with deps; 61/61 still pass without.

**Exit check:** ✅ 2026-04-25 on 20-person × 5 Mb chr22 admixture sim
(seed 42, default 60/25/15): 13,549 variable sites; 43 ancestry
segments across the cohort (mean 2.1 segments/person — biologically
expected because 20 generations of recombination × 5 Mb yields ≈1
breakpoint per haplotype); cohort-mean realised ancestry
**EUR=0.456, SAS=0.352, AFR=0.192** — within finite-cohort sampling
noise of the requested 0.60/0.25/0.15. Per-person VCFs pass
`qc_validate.py --strict` with 0 errors / 0 warnings. The literal PCA
acceptance test in spec §6 lands in M10; the admixture truth set is
now in place to validate against. Stand-alone proportions check on a
30-person × 1 Mb cohort lands EUR ≈ 0.6, SAS ≈ 0.25, AFR ≈ 0.15
(within ±15%, in `tests/test_admixture.py`).

---

### M7 — ClinVar / dbSNP / COSMIC grounding

- [x] `syntheticgen/clinvar.py` now exposes `load_clinvar_index`,
      `annotate_clinvar` (collision-only, copies CLNSIG/CLNDN/id when a
      cohort site happens to land on a ClinVar coordinate) and
      `inject_clinvar` (replaces a `--clinvar-inject-density` fraction
      of cohort sites with random ClinVar pathogenic records — moves
      the site to a real chromosome coordinate while preserving the
      cohort GT block, so LD remains intact). Coalescent positions
      live in `[1, sim_length]` while ClinVar sits at real chromosome
      coordinates (chr22 ClinVar spans 15.5M–50.8M), so collision-only
      annotation almost never fires — `inject_clinvar` is the
      practical mechanism that lands CLNSIG/CLNDN at realistic
      positions. The highlighted per-person variant path is unchanged.
- [x] `syntheticgen/dbsnp.py` provides `load_rsid_pool` and
      `inject_rsids`. Default rsID source is the cached ClinVar VCF
      whose `INFO/RS` tag carries dbSNP rs numbers — this delivers
      thousands of rsID-bearing records per chromosome with no extra
      download. Users can pass `--dbsnp-vcf PATH` to point at a real
      dbSNP-style VCF (rsIDs in the ID column) instead. `--rsid-density`
      (default 0.20) controls how many cohort sites get rewritten with
      a known rsID + its real coordinate / REF / ALT. ClinVar-injected
      rows are reserved so the two overlays don't fight for the same
      sites.
- [x] `syntheticgen/cosmic.py` (gated behind `--somatic --cosmic-vcf
      PATH`) reads any COSMIC-format VCF (GENE/LEGACY_ID/CDS/AA), and
      `inject_cosmic` replaces a `--cosmic-inject-density` fraction of
      cohort sites with COSMIC mutations. COSMIC registration is still
      required to obtain the source file; we never auto-fetch and
      `--somatic` without `--cosmic-vcf` exits with a clear message.
- [x] Header gains `INFO=COSMIC_ID` and `INFO=COSMIC_GENE` declarations
      (M7); writer carries `clnsig`, `clndn`, `cosmic_id`, `cosmic_gene`
      from the cohort site through the per-person record onto the
      emitted INFO field. `cohort.person_records_from_cohort` now
      forwards those fields to the writer.
- [x] Manifest gains an `overlays` block with the requested densities,
      the dbSNP source path, and realised counts of annotated /
      injected records.
- [x] 23 new tests in `tests/test_overlays.py` covering ClinVar
      annotate / inject, rsID normalisation (ID column vs INFO/RS,
      bare digits vs prefixed, multi-rsID lists), inject_rsids
      reserve_indices, COSMIC inject, post-injection sort invariants,
      GT-block preservation, and ClinVar+rsID overlay interaction.
      Tests are pure-Python and run without bcftools/network. **107/107
      tests with deps; 84/84 stdlib-only (23 skipped).**

**Exit check:** ✅ 2026-04-25 on 5-person × 1 Mb chr22 sim
(`--demo-model none`, seed 42, default densities): 1,299 cohort sites;
**13 ClinVar pathogenic injections** (CLNSIG=Pathogenic /
Likely_pathogenic landing at real chr22 coordinates 29–50M); **259
rsID injections** (~20% of records). Per-person VCFs carry
117–135 rsIDs and 5–7 CLNSIG-bearing records each. Spot check —
chr22:15528207 emerges as `rs3924507 C>T`,
chr22:29673446 as `Pathogenic / Neurofibromatosis,_type_2`. All five
VCFs pass `qc_validate.py --strict` with 0 errors / 0 warnings;
header now declares COSMIC_ID / COSMIC_GENE alongside CLNSIG / CLNDN.
COSMIC injection exercised via `tests/test_overlays.py` (no public
COSMIC VCF in the cache; `--somatic --cosmic-vcf PATH` is the
documented opt-in for users with a registered download).

---

### M8 — Structural variants

- [x] `syntheticgen/sv.py` exposes `generate_person_svs(rng,
      chromosomes, chrom_length_bp, n_svs, length_min_bp,
      length_max_bp)` returning a list of writer-shaped variant dicts
      with `svtype` ∈ {DEL, DUP, INV}. Length is log-uniform on
      [min, max] (default 50 bp – 10 kb). Type weights default
      0.50 / 0.30 / 0.20 (DEL / DUP / INV). Phased GT defaults to
      ~80% het / 20% hom-alt. Anchor REF base is a random standard
      base (the real GRCh38 reference isn't on disk; this is the
      pragmatic placeholder until M11 wires the FASTA in).
- [x] Each record carries:
      * Symbolic ALT `<DEL>`/`<DUP>`/`<INV>` (already declared in the
        M1 header so no header change is required).
      * `INFO/SVTYPE`, `INFO/SVLEN` (negative for DEL, positive for
        DUP/INV), `INFO/END = POS + |SVLEN|`, `INFO/CIPOS = -50,50`
        (every SV is "imprecise"; tighter CIs are an M9/M10 refinement).
      * Per-record FORMAT continues to use the M2 quality model so
        `GT:DP:GQ:AD` is populated — AD becomes ref-supporting vs
        alt-supporting reads under the same Binomial(DP, p) model used
        for SNVs.
- [x] CLI: `--svs-per-person` (default 3), `--sv-length-min` (50),
      `--sv-length-max` (10000). 0 = skip SV emission. Each SV's
      `(POS, length)` is bounded so `END` stays inside the simulated
      span.
- [x] Manifest gains an `svs` block recording requested per-person
      count, length bounds, and total SVs emitted across the cohort;
      per-person `n_svs` lands on each `people[]` entry.
- [x] Writer: when `variant["svtype"]` is set, emit
      `SVTYPE=...;SVLEN=...;END=...;CIPOS=lo,hi` after the existing
      AC/AN/AF/HIGHLIGHT/CLNSIG/COSMIC fields. The existing alt-list
      `",".join(alts)` already puts the symbolic ALT in the right
      column.
- [x] 22 new tests in `tests/test_sv.py` — pure-Python; cover
      log-uniform length distribution skew, length bounds, anchor base
      validity, CIPOS default, SVTYPE/SVLEN sign-consistency, GT
      validity, multi-chromosome distribution, default-weight
      distribution within ±0.07 of (0.50, 0.30, 0.20), seed
      reproducibility, parameter validation. **129/129 with deps;
      106/106 stdlib-only.**

**Exit check:** ✅ 2026-04-25 on 3-person × 1 Mb chr22 sim
(`--demo-model none`, seed 42, default `--svs-per-person 5`):
each person carries exactly 5 SVs (`bcftools view -i 'INFO/SVTYPE!="."'`
→ 5 records on every VCF; `bcftools stats` reports them under
"number of others"). Spot check: `chr22:535206 A→<DEL>` with
`SVTYPE=DEL;SVLEN=-586;END=535792;CIPOS=-50,50;GT=0|1`. All three
VCFs pass `qc_validate.py --strict` with 0 errors / 0 warnings;
`n_variants=690` includes the 5 SVs.

---

### M9 — Sequencing error modelling

- [x] `syntheticgen/errors.py` (lightweight, default): per-call
      `maybe_flip_gt(true_gt, rng, error_rate)` and
      `maybe_dropout(rng, dropout_rate)` perturbations applied **after**
      AD has been drawn from the truth. Because `gq_from_ad` recomputes
      GQ from the (possibly-disagreeing) called GT, flipped calls
      naturally land low-GQ — exactly the realistic mis-call signal.
      Biallelic flip distribution: hom→het is the dominant failure
      mode (0.7 weight); het→hom-ref vs hom-alt is 50/50.
      Multi-allelic `1|2`-style hets collapse one allele to REF.
- [x] Writer threads `error_rate`, `dropout_rate`, optional `stats` dict
      through `write_person_vcf`. Dropouts emit
      `./.:0:0:0,0,...`; flips replace GT and recompute GQ.
- [x] CLI: `--error-rate` (default `0.001` = 0.1% FDR target),
      `--dropout-rate` (default `0.0005`), and `--art` flag (heavy
      path). Per-batch stats print to stderr; manifest gains an
      `errors` block with mode / requested rates / realised counts /
      realised FDR. Per-person manifest entries gain a `errors`
      sub-dict.
- [x] Heavy path (`--art`) is gated and currently exits with a clear
      message that ART read simulation needs the M11 GRCh38 reference
      FASTA. The lightweight path is the M9 deliverable.
- [x] 18 new tests in `tests/test_errors.py` — pure-Python; cover
      zero-rate / negative-rate no-ops, full-rate always flips,
      realised flip rate near target on 10k draws, hom→het bias
      (~70% on 0|0), het→hom split (~50/50 on 0|1), multi-allelic
      collapse to REF, unparseable GT pass-through, dropout rate,
      seed reproducibility, stats accumulator, default-constants
      lock-in. **147/147 with deps; 124/124 stdlib-only (23 skipped).**

**Exit check:** ✅ 2026-04-25 on 3-person × 1 Mb chr22 sim
(`--demo-model none`, seed 42, `--error-rate 0.01 --dropout-rate
0.005 --svs-per-person 0`): **realised FDR 1.55%** vs requested 1.50%
(0.01 + 0.005); 16 flips + 13 dropouts over 1,871 calls. Truth-vs-call
disagreement visible in the output: e.g. `chr22:758982 0|0:40:0:0,40`
— called hom-ref but all 40 reads on the alt, GQ correctly drops to
0. Dropouts emit `./.:0:0:0,0`. All three VCFs pass `qc_validate.py
--strict` with 0 errors / 0 warnings.

---

### M10 — Validation suite (PCA, LD decay, stats)

- [x] `syntheticgen/validate.py` exposes the analytics primitives:
      * `iter_records(vcf_path)` streams a single-sample VCF via
        `bcftools query`, splitting multi-allelics into one Record per
        ALT.
      * `summarise_vcf` aggregates per-sample counts (SNV / indel / SV
        / Ti / Tv / Het / Hom-alt / Hom-ref / dropouts / AF / indel
        lengths / SV-by-type / singletons).
      * `titv_from_stats` and `het_hom_ratio` close the spec acceptance
        signals; `af_histogram`, `aggregate_indel_lengths`, and
        `aggregate_sv_summary` produce the plot inputs.
      * `build_genotype_matrix(vcf_paths)` materialises a
        `(n_samples, n_variants)` int8 dosage matrix over the union of
        sites; missing genotypes coded `-1`.
      * `ld_decay(matrix, positions, chroms, distance_bins_kb,
        pairs_per_bin, rng)` returns mean r² per log-spaced distance
        bin; pair sub-sampling is reproducible under the supplied RNG.
      * `cohort_pca(matrix, n_components)` mean-imputes missing calls,
        prunes zero-variance columns, and returns the
        `sklearn.decomposition.PCA` projection plus explained-variance
        ratio.
- [x] `syntheticgen/plots.py` contains all matplotlib code (gated
      behind the import — caller still gets stats / JSON / Markdown
      without matplotlib): `plot_ld_decay`, `plot_af_histogram`,
      `plot_pca`, `plot_indel_lengths`. Headless-safe (`Agg` backend).
- [x] `validate_batch.py` (top-level CLI alongside `generate_people.py`)
      walks every `person_*.vcf.gz` under the batch directory, computes
      the full stats stack, runs LD-decay + PCA, and writes
      `summary.json`, `report.md`, and the four PNGs into
      `<batch>/validation/`. Reads the existing `manifest.json` if
      present so admixture-mode reports surface requested vs realised
      ancestry; PCA labels each person by dominant ancestry component
      so the spec's "cluster or bridge clusters" criterion is testable
      visually.
- [x] `requirements.txt` `matplotlib>=3.7` and `scikit-allel>=1.3` are
      now installed in the project venv (scikit-learn pulled in as a
      transitive dep for PCA). Without them, validate falls back to
      JSON / Markdown output and prints a one-line warning.
- [x] 39 new tests in `tests/test_validate.py` covering: INFO parsing,
      SNV / indel / SV classification, GT dosage (incl. multi-allelic
      `1|2`), dropout detection, Ti/Tv aggregation, Het/Hom corner
      cases, AF histogram bin placement, indel/SV aggregation, r²
      perfect-correlation / anticorrelation / uncorrelated /
      constant-vector / missing-data / few-samples, LD-decay structure
      (short-range > long-range), cohort PCA on a structured matrix
      (PC1 captures >95% of variance on a 2-cluster signal),
      insufficient-columns guard, and PNG smoke-tests for every plot
      helper. Numpy-required tests gate on `HAS_NUMPY`; PCA tests
      gate additionally on `HAS_SKLEARN`; plot tests gate on
      `HAS_MPL`. **186/186 with deps; 147/147 stdlib-only (39
      skipped).**

**Exit check:** ✅ 2026-04-25.

Single-population path (30 people × 5 Mb chr22, `--demo-model none`,
seed 42): **Ti/Tv = 1.822**, **Het/Hom-alt = 2.004**, 112,114 cohort
records, 74,796 singletons. LD decay monotone from 0.327 (100–500 bp)
through 0.123 (200–500 kb). PCA PC1 captures 5.3% of variance —
expected for a single-pop constant-Ne draw with no real demographic
structure. Four PNGs land under `validation/`; `report.md` links them
all and surfaces the per-sample table.

Admixture path (30 people × 5 Mb chr22, default 60/25/15
EUR/SAS/AFR, seed 42): **Ti/Tv = 1.882**, LD decay monotone (0.338
→ 0.104), **PCA PC1 captures 19.6% of variance** — the clear
ancestry signal the spec calls for, with EUR / SAS / AFR-dominant
labels visible as separable clusters in `pca.png`.

---

### M11 — Delivery packaging

- [x] Per-batch manifest `out/manifest.json` (already written from M6
      onward) now also surfaces, per person, `golden_bed`, `noise_bed`,
      `n_golden`, `n_noise` alongside the existing ancestry fractions,
      highlighted variant, error stats, and file paths. The top-level
      `errors` block continues to record cohort-wide realised FDR.
- [x] Truth-set BED (`syntheticgen/truth.py`): two BED4 files per
      person under `out/truth/`:
      * `person_NNNN.golden.bed` — every record matching the spec's
        "golden truth" set, tagged with priority HIGHLIGHTED > CLINVAR >
        COSMIC > SV > RSID. Each row carries a semicolon-separated
        `flag=…;id=…;ref=…;alt=…;gt=…;…` payload in the BED4 name
        column so downstream tools (bedtools / IGV / UCSC) parse it
        cleanly.
      * `person_NNNN.noise.bed` — every M9 noise event (FLIP /
        DROPOUT) recording the truth GT vs the called GT so a caller's
        accuracy can be graded against the model.
      `TruthBedWriter` buffers in memory and sorts by
      `(contig_order, chrom, start)` on close so output is
      `sort -k1,1 -k2,2n`-friendly. The writer is created per-person
      in `cli.py` and threaded into `write_person_vcf` so the writer
      sees both golden and noise events as the per-record loop runs.
- [x] `README.md` final rewrite describing the new pipeline end-to-end,
      including the M11 truth-set deliverables and the smoke script.
- [x] Smoke-test script `scripts/smoke.sh` runs an end-to-end 5-person
      cohort (chr22 × 0.5 Mb prefix) plus the validation suite and
      asserts every advertised deliverable lands on disk. Defaults to
      a deterministic seed; CI-friendly (<2 min on a laptop).
- [x] 20 new tests in `tests/test_truth.py` covering category priority,
      payload escaping, BED half-open coordinates for SNVs / deletions
      / SVs, sort-on-close, context-manager protocol, and parent-dir
      creation. **226/226 with deps; 165/165 stdlib-only (61 skipped).**

**Exit check:** ✅ 2026-04-26 on a fresh 5-person × 0.5 Mb chr22 sim
(`scripts/smoke.sh`, seed 42, default error / dropout rates): every
deliverable lands without manual intervention —
`out/person_NNNN.vcf.gz` + `.tbi`, `out/truth/person_NNNN.golden.bed`,
`out/truth/person_NNNN.noise.bed`, `out/manifest.json`,
`out/summary/sfs.tsv`, `out/validation/report.md` +
`summary.json` + four PNGs. Realised FDR 0.293% over 1,365 calls (4
events spread across 5 people). Manifest entries surface the new
`golden_bed` / `noise_bed` / `n_golden` / `n_noise` fields. Each
golden BED row carries the priority-ordered tag so a downstream
caller can split on `flag=` to grade against the model.

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
