# Test-suite review

Status: **Phase A complete (2026-05-15)**. Phase B–D pending.

This document holds the plan for reviewing the synthetic_people test suite plus its findings as each phase completes. The goal is to assess test quality, identify redundancy / consolidation opportunities, and surface anything obsolete — without losing coverage.

## Why now

626 tests, 191 s wall-time, 36 test files. Big enough that low-quality patterns aren't drowned in noise any more (recent example: 30 redundant ClinVar downloads in a single CI run — PR #90). Small enough that a structured review can finish in a few PRs rather than open-ended.

---

## Plan

### Phase A — Cheap signal first (read-only, ~30 min) ✅

Per-file scorecard so judgments later have data to lean on. No code changes.

Outputs go to the **Phase A findings** section below.

### Phase B — Pattern detection (~1-2 hr) ⏳

From Phase A's scorecard, look for:
- **Redundancy clusters** — tests covering the same code path with trivial input variation, often spread across files because added one feature at a time.
- **Brittle / over-coupled tests** — assertions on implementation details vs. contracts. Telltale: long mock chains, asserting subprocess argv beyond the load-bearing flag.
- **Slow tests not earning their keep** — anything in the >2 s tier that isn't an integration smoke.
- **Skipped tests that shouldn't be skipped** — tests gated on optional deps that CI now has. Per the standing rule "skipped tests are not validation", install the dep + unskip.
- **Tests for deprecated/removed code** — leftover assertions about behaviours that no longer exist.

### Phase C — Deep dive per area (~2-3 hr, prioritized) ⏳

Per-area review. Phase A's scorecard ranks priority. For each area, ask:
1. What's the contract being tested? State it in one sentence.
2. Is each test pulling weight against that contract? Drop or merge the ones that aren't.
3. Does a parameterized form catch the same regressions with less code?
4. Are helpers duplicated?

### Phase D — Action, in scoped PRs ⏳

Never one big PR. Per-area changes so each is reviewable + reversible.

**Discipline:** before deleting or merging a test, verify the case it covers is still asserted somewhere. For each test removed, grep that the function / branch under test is still referenced by another test's assertions. If not, that test was load-bearing and the consolidation lost coverage.

### Out of scope (won't do unless a later phase argues for it)

- Mass coverage instrumentation (coverage.py with HTML output).
- Migrating `unittest.TestCase` → pytest-style functions.
- Parallelising tests (pytest-xdist) — `tests/_shared_cache.py` would need race-safety review first.

---

## Phase A — Findings

### Top-line numbers (2026-05-15)

- **626 tests passed**, 2 skipped, 24 subtests passed.
- **Wall-time: 191 s** (3 min 11 s on the dev machine; CI is comparable).
- **36 test files**, mean 17 tests / file (median 13, max 83).

### Per-file scorecard

Columns: tests = count of `def test_`; LOC = file line count; cli.main = `cli_module.main` invocations (slow integration); setUpClass = number of `setUpClass` hooks (fixture-sharing); subproc = `subprocess.run/Popen/check_*` count (external-tool integration); skipUnless = optional-dep gate count.

| Tests | LOC | cli.main | setUpClass | subproc | skipUnless | File |
|------:|-----:|--------:|----------:|--------:|----------:|------|
| 83 | 1089 | 0 | 0 | 5 | 5 | test_validate.py |
| 46 | 810 | 4 | 2 | 1 | 2 | test_cohort_arrow_cli.py |
| 43 | 806 | 0 | 0 | 0 | 8 | test_config.py |
| 36 | 970 | 3 | 0 | 2 | 7 | test_reference.py |
| 25 | 357 | 0 | 0 | 0 | 1 | test_chunked_simulation.py |
| 23 | 293 | 0 | 0 | 0 | 0 | test_overlays.py |
| 22 | 411 | 0 | 3 | 5 | 3 | test_cohort_parallel_write.py |
| 22 | 207 | 0 | 0 | 0 | 0 | test_sv.py |
| 20 | 353 | 0 | 0 | 2 | 1 | test_mutation_spectrum.py |
| 20 | 202 | 0 | 0 | 0 | 0 | test_truth.py |
| 19 | 458 | 0 | 0 | 0 | 7 | test_cohort_arrow.py |
| 18 | 200 | 0 | 0 | 0 | 0 | test_quality.py |
| 18 | 188 | 0 | 0 | 0 | 0 | test_errors.py |
| 17 | 330 | 0 | 0 | 0 | 0 | test_carriers_sidecar.py |
| 15 | 321 | 0 | 4 | 0 | 6 | test_cohort_derivation.py |
| 14 | 474 | 0 | 0 | 0 | 0 | test_overlay_planners.py |
| 14 | 113 | 0 | 0 | 0 | 0 | test_titv.py |
| 13 | 222 | 0 | 0 | 0 | 4 | test_admixture.py |
| 13 | 81 | 0 | 0 | 0 | 0 | test_cli_chromosomes.py |
| 12 | 593 | 0 | 0 | 0 | 6 | test_streaming_cohort.py |
| 12 | 190 | 0 | 0 | 0 | 0 | test_cohort.py |
| 12 | 119 | 0 | 0 | 0 | 0 | test_sfs.py |
| 11 | 382 | 0 | 0 | 5 | 12 | test_cohort_arrow_bcf_writer.py |
| 10 | 272 | 3 | 3 | 2 | 3 | test_cli_modes.py |
| 10 | 158 | 0 | 0 | 0 | 2 | test_coalescent.py |
| 10 | 130 | 0 | 0 | 0 | 0 | test_sample_ids.py |
| 9 | 259 | 3 | 1 | 0 | 1 | test_resume.py |
| 9 | 212 | 0 | 2 | 4 | 3 | test_bcf_writer.py |
| 9 | 198 | 0 | 0 | 2 | 4 | test_phase1_concurrency.py |
| 7 | 526 | 1 | 1 | 1 | 7 | test_performance_budgets.py |
| 7 | 320 | 3 | 3 | 6 | 4 | test_cohort_streaming.py |
| 7 | 232 | 0 | 0 | 0 | 5 | test_memprofile.py |
| 7 | 218 | 0 | 0 | 0 | 0 | test_phase2_prefetch.py |
| 7 | 68 | 0 | 0 | 0 | 0 | test_progress_formatting.py |
| 5 | 93 | 0 | 0 | 0 | 0 | test_multiallelic.py |
| 3 | 113 | 0 | 1 | 2 | 1 | test_overlay_loaders.py |

### Slowest 16 tests = 84% of wall-time

All are `cli.main` integration tests. Below this tier, durations drop to <1 s and the long tail of ~610 fast unit tests collectively takes ~31 s.

| Time | Phase | Test |
|-----:|-------|------|
| 17.59 s | setup | test_cohort_arrow_cli :: `CohortModeArrowStreamingParityTest::test_byte_identical_cohort_bcf` |
| 17.35 s | setup | test_resume :: `ResumeEndToEndTest::test_chr21_was_regenerated` |
| 17.28 s | setup | test_cohort_arrow_cli :: `CohortModeArrowParityTest::test_arrow_scratch_cleaned_up_on_success` |
| 17.09 s | setup | test_cohort_streaming :: `CohortStreamedDeterminismTest::test_per_record_equivalence_across_runs` |
| 13.30 s | call | test_performance_budgets :: `test_arrow_streaming_within_budget` |
| 9.83 s | setup | test_cohort_streaming :: `CohortStreamedStdpopsimMutationModelTest` |
| 8.96 s | setup | test_cli_modes :: `CliModeCohortTest` |
| 8.93 s | call | test_performance_budgets :: `test_arrow_within_budget` |
| 8.89 s | call | test_reference :: `AdmixtureCliFastaTest::test_cli_admixture_emits_real_refs` |
| 8.86 s | setup | test_cli_modes :: `CliModePerPersonTest` |
| 8.85 s | setup | test_cohort_streaming :: `CohortStreamedMultiChromTest` |
| 8.79 s | call | test_performance_budgets :: `test_sites_list_within_budget` |
| 8.74 s | call | test_resume :: `test_no_resume_wipes_and_redoes_everything` |
| 8.70 s | setup | test_cli_modes :: `CliModeBothTest` |
| 8.59 s | call | test_reference :: `ReferenceEndToEndTest::test_emitted_ref_matches_fasta` |
| 8.44 s | call | test_reference :: `AdmixtureCliFastaTest::test_cli_admixture_rejects_missing_fasta_early` |

### Skipped tests

Both skipped tests are intentional platform / size gates, not deps-missing:

| Test | Reason |
|------|--------|
| `test_phase1_concurrency :: test_non_linux_caps_to_one` | Skipped on Linux (no cap to test there) |
| `test_sample_ids :: test_at_pool_size_works` | Skipped when pool > 2 M to keep the test fast |

No action.

### Duplicated helpers (consolidation candidates)

Verbatim or near-verbatim copies across files:

| Helper | Files | Notes |
|--------|-------|-------|
| `_site` | test_bcf_writer, test_cohort_parallel_write, test_cohort_derivation, test_cohort_arrow_bcf_writer, test_overlay_planners, test_overlays | 6 definitions; at least 3 are byte-identical first lines |
| `_bcf_data_md5` | test_cohort_parallel_write, test_cohort_arrow_bcf_writer, test_cohort_arrow_cli | All 3 copies declared as the "same helper as test_… ." per their own docstrings |
| `_bcf_sample_columns` | Same 3 files as `_bcf_data_md5` | Pair-bound with the helper above |
| `_write_fasta` | test_reference, test_mutation_spectrum | Identical signatures |
| `_clinvar_rec` | 2 files (TBD) | |

### Initial observations (carry into Phase B)

1. **cli.main integration tests dominate the runtime.** 16 of 626 tests account for ~84 % of wall-time. Any consolidation that doesn't touch those is a rounding error on CI time — but those tests are also the ones with the highest coverage value (real cli flows). Don't cut them aggressively; do look for cases where the **same fixture** could serve multiple test classes (e.g. the three `test_cli_modes` classes each do their own `cli.main` setUpClass).
2. **`test_validate.py` is an outlier** at 83 tests / 1089 LOC and uses no `setUpClass` / no integration. That's probably the right shape for a validator-internals file, but worth confirming in Phase C whether the 83 tests are all pulling weight or whether some are redundant overlays.
3. **`test_config.py` is the second-largest file (43 tests, 806 LOC) with 8 `skipUnless` gates and zero integration**. Config logic shouldn't need that many optional-dep gates — Phase C should check whether the gates are still necessary or are vestigial.
4. **`test_reference.py` has 7 `skipUnless` gates for pysam** (recently grew rapidly through M12 work). pysam is now a hard requirement (PR #87), so several of these gates may be vestigial. Easy Phase B win.
5. **`_bcf_data_md5` + `_bcf_sample_columns` in 3 cohort-parity files** is the cleanest consolidation candidate — same helpers, same purpose, three copies. Hoist to `tests/_helpers.py` is a small first PR.
6. **No `setUpClass` in 26 of 36 files** — fixture-per-test pattern dominates. For the fast tests that's fine; the integration-heavy files already use `setUpClass` correctly.

---

## Phase B — Pattern detection ⏳

Pending. Will fill in:
- Redundancy clusters (specific test groups)
- Brittle-test sites
- Re-evaluation of `skipUnless` gates (which are still needed)
- Tests for code that no longer exists

## Phase C — Deep dive ⏳

Pending. Per-area starting prioritization (subject to revision after Phase B):
1. `test_cohort_arrow_*` family (3 files, ~79 tests total, overlapping parity matrices)
2. `test_cohort_streaming.py` + `test_cli_modes.py` (3 cli.main classes each, similar fixtures)
3. `test_validate.py` (83 tests in one file — internal redundancy check)
4. `test_reference.py` (recently grown; old skipUnless gates)
5. `test_config.py` (43 tests, optional-dep gates)

## Phase D — Action PRs ⏳

Pending. Likely shape (rough sketch, will firm up after Phase B/C):
- PR: hoist `_bcf_data_md5`, `_bcf_sample_columns`, `_write_fasta`, `_site` to `tests/_helpers.py`
- PR: remove vestigial pysam `skipUnless` gates (pysam now hard-required since PR #87)
- PR per per-area consolidation (e.g. cohort parity matrix → parameterized)
- PR for any unskip-and-install-dep-in-CI work
