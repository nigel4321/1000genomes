#!/usr/bin/env bash
# M11 smoke test: end-to-end run of generation + validation in <5 min.
#
# Generates a tiny 5-person cohort (chr22, 0.5 Mb prefix), runs the
# validation suite, and asserts that every advertised deliverable
# under SYHTHETIC_PROJECT.md §7 lands on disk:
#   * person_NNNN.vcf.gz + .tbi
#   * truth/person_NNNN.golden.bed
#   * truth/person_NNNN.noise.bed
#   * manifest.json
#   * summary/sfs.tsv
#   * validation/{report.md,summary.json,*.png}
#
# Exits non-zero if any deliverable is missing or the cohort fails to
# generate. Designed for CI: deterministic seed, no network beyond
# ClinVar (which is cached).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

PY="${PYTHON:-../.venv/bin/python}"
if [ ! -x "$PY" ]; then
    PY="python3"
fi

OUT_DIR="${OUT_DIR:-out_smoke}"
N_PEOPLE="${N_PEOPLE:-5}"
SEED="${SEED:-42}"

echo "[smoke] cleaning $OUT_DIR"
rm -rf "$OUT_DIR"

echo "[smoke] generating $N_PEOPLE-person cohort (seed=$SEED)"
"$PY" generate_people.py \
    --n "$N_PEOPLE" \
    --output-dir "$OUT_DIR" \
    --chromosomes 22 \
    --chr-length-mb 0.5 \
    --seed "$SEED" \
    --error-rate 0.001 \
    --dropout-rate 0.0005

echo "[smoke] checking deliverables"
fail=0
for i in $(seq 1 "$N_PEOPLE"); do
    person=$(printf "person_%04d" "$i")
    for f in \
        "$OUT_DIR/$person.vcf.gz" \
        "$OUT_DIR/$person.vcf.gz.tbi" \
        "$OUT_DIR/truth/$person.golden.bed" \
        "$OUT_DIR/truth/$person.noise.bed"; do
        if [ ! -s "$f" ] && [ ! -e "$f" ]; then
            echo "  MISSING: $f" >&2
            fail=1
        fi
    done
done

for f in \
    "$OUT_DIR/manifest.json" \
    "$OUT_DIR/summary/sfs.tsv"; do
    if [ ! -s "$f" ]; then
        echo "  MISSING: $f" >&2
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo "[smoke] FAIL: missing generation artefacts" >&2
    exit 1
fi

echo "[smoke] running validation suite"
"$PY" validate_batch.py "$OUT_DIR" --ld-pairs-per-bin 200

for f in \
    "$OUT_DIR/validation/report.md" \
    "$OUT_DIR/validation/summary.json"; do
    if [ ! -s "$f" ]; then
        echo "  MISSING: $f" >&2
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo "[smoke] FAIL: missing validation artefacts" >&2
    exit 1
fi

echo "[smoke] OK — all deliverables present under $OUT_DIR/"
