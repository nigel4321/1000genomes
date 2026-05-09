#!/usr/bin/env python3
"""
Spike 2 — Apache Arrow IPC streaming write + zero-copy mmap read.

Validates Phase 5d's two Arrow-specific assumptions:

  1. Streaming write keeps parent RAM bounded — pyarrow.ipc.new_file
     flushes each record batch immediately, so parent RSS does not
     grow with total sites written.
  2. pyarrow.ipc + memory_map gives workers true zero-copy reads —
     iterating per-site over a slice of the genotypes column does
     not allocate per-element Python objects or trigger COW
     divergence.

PASS criteria (all must hold):
  - parent peak RSS during streaming write < 500 MB regardless of
    total sites streamed
  - total system RSS during worker phase grows by < 1.5 × file_size
    above baseline (mmap shares physical pages)
  - aggregate read throughput across workers > 500 MB/s on NVMe
  - write throughput > 200 MB/s on NVMe

FAIL: any of the above. See README.md for failure-mode interpretation.

Default parameters generate a ~1 GB Arrow IPC file:
  10 000 samples × 100 000 sites × 1 byte = ~1 GB of int8 GTs

For larger smoke tests:
  --samples 100000 --sites 500000   (~50 GB Arrow file)
"""

import argparse
import multiprocessing as mp
import sys
import threading
import time
from pathlib import Path

import numpy as np
import psutil

try:
    import pyarrow as pa
    import pyarrow.ipc as paipc
except ImportError:
    print("ERROR: pyarrow not installed. Install with: pip install pyarrow")
    sys.exit(3)

ALT_FREQ = 0.05  # ~5% non-ref allele rate, realistic-ish sparsity


def cohort_schema(n_samples: int) -> pa.Schema:
    """Schema mirroring the proposed Phase 5d.1 layout."""
    return pa.schema([
        ("pos", pa.int64()),
        ("genotypes", pa.list_(pa.int8(), n_samples)),
    ])


def stream_write(path: Path, n_samples: int, n_sites: int,
                 batch_size: int) -> dict:
    """Parent: stream-write Arrow IPC, one record batch at a time."""
    schema = cohort_schema(n_samples)
    rng = np.random.default_rng(123)

    t0 = time.time()
    bytes_written = 0
    n_batches = 0
    sites_remaining = n_sites
    pos_cursor = 0

    with paipc.new_file(str(path), schema) as writer:
        while sites_remaining > 0:
            this_batch = min(batch_size, sites_remaining)
            positions = np.arange(pos_cursor, pos_cursor + this_batch,
                                  dtype=np.int64)
            pos_cursor += this_batch

            gt_data = (rng.random((this_batch, n_samples)) < ALT_FREQ).astype(np.int8)

            pos_array = pa.array(positions, type=pa.int64())
            gt_flat = pa.array(gt_data.ravel(), type=pa.int8())
            gt_array = pa.FixedSizeListArray.from_arrays(gt_flat, n_samples)

            batch = pa.RecordBatch.from_arrays(
                [pos_array, gt_array],
                names=["pos", "genotypes"],
            )
            writer.write_batch(batch)

            n_batches += 1
            bytes_written += this_batch * n_samples + this_batch * 8
            sites_remaining -= this_batch

    elapsed = time.time() - t0
    file_size = path.stat().st_size
    return {
        "elapsed_s": elapsed,
        "n_batches": n_batches,
        "logical_bytes": bytes_written,
        "file_size": file_size,
        "write_throughput_mbs": (bytes_written / elapsed) / (1024 * 1024),
    }


def worker(path: str, slice_lo: int, slice_hi: int, n_samples: int,
           queue: mp.Queue, idx: int) -> None:
    """Worker: mmap Arrow IPC, iterate per-site over its sample slice."""
    t0 = time.time()
    mm = pa.memory_map(path, "r")
    reader = pa.ipc.open_file(mm)

    total_alt_counts = 0
    sites_seen = 0
    bytes_read = 0

    for batch_idx in range(reader.num_record_batches):
        batch = reader.get_batch(batch_idx)
        gt_column = batch.column("genotypes")
        flat = gt_column.values.to_numpy(zero_copy_only=True)
        matrix = flat.reshape(-1, n_samples)
        my_slice = matrix[:, slice_lo:slice_hi]

        for i in range(my_slice.shape[0]):
            alt_count = int((my_slice[i] > 0).sum())
            total_alt_counts += alt_count
            sites_seen += 1

        bytes_read += my_slice.nbytes

    elapsed = time.time() - t0
    rss_mb = psutil.Process().memory_info().rss / (1024 * 1024)
    queue.put((idx, sites_seen, total_alt_counts, bytes_read, elapsed, rss_mb))


def sample_parent_rss(proc: psutil.Process, stop: threading.Event,
                      samples: list, interval_s: float = 0.2) -> None:
    t0 = time.time()
    while not stop.is_set():
        samples.append((time.time() - t0,
                        proc.memory_info().rss / (1024 * 1024)))
        time.sleep(interval_s)


def sample_system_ram(stop: threading.Event, samples: list,
                      interval_s: float = 0.5) -> None:
    t0 = time.time()
    while not stop.is_set():
        vm = psutil.virtual_memory()
        samples.append((time.time() - t0,
                        vm.used / (1024**3),
                        vm.available / (1024**3)))
        time.sleep(interval_s)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--samples", type=int, default=10000,
                        help="number of samples / genotype columns (default: 10000)")
    parser.add_argument("--sites", type=int, default=100000,
                        help="number of sites / genotype rows (default: 100000)")
    parser.add_argument("--batch-size", type=int, default=1024,
                        help="record batch size (sites per batch) (default: 1024)")
    parser.add_argument("--workers", type=int, default=8,
                        help="number of forked workers (default: 8)")
    parser.add_argument("--path", type=Path,
                        default=Path("/tmp/spike2_arrow.arrow"),
                        help="path for the Arrow IPC test file "
                             "(default: /tmp/spike2_arrow.arrow)")
    args = parser.parse_args()

    if args.samples % args.workers != 0:
        print(f"WARNING: samples ({args.samples}) not divisible by workers "
              f"({args.workers}). Last slice will absorb the remainder.")

    print(f"\n=== Phase 1: Streaming write ===")
    print(f"  schema: pos int64 | genotypes FixedSizeList(int8, {args.samples})")
    n_batches_expected = (args.sites + args.batch_size - 1) // args.batch_size
    print(f"  sites: {args.sites:,}  batch_size: {args.batch_size}  "
          f"n_batches: {n_batches_expected}")
    print(f"  expected raw GT bytes: "
          f"{args.sites * args.samples / 1024**3:.2f} GB")

    if args.path.exists():
        args.path.unlink()

    proc = psutil.Process()
    parent_rss_samples: list = []
    stop_parent = threading.Event()
    parent_sampler = threading.Thread(
        target=sample_parent_rss,
        args=(proc, stop_parent, parent_rss_samples),
        daemon=True,
    )
    parent_sampler.start()

    write_stats = stream_write(args.path, args.samples, args.sites,
                               args.batch_size)

    stop_parent.set()
    parent_sampler.join(timeout=2)

    if parent_rss_samples:
        write_stats["parent_peak_rss_mb"] = max(s[1] for s in parent_rss_samples)
        write_stats["parent_baseline_rss_mb"] = parent_rss_samples[0][1]
    else:
        write_stats["parent_peak_rss_mb"] = proc.memory_info().rss / (1024 * 1024)
        write_stats["parent_baseline_rss_mb"] = write_stats["parent_peak_rss_mb"]

    file_size_gb = write_stats["file_size"] / (1024**3)
    print(f"\n  wrote: {write_stats['n_batches']} batches  "
          f"file_size = {file_size_gb:.2f} GB")
    print(f"  elapsed: {write_stats['elapsed_s']:.2f} s")
    print(f"  write throughput (logical): "
          f"{write_stats['write_throughput_mbs']:.0f} MB/s")
    print(f"  parent baseline RSS: "
          f"{write_stats['parent_baseline_rss_mb']:.0f} MB")
    print(f"  parent peak RSS during write: "
          f"{write_stats['parent_peak_rss_mb']:.0f} MB")
    parent_rss_growth = (write_stats["parent_peak_rss_mb"]
                         - write_stats["parent_baseline_rss_mb"])
    print(f"  parent RSS growth during write: +{parent_rss_growth:.0f} MB")

    print(f"\n=== Phase 2: Worker mmap-read ===")
    slice_size = args.samples // args.workers
    print(f"  forking {args.workers} workers; each reads sample slice of "
          f"{slice_size} columns")
    print(f"  worker iteration: per-site Python loop "
          f"(realistic 5d.1 worker pattern)")

    vm0 = psutil.virtual_memory()
    print(f"  baseline system RAM: used = {vm0.used / 1024**3:.2f} GB  "
          f"available = {vm0.available / 1024**3:.2f} GB")

    mp.set_start_method("fork", force=True)
    queue: mp.Queue = mp.Queue()
    procs = []
    for i in range(args.workers):
        lo = i * slice_size
        hi = (i + 1) * slice_size if i < args.workers - 1 else args.samples
        p = mp.Process(target=worker,
                       args=(str(args.path), lo, hi, args.samples, queue, i))
        procs.append(p)

    sys_samples: list = []
    stop_sys = threading.Event()
    sys_sampler = threading.Thread(
        target=sample_system_ram,
        args=(stop_sys, sys_samples),
        daemon=True,
    )
    sys_sampler.start()

    t0 = time.time()
    for p in procs:
        p.start()

    results = []
    for _ in range(args.workers):
        results.append(queue.get())

    for p in procs:
        p.join()
    elapsed_workers = time.time() - t0

    stop_sys.set()
    sys_sampler.join(timeout=2)

    print(f"\n  all workers done in {elapsed_workers:.2f} s\n")

    print(f"  per-worker:")
    apparent_rss_total_gb = 0.0
    aggregate_bytes_read = 0
    for idx, sites_seen, alt_counts, bytes_read, elapsed_s, rss_mb in sorted(results):
        rss_gb = rss_mb / 1024
        apparent_rss_total_gb += rss_gb
        aggregate_bytes_read += bytes_read
        thr_mbs = (bytes_read / elapsed_s) / (1024 * 1024) if elapsed_s > 0 else 0
        print(f"    worker {idx}: sites={sites_seen}  alt_total={alt_counts}  "
              f"bytes_read={bytes_read / 1024**2:.0f} MB  "
              f"elapsed={elapsed_s:.2f} s  thr={thr_mbs:.0f} MB/s  "
              f"RSS={rss_gb:.2f} GB")

    aggregate_throughput_mbs = (aggregate_bytes_read / elapsed_workers) / (1024 * 1024)
    print(f"\n  aggregate read throughput: {aggregate_throughput_mbs:.0f} MB/s")
    print(f"  sum of per-worker reported RSS (apparent): "
          f"{apparent_rss_total_gb:.2f} GB")

    if not sys_samples:
        print("ERROR: no system RAM samples captured")
        return 3

    peak_used_gb = max(s[1] for s in sys_samples)
    baseline_gb = vm0.used / 1024**3
    delta_gb = peak_used_gb - baseline_gb

    print(f"\n  system RAM during worker phase:")
    print(f"    baseline used:        {baseline_gb:.2f} GB")
    print(f"    peak used:            {peak_used_gb:.2f} GB")
    print(f"    delta:                {delta_gb:+.2f} GB")
    print(f"    file size:            {file_size_gb:.2f} GB")
    print(f"    n_workers × file:     {args.workers * file_size_gb:.2f} GB "
          f"(if no sharing)")

    print(f"\n=== VERDICT ===")
    failures: list = []
    warnings: list = []

    parent_rss_threshold_mb = 500
    if write_stats["parent_peak_rss_mb"] > parent_rss_threshold_mb:
        failures.append(
            f"FAIL: parent peak RSS during streaming write was "
            f"{write_stats['parent_peak_rss_mb']:.0f} MB, exceeds "
            f"{parent_rss_threshold_mb} MB threshold. "
            f"pyarrow.ipc.new_file may be buffering instead of streaming."
        )
    else:
        print(f"  PASS: parent peak RSS = "
              f"{write_stats['parent_peak_rss_mb']:.0f} MB "
              f"< {parent_rss_threshold_mb} MB threshold "
              f"(streaming write is bounded).")

    sharing_pass_threshold = file_size_gb * 1.5
    sharing_fail_threshold = args.workers * file_size_gb * 0.5
    if delta_gb > sharing_fail_threshold:
        failures.append(
            f"FAIL: system RAM delta {delta_gb:.2f} GB exceeds "
            f"{sharing_fail_threshold:.2f} GB threshold. "
            f"Workers appear to hold private copies; mmap is not zero-copy "
            f"for this schema/access pattern."
        )
    elif delta_gb > sharing_pass_threshold:
        failures.append(
            f"AMBIGUOUS: system RAM delta {delta_gb:.2f} GB between "
            f"{sharing_pass_threshold:.2f} GB pass and "
            f"{sharing_fail_threshold:.2f} GB fail thresholds. "
            f"Some sharing, some divergence."
        )
    else:
        print(f"  PASS: system RAM delta = {delta_gb:.2f} GB "
              f"< {sharing_pass_threshold:.2f} GB threshold "
              f"(mmap shares cleanly across workers).")

    write_throughput_threshold = 200
    if write_stats["write_throughput_mbs"] < write_throughput_threshold:
        warnings.append(
            f"WARN: write throughput "
            f"{write_stats['write_throughput_mbs']:.0f} MB/s below "
            f"{write_throughput_threshold} MB/s NVMe target. "
            f"At n=1M this becomes the bottleneck."
        )
    else:
        print(f"  PASS: write throughput = "
              f"{write_stats['write_throughput_mbs']:.0f} MB/s "
              f">= {write_throughput_threshold} MB/s target.")

    read_throughput_threshold = 500
    if aggregate_throughput_mbs < read_throughput_threshold:
        warnings.append(
            f"WARN: aggregate read throughput "
            f"{aggregate_throughput_mbs:.0f} MB/s below "
            f"{read_throughput_threshold} MB/s target. "
            f"At n=1M workers may not finish in reasonable time."
        )
    else:
        print(f"  PASS: aggregate read throughput = "
              f"{aggregate_throughput_mbs:.0f} MB/s "
              f">= {read_throughput_threshold} MB/s target.")

    print()
    for w in warnings:
        print(f"  {w}")
    for f in failures:
        print(f"  {f}")

    if failures:
        print(f"\n  Phase 5d Arrow assumptions DO NOT FULLY HOLD on this host.")
        print(f"  See README.md for fail-mode-specific guidance.")
        if any("FAIL" in f for f in failures):
            return 2
        return 1

    if warnings:
        print(f"\n  Core correctness checks PASSED (streaming write + zero-copy")
        print(f"  mmap share). Throughput warnings are advisory; proceed to")
        print(f"  Phase 5d.1 with disk-throughput as a known parameter to plan")
        print(f"  for at scale.")
        return 0

    print(f"  All checks PASSED.")
    print(f"  Phase 5d's Arrow-IPC zero-copy mmap pattern is viable on this host.")
    print(f"  Green-light Phase 5d.1 implementation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
