"""Matplotlib plotting helpers for the M10 validation suite.

Every function here is a thin wrapper that takes already-computed
numbers and writes a PNG. The validate module does the math and stays
matplotlib-free; this module is the only place that imports
matplotlib, so a caller without matplotlib still gets the structured
JSON / Markdown artefacts.
"""

from __future__ import annotations

import math
from pathlib import Path


def _setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")  # headless box safe
    import matplotlib.pyplot as plt
    return plt


def plot_ld_decay(bins, out_path: Path) -> Path:
    """LD decay r² vs distance, log-x. `bins` is the list of dicts
    returned by `validate.ld_decay`."""
    plt = _setup_matplotlib()

    centres = []
    means = []
    counts = []
    for b in bins:
        if math.isnan(b["mean_r2"]):
            continue
        centres.append(math.sqrt(b["low_kb"] * b["high_kb"]))
        means.append(b["mean_r2"])
        counts.append(b["n_pairs"])

    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    if centres:
        ax.plot(centres, means, marker="o", linewidth=1.5)
        for x, y, n in zip(centres, means, counts):
            ax.annotate(f"n={n}", (x, y), textcoords="offset points",
                        xytext=(4, 4), fontsize=7, alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Distance between SNPs (kb, log scale)")
    ax.set_ylabel(r"Mean $r^2$")
    ax.set_title("Linkage disequilibrium decay")
    ax.set_ylim(0, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_af_histogram(edges, counts, out_path: Path) -> Path:
    plt = _setup_matplotlib()

    centres = [(edges[i] + edges[i + 1]) / 2.0 for i in range(len(counts))]
    widths = [edges[i + 1] - edges[i] for i in range(len(counts))]

    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    ax.bar(centres, counts, width=widths, edgecolor="black", linewidth=0.4)
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Records")
    ax.set_title("Allele frequency distribution")
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_pca(transformed, labels, out_path: Path,
             explained=None) -> Path:
    """`transformed` shape `(n_samples, 2)`; `labels` is a per-sample
    iterable (sample name or ancestry tag)."""
    plt = _setup_matplotlib()

    fig, ax = plt.subplots(figsize=(6, 5), dpi=120)
    if transformed is None or len(transformed) == 0:
        ax.text(0.5, 0.5, "PCA unavailable", ha="center", va="center",
                transform=ax.transAxes)
    else:
        unique = list(dict.fromkeys(labels))  # preserve order
        cmap = plt.get_cmap("tab10")
        for i, lab in enumerate(unique):
            mask = [lab == L for L in labels]
            xs = [transformed[k][0] for k, m in enumerate(mask) if m]
            ys = [transformed[k][1] for k, m in enumerate(mask) if m]
            ax.scatter(xs, ys, label=str(lab), s=30, alpha=0.8,
                       color=cmap(i % 10))
        if len(unique) > 1:
            ax.legend(fontsize=8, frameon=True, loc="best")
    xl, yl = "PC1", "PC2"
    if explained is not None and len(explained) >= 2:
        xl = f"PC1 ({explained[0] * 100:.1f}%)"
        yl = f"PC2 ({explained[1] * 100:.1f}%)"
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title("Cohort PCA")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_indel_lengths(length_counts: dict, out_path: Path) -> Path:
    """`length_counts` keyed by signed length-delta (insertions positive,
    deletions negative)."""
    plt = _setup_matplotlib()

    if not length_counts:
        keys = [0]
        vals = [0]
    else:
        keys = sorted(length_counts)
        vals = [length_counts[k] for k in keys]

    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    colours = ["#d95f02" if k < 0 else "#1b9e77" for k in keys]
    ax.bar(keys, vals, color=colours, edgecolor="black", linewidth=0.4)
    ax.axvline(0, linestyle="--", alpha=0.3, color="black")
    ax.set_xlabel("Indel length (bp; - = deletion, + = insertion)")
    ax.set_ylabel("Records")
    ax.set_title("Indel length distribution")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
