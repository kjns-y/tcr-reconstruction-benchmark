#!/usr/bin/env python3
"""Plot the minimal BC09 Figure-2-style real-data benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
METRICS = ["CDR3", "V", "J", "AsTCR"]
METHOD_COLORS = {"TRUST4": "#3569A8", "MiXCR": "#D58A24"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, nargs="*", default=[
        ROOT / "results/benchmark/BC09_TRUST4_metrics.tsv",
        ROOT / "results/benchmark/BC09_MiXCR_metrics.tsv",
    ])
    parser.add_argument("--output-prefix", type=Path, default=ROOT / "results/figures/BC09_realdata_benchmark")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frames = [pd.read_csv(path, sep="\t") for path in args.metrics if path.is_file()]
    if not frames:
        raise FileNotFoundError("No benchmark metrics files are available.")
    data = pd.concat(frames, ignore_index=True)
    required = {"method", "chain", "metric", "accuracy", "sensitivity"}
    if missing := required - set(data.columns):
        raise ValueError(f"Metrics tables lack columns: {sorted(missing)}")
    methods = [m for m in ["TRUST4", "MiXCR"] if m in set(data.method)]
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.5), sharey=True, constrained_layout=True)
    for col, metric in enumerate(METRICS):
        for row, measure in enumerate(["accuracy", "sensitivity"]):
            ax = axes[row, col]
            subset = data[data.metric == metric]
            positions, labels, values, colors, hatches = [], [], [], [], []
            cursor = 0
            for chain in ["TRA", "TRB"]:
                for method in methods:
                    match = subset[(subset.chain == chain) & (subset.method == method)]
                    if match.empty:
                        continue
                    positions.append(cursor)
                    labels.append(f"{chain}\n{method}")
                    values.append(float(match.iloc[0][measure]))
                    colors.append(METHOD_COLORS.get(method, "#777777"))
                    hatches.append("" if chain == "TRA" else "//")
                    cursor += 1
                cursor += 0.5
            bars = ax.bar(positions, values, color=colors, edgecolor="#263238", linewidth=0.6)
            for bar, hatch in zip(bars, hatches):
                bar.set_hatch(hatch)
            ax.set_xticks(positions, labels, rotation=45, ha="right", fontsize=7)
            ax.set_ylim(0, 1)
            ax.grid(axis="y", color="#D9DEE3", linewidth=0.6)
            ax.set_axisbelow(True)
            ax.set_title(metric)
            if col == 0:
                ax.set_ylabel(measure.capitalize())
    fig.suptitle("BC09_TUMOR1: cell-level TCR reconstruction benchmark", fontsize=14)
    method_handles = [
        Patch(facecolor=METHOD_COLORS[method], edgecolor="#263238", label=method)
        for method in methods
    ]
    chain_handles = [
        Patch(facecolor="#E8ECEF", edgecolor="#263238", label="TRA"),
        Patch(facecolor="#E8ECEF", edgecolor="#263238", hatch="//", label="TRB"),
    ]
    fig.legend(
        handles=method_handles + chain_handles,
        loc="outside lower center",
        ncols=len(method_handles) + len(chain_handles),
        frameon=False,
        fontsize=8,
    )
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {args.output_prefix.with_suffix('.png')}")
    print(f"Saved: {args.output_prefix.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
