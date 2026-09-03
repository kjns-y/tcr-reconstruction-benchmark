#!/usr/bin/env python3
"""Plot four fixed YASIM-scTCR depth design points for AsTCR performance."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COLORS = {"TRUST4": "#3569A8", "MiXCR": "#D58A24"}
MARKERS = {"TRA": "o", "TRB": "s"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "results/benchmark/simulation_metrics.tsv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/figures")
    return parser.parse_args()


def plot_measure(data: pd.DataFrame, measure: str, output: Path) -> None:
    subset = data[data.metric == "AsTCR"].copy()
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    for (method, chain), group in subset.groupby(["method", "chain"], sort=False):
        group = group.sort_values("depth")
        ax.plot(
            group.depth, group[measure], label=f"{method} {chain}",
            color=COLORS.get(method, "#777777"), marker=MARKERS.get(chain, "o"),
            linestyle="-" if chain == "TRA" else "--", linewidth=2, markersize=6,
        )
    ax.set_xscale("log")
    ax.set_xticks([2, 10, 50, 100], ["2×", "10×", "50×", "100×"])
    ax.set_ylim(0, 1)
    ax.set_xlabel("TCR sequencing depth (four fixed design points)")
    ax.set_ylabel(measure.capitalize())
    ax.set_title(f"YASIM-scTCR depth experiment: AsTCR {measure}")
    ax.grid(color="#D9DEE3", linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.input, sep="\t")
    required = {"method", "depth", "chain", "metric", "sensitivity", "accuracy"}
    if missing := required - set(data.columns):
        raise ValueError(f"Simulation metrics lack columns: {sorted(missing)}")
    data["depth"] = pd.to_numeric(data["depth"], errors="raise")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_measure(data, "sensitivity", args.output_dir / "simulation_depth_sensitivity.pdf")
    plot_measure(data, "accuracy", args.output_dir / "simulation_depth_accuracy.pdf")
    print(f"Saved depth curves in: {args.output_dir}")


if __name__ == "__main__":
    main()
