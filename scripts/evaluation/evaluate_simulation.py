#!/usr/bin/env python3
"""Evaluate every available method/depth prediction using the common evaluator."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from evaluate import evaluate_maps, load_receptors

ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, default=ROOT / "simulation/truth/simulation_ground_truth.tsv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/benchmark/simulation_metrics.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    truth = load_receptors(args.truth)
    frames: List[pd.DataFrame] = []
    for method, slug in [("TRUST4", "trust4"), ("MiXCR", "mixcr")]:
        for depth in [2, 10, 50, 100]:
            prediction_path = ROOT / f"data/processed/simulation_{depth}x_{slug}_prediction.tsv"
            if not prediction_path.is_file():
                print(f"SKIP missing: {prediction_path}")
                continue
            metrics, details = evaluate_maps(truth, load_receptors(prediction_path), method)
            metrics.insert(1, "depth", depth)
            frames.append(metrics)
            details_path = ROOT / f"results/benchmark/simulation_{depth}x_{method}_cell_details.tsv"
            details.to_csv(details_path, sep="\t", index=False)
    if not frames:
        raise FileNotFoundError("No simulation prediction tables are available.")
    output = pd.concat(frames, ignore_index=True)
    output = output[["method", "depth", "chain", "metric", "all_cells", "result_cells", "true_cells", "sensitivity", "accuracy"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False, float_format="%.6f")
    print(output.to_string(index=False))
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

