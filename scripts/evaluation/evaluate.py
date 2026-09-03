#!/usr/bin/env python3
"""Compute paper-style cell-level TCR accuracy and sensitivity."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, DefaultDict, Dict, Iterable, List, Set, Tuple

import pandas as pd

from normalization import normalize_barcode, normalize_cdr3aa, normalize_gene

CHAINS = ("TRA", "TRB")
METRICS = ("CDR3", "V", "J", "AsTCR")
REQUIRED = {"barcode", "chain", "v_gene", "j_gene", "cdr3aa"}


@dataclass(frozen=True, order=True)
class Receptor:
    v_gene: str
    j_gene: str
    cdr3aa: str

    def display(self) -> str:
        return f"{self.v_gene}|{self.j_gene}|{self.cdr3aa}"


ReceptorMap = DefaultDict[str, DefaultDict[str, Set[Receptor]]]


def load_receptors(path: Path) -> ReceptorMap:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = REQUIRED - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    receptors: ReceptorMap = defaultdict(lambda: defaultdict(set))
    for row in frame.itertuples(index=False):
        barcode = normalize_barcode(getattr(row, "barcode"))
        chain = str(getattr(row, "chain")).strip().upper()
        if not barcode or chain not in CHAINS:
            continue
        receptor = Receptor(
            normalize_gene(getattr(row, "v_gene")),
            normalize_gene(getattr(row, "j_gene")),
            normalize_cdr3aa(getattr(row, "cdr3aa")),
        )
        receptors[barcode][chain].add(receptor)
    return receptors


def receptor_match(left: Receptor, right: Receptor, metric: str) -> bool:
    if metric == "CDR3":
        return bool(left.cdr3aa) and left.cdr3aa == right.cdr3aa
    if metric == "V":
        return bool(left.v_gene) and left.v_gene == right.v_gene
    if metric == "J":
        return bool(left.j_gene) and left.j_gene == right.j_gene
    if metric == "AsTCR":
        return (
            bool(left.cdr3aa and left.v_gene and left.j_gene)
            and left.cdr3aa == right.cdr3aa
            and left.v_gene == right.v_gene
            and left.j_gene == right.j_gene
        )
    raise ValueError(f"Unknown metric: {metric}")


def any_match(truth: Iterable[Receptor], prediction: Iterable[Receptor], metric: str) -> bool:
    return any(receptor_match(t, p, metric) for t in truth for p in prediction)


def evaluate_maps(truth: ReceptorMap, prediction: ReceptorMap, method: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: List[Dict[str, object]] = []
    detail_rows: List[Dict[str, object]] = []
    barcodes = sorted(set(truth) | set(prediction))
    for chain in CHAINS:
        for barcode in barcodes:
            truth_set = truth[barcode][chain]
            prediction_set = prediction[barcode][chain]
            flags = {metric: any_match(truth_set, prediction_set, metric) for metric in METRICS}
            detail_rows.append(
                {
                    "barcode": barcode,
                    "chain": chain,
                    "has_truth": bool(truth_set),
                    "has_prediction": bool(prediction_set),
                    "cdr3_correct": flags["CDR3"],
                    "v_correct": flags["V"],
                    "j_correct": flags["J"],
                    "astcr_correct": flags["AsTCR"],
                    "truth_receptors": ";".join(r.display() for r in sorted(truth_set)),
                    "predicted_receptors": ";".join(r.display() for r in sorted(prediction_set)),
                }
            )
        all_cells_set = {barcode for barcode in truth if truth[barcode][chain]}
        result_cells_set = {barcode for barcode in prediction if prediction[barcode][chain]}
        for metric in METRICS:
            true_cells_set = {
                barcode
                for barcode in result_cells_set
                if barcode in all_cells_set and any_match(truth[barcode][chain], prediction[barcode][chain], metric)
            }
            all_cells = len(all_cells_set)
            result_cells = len(result_cells_set)
            true_cells = len(true_cells_set)
            sensitivity = true_cells / all_cells if all_cells else 0.0
            accuracy = true_cells / result_cells if result_cells else 0.0
            if not (true_cells <= all_cells and true_cells <= result_cells):
                raise AssertionError("Cell counts violate benchmark set relationships.")
            if not (0.0 <= sensitivity <= 1.0 and 0.0 <= accuracy <= 1.0):
                raise AssertionError("A benchmark rate is outside [0, 1].")
            metric_rows.append(
                {
                    "method": method,
                    "chain": chain,
                    "metric": metric,
                    "all_cells": all_cells,
                    "result_cells": result_cells,
                    "true_cells": true_cells,
                    "sensitivity": sensitivity,
                    "accuracy": accuracy,
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(detail_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--details", type=Path, help="Default: METHOD output name with _cell_details.tsv")
    return parser.parse_args()


def default_details(output: Path) -> Path:
    name = output.name.replace("_metrics.tsv", "_cell_details.tsv")
    if name == output.name:
        name = output.stem + "_cell_details.tsv"
    return output.with_name(name)


def main() -> None:
    args = parse_args()
    for path in (args.truth, args.prediction):
        if not path.is_file():
            raise FileNotFoundError(path)
    metrics, details = evaluate_maps(load_receptors(args.truth), load_receptors(args.prediction), args.method)
    details_path = args.details or default_details(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    details_path.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output, sep="\t", index=False, float_format="%.6f")
    details.to_csv(details_path, sep="\t", index=False)
    print(metrics.to_string(index=False))
    print(f"Saved: {args.output}")
    print(f"Saved: {details_path}")


if __name__ == "__main__":
    main()

