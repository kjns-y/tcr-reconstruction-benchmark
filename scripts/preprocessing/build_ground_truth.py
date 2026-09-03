#!/usr/bin/env python3
"""Build the BC09 evaluation truth table from official 10x VDJ annotations."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/evaluation"))
from normalization import normalize_barcode, normalize_cdr3aa, normalize_gene  # noqa: E402

REQUIRED_INPUT = {"barcode", "productive", "full_length", "chain", "v_gene", "j_gene", "cdr3"}
REQUIRED_OUTPUT = ["barcode", "chain", "v_gene", "j_gene", "cdr3aa"]


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    available = set(columns)
    return next((name for name in candidates if name in available), None)


def build_ground_truth(input_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(input_path, compression="infer", low_memory=False)
    missing = REQUIRED_INPUT - set(frame.columns)
    if missing:
        raise ValueError(f"Input is missing required 10x columns: {sorted(missing)}")
    keep = as_bool(frame["productive"]) & as_bool(frame["full_length"]) & frame["chain"].isin(["TRA", "TRB"])
    filtered = frame.loc[keep].copy()
    filtered["barcode"] = filtered["barcode"].map(normalize_barcode)
    filtered["chain"] = filtered["chain"].astype(str).str.upper()
    filtered["v_gene"] = filtered["v_gene"].map(normalize_gene)
    filtered["j_gene"] = filtered["j_gene"].map(normalize_gene)
    filtered["cdr3aa"] = filtered["cdr3"].map(normalize_cdr3aa)

    optional_map = {
        "cdr3_nt": first_existing(frame.columns, ["cdr3_nt"]),
        "contig_id": first_existing(frame.columns, ["contig_id"]),
        "umis": first_existing(frame.columns, ["umis"]),
        "reads": first_existing(frame.columns, ["reads"]),
        "raw_clonotype_id": first_existing(frame.columns, ["raw_clonotype_id", "clonotype_id"]),
    }
    columns = REQUIRED_OUTPUT.copy()
    for output_name, input_name in optional_map.items():
        if input_name:
            if output_name != input_name:
                filtered[output_name] = filtered[input_name]
            columns.append(output_name)
    filtered = filtered[columns]
    filtered = filtered[(filtered["barcode"] != "") & (filtered["cdr3aa"] != "")]
    filtered = filtered.drop_duplicates().sort_values(["barcode", "chain", "cdr3aa", "v_gene", "j_gene"])
    return filtered.reset_index(drop=True)


def validation_summary(frame: pd.DataFrame) -> dict[str, object]:
    per_cell = frame.groupby(["barcode", "chain"]).size()
    return {
        "number_of_cells": int(frame["barcode"].nunique()),
        "number_of_TRA_receptors": int((frame["chain"] == "TRA").sum()),
        "number_of_TRB_receptors": int((frame["chain"] == "TRB").sum()),
        "cells_with_TRA": int(frame.loc[frame["chain"] == "TRA", "barcode"].nunique()),
        "cells_with_TRB": int(frame.loc[frame["chain"] == "TRB", "barcode"].nunique()),
        "productive_receptors_per_cell_chain": {
            str(int(k)): int(v) for k, v in per_cell.value_counts().sort_index().items()
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=ROOT / "data/ground_truth/raw/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv.gz",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "data/ground_truth/BC09_ground_truth.tsv")
    parser.add_argument("--summary", type=Path, default=ROOT / "data/metadata/BC09_ground_truth_summary.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(f"Ground-truth source not found: {args.input}")
    result = build_ground_truth(args.input)
    if result.empty:
        raise RuntimeError("All rows were filtered out; inspect 10x boolean and chain columns.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    summary = validation_summary(result)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

