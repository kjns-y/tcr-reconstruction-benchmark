#!/usr/bin/env python3
"""Select a few interpretable BC09 cells representing cross-method outcomes."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def read_details(path: Path, prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"barcode", "chain", "has_prediction", "astcr_correct", "truth_receptors", "predicted_receptors"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {sorted(missing)}")
    for column in ["has_prediction", "astcr_correct"]:
        frame[column] = frame[column].astype(str).str.lower().eq("true")
    return frame.rename(
        columns={
            "has_prediction": f"{prefix}_has_prediction",
            "astcr_correct": f"{prefix}_astcr_correct",
            "predicted_receptors": f"{prefix}_predicted_receptors",
        }
    )


def append_examples(
    rows: List[Dict[str, object]],
    joined: pd.DataFrame,
    category: str,
    mask: pd.Series,
    limit: int,
    include_both_chains: bool = False,
) -> None:
    """Append relevant rows for at most ``limit`` unique cell barcodes."""
    matching = joined.loc[mask].sort_values(["barcode", "chain"])
    selected_barcodes = matching["barcode"].drop_duplicates().head(limit)
    selected = joined[joined["barcode"].isin(selected_barcodes)] if include_both_chains else matching[
        matching["barcode"].isin(selected_barcodes)
    ]
    for record in selected.sort_values(["barcode", "chain"]).to_dict(orient="records"):
        rows.append({"category": category, **record})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trust4-details", type=Path, default=ROOT / "results/benchmark/BC09_TRUST4_cell_details.tsv")
    parser.add_argument("--mixcr-details", type=Path, default=ROOT / "results/benchmark/BC09_MiXCR_cell_details.tsv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/benchmark/example_cells.tsv")
    parser.add_argument("--max-per-category", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    trust = read_details(args.trust4_details, "trust4")
    mixcr = read_details(args.mixcr_details, "mixcr")
    shared = ["barcode", "chain", "has_truth", "truth_receptors"]
    joined = trust.merge(mixcr, on=shared, how="outer")
    boolean_columns = [
        "trust4_has_prediction",
        "trust4_astcr_correct",
        "mixcr_has_prediction",
        "mixcr_astcr_correct",
    ]
    for column in boolean_columns:
        joined[column] = joined[column].map(lambda value: bool(value) if pd.notna(value) else False).astype(bool)
    for column in ["trust4_predicted_receptors", "mixcr_predicted_receptors"]:
        joined[column] = joined[column].fillna("").astype(str)
    rows: List[Dict[str, object]] = []
    append_examples(rows, joined, "both_correct", joined.trust4_astcr_correct & joined.mixcr_astcr_correct, args.max_per_category)
    append_examples(
        rows, joined, "trust4_correct_mixcr_no_result",
        joined.trust4_astcr_correct & ~joined.mixcr_has_prediction, args.max_per_category,
    )
    append_examples(
        rows, joined, "mixcr_correct_trust4_no_result",
        joined.mixcr_astcr_correct & ~joined.trust4_has_prediction, args.max_per_category,
    )
    append_examples(
        rows, joined, "both_result_at_least_one_wrong",
        joined.trust4_has_prediction & joined.mixcr_has_prediction
        & ~(joined.trust4_astcr_correct & joined.mixcr_astcr_correct), args.max_per_category,
    )

    by_cell = joined.set_index(["barcode", "chain"])
    differing: Set[str] = set()
    for barcode in joined["barcode"].unique():
        subset = by_cell.loc[barcode] if barcode in by_cell.index.get_level_values(0) else pd.DataFrame()
        if isinstance(subset, pd.Series) or not {"TRA", "TRB"}.issubset(set(subset.index)):
            continue
        for prefix in ["trust4", "mixcr"]:
            if bool(subset.loc["TRA", f"{prefix}_astcr_correct"]) != bool(subset.loc["TRB", f"{prefix}_astcr_correct"]):
                differing.add(barcode)
                break
    append_examples(
        rows,
        joined,
        "TRA_TRB_performance_differs",
        joined.barcode.isin(sorted(differing)),
        args.max_per_category,
        include_both_chains=True,
    )
    columns = ["category"] + list(joined.columns)
    output = pd.DataFrame(rows, columns=columns).drop_duplicates(["category", "barcode", "chain"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(output.groupby("category").size().to_string() if not output.empty else "No examples available.")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
