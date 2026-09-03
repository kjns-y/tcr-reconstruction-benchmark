#!/usr/bin/env python3
"""Join YASIM receptor stats to cell barcodes and emit common truth format."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/evaluation"))
from normalization import normalize_barcode, normalize_cdr3aa, normalize_gene  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcr-stats", type=Path, required=True)
    parser.add_argument("--barcode-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "simulation/truth/simulation_ground_truth.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stats = pd.read_csv(args.tcr_stats, sep="\t", dtype=str, keep_default_na=False, quotechar="'")
    mapping = pd.read_csv(args.barcode_map, sep="\t", dtype=str, keep_default_na=False, quotechar="'")
    required_stats = {"UUID", "TRAV", "TRAJ", "TRBV", "TRBJ", "ACDR3_AA", "BCDR3_AA"}
    if missing := required_stats - set(stats.columns):
        raise ValueError(f"YASIM stats lacks: {sorted(missing)}")
    if missing := {"UUID", "barcode"} - set(mapping.columns):
        raise ValueError(f"YASIM barcode map lacks: {sorted(missing)}")
    joined = mapping.merge(stats, on="UUID", how="left", validate="many_to_one")
    if joined["TRAV"].isna().any():
        raise ValueError("Some barcode UUIDs did not join to the receptor truth.")
    rows = []
    for record in joined.to_dict(orient="records"):
        barcode = normalize_barcode(record["barcode"])
        rows.extend(
            [
                {"barcode": barcode, "chain": "TRA", "v_gene": normalize_gene(record["TRAV"]),
                 "j_gene": normalize_gene(record["TRAJ"]), "cdr3aa": normalize_cdr3aa(record["ACDR3_AA"])},
                {"barcode": barcode, "chain": "TRB", "v_gene": normalize_gene(record["TRBV"]),
                 "j_gene": normalize_gene(record["TRBJ"]), "cdr3aa": normalize_cdr3aa(record["BCDR3_AA"])},
            ]
        )
    output = pd.DataFrame(rows).drop_duplicates().sort_values(["barcode", "chain"])
    if output["barcode"].nunique() != mapping["barcode"].nunique():
        raise AssertionError("Simulation truth lost barcodes during conversion.")
    if not output["barcode"].str.fullmatch(r"[ACGT]{16}").all():
        raise ValueError("Simulation truth contains malformed or quoted barcodes.")
    if not output["v_gene"].str.match(r"^TR[AB]V").all() or not output["j_gene"].str.match(r"^TR[AB]J").all():
        raise ValueError("Simulation truth contains malformed V/J gene names.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(output.groupby("chain")["barcode"].nunique().to_string())
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
