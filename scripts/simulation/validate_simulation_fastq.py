#!/usr/bin/env python3
"""Validate read counts, lengths, and pairing for all simulated depth FASTQs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

DEPTHS = (2, 10, 50, 100)
EXPECTED_CELLS = 500
RECEPTORS_PER_CELL = 2
NAME_RE = re.compile(r"BC09SIM_(2|10|50|100)x_([12])\.fastq\.gz$")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", type=Path, required=True, help="TSV produced by seqkit stats --all")
    return parser.parse_args()


def numeric(row: pd.Series, *names: str) -> int:
    """Read an integer from one of the supported seqkit column names."""
    for name in names:
        if name in row.index:
            return int(float(str(row[name]).replace(",", "")))
    raise KeyError(f"Missing seqkit column; expected one of {names}")


def main() -> None:
    """Fail unless every depth has the exact expected 10x-like read pair."""
    args = parse_args()
    if not args.stats.is_file() or args.stats.stat().st_size == 0:
        raise FileNotFoundError(args.stats)
    frame = pd.read_csv(args.stats, sep="\t")
    file_column = "file" if "file" in frame.columns else "file name"
    observed: Dict[Tuple[int, int], pd.Series] = {}
    for _, row in frame.iterrows():
        match = NAME_RE.search(Path(str(row[file_column])).name)
        if match:
            observed[(int(match.group(1)), int(match.group(2)))] = row

    for depth in DEPTHS:
        expected_reads = EXPECTED_CELLS * RECEPTORS_PER_CELL * depth
        for mate, expected_length in ((1, 26), (2, 150)):
            key = (depth, mate)
            if key not in observed:
                raise ValueError(f"Missing stats row for depth={depth}x R{mate}")
            row = observed[key]
            count = numeric(row, "num_seqs", "num sequences")
            minimum = numeric(row, "min_len", "min length")
            maximum = numeric(row, "max_len", "max length")
            if count != expected_reads:
                raise ValueError(
                    f"{depth}x R{mate}: expected {expected_reads} reads, observed {count}"
                )
            if minimum != expected_length or maximum != expected_length:
                raise ValueError(
                    f"{depth}x R{mate}: expected fixed {expected_length} bp, "
                    f"observed min/max {minimum}/{maximum}"
                )
        count_r1 = numeric(observed[(depth, 1)], "num_seqs", "num sequences")
        count_r2 = numeric(observed[(depth, 2)], "num_seqs", "num sequences")
        if count_r1 != count_r2:
            raise ValueError(f"{depth}x mate count mismatch: R1={count_r1}, R2={count_r2}")
        print(f"PASS {depth}x: {count_r1} pairs, R1=26 bp, R2=150 bp")


if __name__ == "__main__":
    main()
