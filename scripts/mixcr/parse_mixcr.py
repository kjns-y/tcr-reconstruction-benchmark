#!/usr/bin/env python3
"""Convert a MiXCR cell-split export to the common receptor TSV."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/evaluation"))
from normalization import infer_chain, normalize_barcode, normalize_cdr3aa, normalize_gene  # noqa: E402


def choose(columns: Iterable[str], candidates: Iterable[str], label: str) -> str:
    available = list(columns)
    lowered = {column.lower(): column for column in available}
    for candidate in candidates:
        if candidate in available:
            return candidate
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    raise ValueError(f"MiXCR export lacks {label}; available columns: {available}")


def clean_hit(value: object) -> str:
    text = str(value).split(",")[0]
    text = re.sub(r"\([^)]*\)$", "", text)
    return normalize_gene(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/BC09_mixcr_prediction.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    barcode_col = choose(frame.columns, ["cellId", "tagValueCELL", "barcode", "cell_id"], "cell barcode")
    v_col = choose(frame.columns, ["bestVGene", "vGene", "allVHitsWithScore", "allVHits"], "V gene")
    j_col = choose(frame.columns, ["bestJGene", "jGene", "allJHitsWithScore", "allJHits"], "J gene")
    cdr3_col = choose(frame.columns, ["aaSeqCDR3", "cdr3aa", "aaFeatureCDR3"], "CDR3 amino-acid sequence")
    chain_col = next((c for c in ["chain", "locus", "chains"] if c in frame.columns), None)
    output = pd.DataFrame()
    output["barcode"] = frame[barcode_col].map(normalize_barcode)
    output["v_gene"] = frame[v_col].map(clean_hit)
    output["j_gene"] = frame[j_col].map(clean_hit)
    output["cdr3aa"] = frame[cdr3_col].map(normalize_cdr3aa)
    output["chain"] = [
        infer_chain(frame.at[i, chain_col] if chain_col else "", output.at[i, "v_gene"], output.at[i, "j_gene"])
        for i in frame.index
    ]
    output = output[output["chain"].isin(["TRA", "TRB"]) & (output["barcode"] != "") & (output["cdr3aa"] != "")]
    output = output[["barcode", "chain", "v_gene", "j_gene", "cdr3aa"]].drop_duplicates()
    output = output.sort_values(["barcode", "chain", "cdr3aa", "v_gene", "j_gene"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    counts = output.groupby("chain")["barcode"].nunique().to_dict()
    print(f"cells with reconstructed TRA: {counts.get('TRA', 0)}")
    print(f"cells with reconstructed TRB: {counts.get('TRB', 0)}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

