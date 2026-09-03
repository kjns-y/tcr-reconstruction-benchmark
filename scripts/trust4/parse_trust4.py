#!/usr/bin/env python3
"""Convert TRUST4 cell-level AIRR/barcode reports to the common receptor TSV."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/evaluation"))
from normalization import infer_chain, normalize_barcode, normalize_cdr3aa, normalize_gene, normalize_text  # noqa: E402


def parse_airr(frame: pd.DataFrame) -> List[Dict[str, str]]:
    barcode_col = next((c for c in ["cell_id", "barcode", "cell", "sequence_id"] if c in frame.columns), None)
    chain_col = next((c for c in ["locus", "chain"] if c in frame.columns), None)
    cdr3_col = next((c for c in ["junction_aa", "cdr3_aa", "CDR3aa"] if c in frame.columns), None)
    if not barcode_col or not cdr3_col or "v_call" not in frame or "j_call" not in frame:
        raise ValueError("AIRR table lacks cell/barcode, v_call, j_call, or junction_aa columns.")
    rows: List[Dict[str, str]] = []
    for row in frame.to_dict(orient="records"):
        chain = infer_chain(row.get(chain_col, "") if chain_col else "", row.get("v_call", ""), row.get("j_call", ""))
        if chain not in {"TRA", "TRB"}:
            continue
        rows.append(
            {
                "barcode": normalize_barcode(row[barcode_col]),
                "chain": chain,
                "v_gene": normalize_gene(str(row["v_call"]).split(",")[0]),
                "j_gene": normalize_gene(str(row["j_call"]).split(",")[0]),
                "cdr3aa": normalize_cdr3aa(row[cdr3_col]),
            }
        )
    return rows


def parse_chain_info(barcode: str, value: object) -> Dict[str, str] | None:
    text = normalize_text(value)
    if not text or text.lower() in {"none", "null"}:
        return None
    fields = next(csv.reader([text]))
    if len(fields) < 6:
        raise ValueError(f"Unexpected TRUST4 chain-info field for barcode {barcode}: {text}")
    v_gene, _d_gene, j_gene, c_gene, _cdr3nt, cdr3aa = fields[:6]
    chain = infer_chain(v_gene, j_gene, c_gene)
    if chain not in {"TRA", "TRB"}:
        return None
    return {
        "barcode": normalize_barcode(barcode),
        "chain": chain,
        "v_gene": normalize_gene(v_gene),
        "j_gene": normalize_gene(j_gene),
        "cdr3aa": normalize_cdr3aa(cdr3aa),
    }


def parse_barcode_report(frame: pd.DataFrame) -> List[Dict[str, str]]:
    barcode_col = next((c for c in frame.columns if c.lower() == "barcode"), frame.columns[0])
    ignored = {barcode_col, "cell_type", "cellType"}
    chain_columns = [c for c in frame.columns if c not in ignored]
    if not chain_columns:
        raise ValueError("TRUST4 barcode report has no chain-information columns.")
    rows: List[Dict[str, str]] = []
    for record in frame.to_dict(orient="records"):
        for column in chain_columns:
            parsed = parse_chain_info(str(record[barcode_col]), record.get(column, ""))
            if parsed:
                rows.append(parsed)
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/BC09_trust4_prediction.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    frame = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    is_airr = {"v_call", "j_call"}.issubset(frame.columns)
    rows = parse_airr(frame) if is_airr else parse_barcode_report(frame)
    output = pd.DataFrame(rows, columns=["barcode", "chain", "v_gene", "j_gene", "cdr3aa"])
    output = output[(output["barcode"] != "") & (output["cdr3aa"] != "")].drop_duplicates()
    output = output.sort_values(["barcode", "chain", "cdr3aa", "v_gene", "j_gene"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    counts = output.groupby("chain")["barcode"].nunique().to_dict()
    print(f"cells with reconstructed TRA: {counts.get('TRA', 0)}")
    print(f"cells with reconstructed TRB: {counts.get('TRB', 0)}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

