#!/usr/bin/env python3
"""Wrap YASIM single-end cDNA reads in deterministic 10x-like R1/R2 pairs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
from pathlib import Path
from typing import Iterator, TextIO, Tuple

BARCODE_RE = re.compile(r"([ACGT]{16})_[AB]")


def open_text(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode + "t") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def fastq_records(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with open_text(path, "r") as handle:
        while True:
            record = tuple(handle.readline().rstrip("\n") for _ in range(4))
            if not record[0]:
                return
            if len(record) != 4 or not record[0].startswith("@") or not record[2].startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}: {record[0]}")
            yield record  # type: ignore[misc]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    args = parser.parse_args()
    args.r1.parent.mkdir(parents=True, exist_ok=True)
    args.r2.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(args.r1, "w") as r1_handle, open_text(args.r2, "w") as r2_handle:
        for header, sequence, plus, quality in fastq_records(args.input):
            match = BARCODE_RE.search(header)
            if not match:
                raise ValueError(f"Cannot recover 16 bp cell barcode from YASIM header: {header}")
            barcode = match.group(1)
            umi = "".join("ACGT"[byte % 4] for byte in hashlib.sha256(header.encode()).digest()[:10])
            r1_sequence = barcode + umi
            r1_handle.write(f"{header} 1:N:0:1\n{r1_sequence}\n+\n{'I' * len(r1_sequence)}\n")
            r2_handle.write(f"{header} 2:N:0:1\n{sequence}\n+\n{quality}\n")
            count += 1
    if count == 0:
        raise RuntimeError("No reads were converted.")
    print(f"Converted {count} YASIM reads to 10x-like R1/R2 pairs")


if __name__ == "__main__":
    main()

