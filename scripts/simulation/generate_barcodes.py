#!/usr/bin/env python3
"""Generate deterministic A/C/G/T-only 16 bp cell barcodes."""

from __future__ import annotations

import argparse
from pathlib import Path

ALPHABET = "ACGT"


def encode(index: int, length: int = 16) -> str:
    chars = ["A"] * length
    value = index
    for position in range(length - 1, -1, -1):
        chars[position] = ALPHABET[value % 4]
        value //= 4
    return "".join(chars)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.cells <= 0 or args.cells > 4**16:
        raise ValueError("cells must be between 1 and 4^16")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(f"{encode(i)}\n" for i in range(args.cells)), encoding="utf-8")
    print(f"Saved {args.cells} deterministic barcodes: {args.output}")


if __name__ == "__main__":
    main()

