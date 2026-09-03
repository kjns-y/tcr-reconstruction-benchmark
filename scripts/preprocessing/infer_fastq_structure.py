#!/usr/bin/env python3
"""Inspect a paired FASTQ sample and infer the 10x barcode/UMI read structure."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path
from typing import Iterator, Optional, TextIO, Tuple


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read1", required=True, type=Path, help="R1 FASTQ or FASTQ.gz")
    parser.add_argument("--read2", required=True, type=Path, help="R2 FASTQ or FASTQ.gz")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV")
    parser.add_argument(
        "--known-barcodes",
        type=Path,
        help="Optional TSV containing a barcode column, used to verify barcode coordinates",
    )
    parser.add_argument(
        "--sample-reads",
        type=int,
        default=10000,
        help="Maximum read pairs to inspect (default: 10000)",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    """Open plain or gzip-compressed text."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def iter_fastq(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield FASTQ identifier and sequence, failing on malformed records."""
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip("\n\r")
            separator = handle.readline()
            quality = handle.readline().rstrip("\n\r")
            if not sequence or not separator or not quality:
                raise ValueError(f"Truncated FASTQ record in {path}")
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}: {header.rstrip()}")
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch in {path}: {header.rstrip()}")
            yield header[1:].split()[0], sequence.upper()


def normalize_read_id(read_id: str) -> str:
    """Remove the conventional /1 or /2 mate suffix for comparison."""
    return read_id[:-2] if read_id.endswith(("/1", "/2")) else read_id


def load_known_barcodes(path: Optional[Path]) -> set[str]:
    """Load normalized 10x barcodes from a TSV when available."""
    if path is None:
        return set()
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: missing or empty known-barcode TSV: {path}")
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "barcode" not in reader.fieldnames:
            raise ValueError(f"Known-barcode TSV has no barcode column: {path}")
        return {row["barcode"].removesuffix("-1") for row in reader if row.get("barcode")}


def main() -> None:
    """Measure pair consistency and report evidence for R1 barcode/UMI, R2 cDNA."""
    args = parse_args()
    if args.sample_reads <= 0:
        raise SystemExit("ERROR: --sample-reads must be positive")
    for path in (args.read1, args.read2):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing or empty FASTQ: {path}")

    lengths_r1: Counter[int] = Counter()
    lengths_r2: Counter[int] = Counter()
    barcode16: Counter[str] = Counter()
    umi10: Counter[str] = Counter()
    pairs = 0
    mismatched_ids = 0
    n_bases_r1 = 0
    n_bases_r2 = 0
    known_barcodes = load_known_barcodes(args.known_barcodes)
    barcode_window_hits: Counter[int] = Counter()
    tail_bases: Counter[str] = Counter()

    iterator1 = iter_fastq(args.read1)
    iterator2 = iter_fastq(args.read2)
    for _ in range(args.sample_reads):
        record1 = next(iterator1, None)
        record2 = next(iterator2, None)
        if record1 is None and record2 is None:
            break
        if record1 is None or record2 is None:
            raise ValueError("R1/R2 contain different numbers of records in the inspected range")
        id1, seq1 = record1
        id2, seq2 = record2
        pairs += 1
        mismatched_ids += normalize_read_id(id1) != normalize_read_id(id2)
        lengths_r1[len(seq1)] += 1
        lengths_r2[len(seq2)] += 1
        n_bases_r1 += seq1.count("N")
        n_bases_r2 += seq2.count("N")
        if len(seq1) >= 26:
            barcode16[seq1[:16]] += 1
            umi10[seq1[16:26]] += 1
            tail_bases.update(seq1[26:])
        if known_barcodes:
            for start in range(max(0, len(seq1) - 15)):
                barcode_window_hits[start] += seq1[start : start + 16] in known_barcodes

    if pairs == 0:
        raise ValueError("No read pairs were found")
    if mismatched_ids:
        raise ValueError(f"{mismatched_ids} inspected read pairs have mismatching identifiers")

    r1_min = min(lengths_r1)
    r1_max = max(lengths_r1)
    has_10x_core = r1_min >= 26 and r1_max <= 28
    long_r2 = min(lengths_r2) >= 50
    barcode_evidence = True
    if known_barcodes:
        best_start, best_hits = max(barcode_window_hits.items(), key=lambda item: item[1])
        second_hits = max((hits for start, hits in barcode_window_hits.items() if start != best_start), default=0)
        barcode_evidence = best_start == 0 and best_hits > max(10, 10 * second_hits)
    else:
        best_start, best_hits, second_hits = -1, 0, 0
    inference = (
        f"R1=16bp_cell_barcode+10bp_UMI+{r1_max - 26}bp_non_UMI_tail;R2=cDNA"
        if has_10x_core and long_r2 and barcode_evidence
        else "AMBIGUOUS"
    )

    rows = [
        ("read_pairs_inspected", pairs),
        ("mate_id_mismatches", mismatched_ids),
        ("r1_length_counts", ";".join(f"{length}:{count}" for length, count in sorted(lengths_r1.items()))),
        ("r2_length_counts", ";".join(f"{length}:{count}" for length, count in sorted(lengths_r2.items()))),
        ("r1_unique_first16", len(barcode16)),
        ("r1_unique_bases17_26", len(umi10)),
        ("known_barcodes", len(known_barcodes)),
        ("best_16bp_barcode_window_1based", best_start + 1 if best_start >= 0 else "NA"),
        ("best_barcode_window_hits", best_hits),
        ("second_best_barcode_window_hits", second_hits),
        ("r1_tail_after_base26", ";".join(f"{base}:{count}" for base, count in sorted(tail_bases.items())) or "none"),
        ("r1_n_bases", n_bases_r1),
        ("r2_n_bases", n_bases_r2),
        ("inference", inference),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("property\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")

    print(f"Inspected {pairs} paired reads; R1 lengths={dict(lengths_r1)}, R2 lengths={dict(lengths_r2)}")
    print(f"Mate identifiers: MATCH; structure inference: {inference}")
    if known_barcodes:
        print(
            f"Known-barcode coordinate check: bases {best_start + 1}-{best_start + 16} "
            f"matched {best_hits}/{pairs} reads; second-best window matched {second_hits}"
        )


if __name__ == "__main__":
    main()
