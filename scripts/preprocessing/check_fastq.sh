#!/usr/bin/env bash
# Measure BC09 FASTQs before assigning barcode/UMI and cDNA roles.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
FASTQ_DIR="${1:-$ROOT/data/raw/BC09_GEX}"
OUT="$ROOT/data/metadata/BC09_fastq_stats.tsv"
LOG="$ROOT/logs/check_fastq.log"
require_command seqkit
mapfile -t files < <(find "$FASTQ_DIR" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) | sort)
(( ${#files[@]} > 0 )) || { echo "ERROR: no FASTQ files found in $FASTQ_DIR" >&2; exit 2; }

seqkit stats --tabular --all "${files[@]}" > "$OUT" 2> "$LOG"

python3 - "$OUT" <<'PY' | tee -a "$LOG"
from __future__ import annotations
import csv, re, sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows:
    raise SystemExit("ERROR: seqkit returned no rows")

def number(row, names):
    for name in names:
        if name in row:
            return float(row[name].replace(",", ""))
    raise KeyError(f"Could not find any of {names} in seqkit columns: {list(row)}")

pairs = {}
for row in rows:
    name = Path(row.get("file", row.get("file name", ""))).name
    match = re.search(r"(.+?)(?:_S\d+)?_L\d{3}_R([12])_\d{3}|(.+?)_([12])(?:\.f(?:ast)?q)", name)
    if not match:
        continue
    stem = match.group(1) or match.group(3)
    read = match.group(2) or match.group(4)
    pairs.setdefault(stem, {})[read] = row

print("FASTQ role inference based on measured lengths (not filenames alone):")
for stem, pair in sorted(pairs.items()):
    if set(pair) != {"1", "2"}:
        print(f"WARNING {stem}: missing mate; found reads {sorted(pair)}")
        continue
    n1 = int(number(pair["1"], ["num_seqs", "num sequences"]))
    n2 = int(number(pair["2"], ["num_seqs", "num sequences"]))
    a1 = number(pair["1"], ["avg_len", "avg length"])
    a2 = number(pair["2"], ["avg_len", "avg length"])
    status = "MATCH" if n1 == n2 else "MISMATCH"
    if a1 <= 40 and a2 >= 50:
        roles = "R1=barcode/UMI candidate; R2=cDNA candidate"
    elif a2 <= 40 and a1 >= 50:
        roles = "R2=barcode/UMI candidate; R1=cDNA candidate"
    else:
        roles = "AMBIGUOUS: inspect library documentation and read content"
    print(f"{stem}: reads {n1}/{n2} {status}; mean lengths {a1:g}/{a2:g}; {roles}")
    if n1 != n2:
        raise SystemExit(f"ERROR: paired read counts differ for {stem}")
PY

echo "Saved: $OUT" | tee -a "$LOG"
