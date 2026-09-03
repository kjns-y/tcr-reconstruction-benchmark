#!/usr/bin/env bash
# Run TRUST4 in native barcode-aware mode after FASTQ QC confirms read roles.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
THREADS="${1:-8}"
FASTQ_DIR="${BC09_FASTQ_DIR:-$ROOT/data/raw/BC09_GEX}"
OUTDIR="${BC09_TRUST4_OUTPUT_DIR:-$ROOT/results/trust4/BC09}"
PREDICTION="${BC09_TRUST4_PREDICTION:-$ROOT/data/processed/BC09_trust4_prediction.tsv}"
COMMAND_LOG="${BC09_TRUST4_COMMAND_LOG:-$ROOT/logs/trust4_bc09_command.txt}"
RUN_LOG="${BC09_TRUST4_RUN_LOG:-$ROOT/logs/trust4_bc09.log}"
QC="${BC09_FASTQ_QC:-$ROOT/data/metadata/BC09_fastq_stats.tsv}"
OUTPUT_PREFIX="${BC09_TRUST4_PREFIX:-TRUST4_BC09}"
require_command run-trust4
[[ -s "$QC" ]] || { echo "ERROR: run check_fastq.sh first; missing $QC" >&2; exit 2; }

find_reference() {
    local name="$1" candidate
    for candidate in \
        "$ROOT/repo/TRUST4/$name" \
        "${TRUST4_REFERENCE_DIR:-}/$name" \
        "${CONDA_PREFIX:-}/share/trust4/$name" \
        "${CONDA_PREFIX:-}/share/trust4-"*/"$name" \
        "$(dirname "$(command -v run-trust4)")/$name" \
        "$(dirname "$(command -v run-trust4)")/../share/trust4/"*/"$name"; do
        [[ -f "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
    done
    return 1
}

imgt="$(find_reference 'human_IMGT+C.fa' || true)"
[[ -n "$imgt" ]] || {
    echo "ERROR: human_IMGT+C.fa not found. Set TRUST4_REFERENCE_DIR to the directory containing TRUST4 references." >&2
    exit 3
}
mapfile -t r1 < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_1.fastq.gz' | sort)
mapfile -t r2 < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_2.fastq.gz' | sort)
(( ${#r1[@]} > 0 && ${#r1[@]} == ${#r2[@]} )) || { echo "ERROR: paired FASTQ sets are absent or unequal." >&2; exit 4; }

# For this SRA layout the actual QC must show a short R1 and longer R2. The
# barcode/UMI coordinates are then the measured 10x 16+10 layout.
python3 - "$QC" <<'PY'
import csv, sys
rows=list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
def avg(r):
    for k in ('avg_len','avg length'):
        if k in r: return float(r[k].replace(',',''))
    raise SystemExit('ERROR: seqkit average-length column not found')
r1=[avg(r) for r in rows if '_1.fastq' in r.get('file','')]
r2=[avg(r) for r in rows if '_2.fastq' in r.get('file','')]
if not r1 or not r2 or max(r1) > 40 or min(r2) < 50:
    raise SystemExit('ERROR: QC does not support R1=barcode/UMI and R2=cDNA; inspect logs/check_fastq.log')
if min(r1) < 26:
    raise SystemExit('ERROR: R1 is shorter than the required 16 bp barcode + 10 bp UMI')
PY

mkdir -p "$OUTDIR"
cmd=(run-trust4 -f "$imgt" --ref "$imgt" -u "$FASTQ_DIR/*_2.fastq.gz" \
     --barcode "$FASTQ_DIR/*_1.fastq.gz" --UMI "$FASTQ_DIR/*_1.fastq.gz" \
     --readFormat 'bc:0:15,um:16:25' -t "$THREADS" -o "$OUTPUT_PREFIX" --od "$OUTDIR")
printf '%q ' "${cmd[@]}" > "$COMMAND_LOG"; printf '\n' >> "$COMMAND_LOG"
"${cmd[@]}" > "$RUN_LOG" 2>&1

report=""
for candidate in "$OUTDIR"/*_barcode_airr.tsv "$OUTDIR"/*_barcode_report.tsv; do
    [[ -s "$candidate" ]] && { report="$candidate"; break; }
done
[[ -n "$report" ]] || { echo "ERROR: TRUST4 finished without a barcode AIRR/report table; see $RUN_LOG" >&2; exit 5; }
python3 "$SCRIPT_DIR/parse_trust4.py" --input "$report" --output "$PREDICTION"
