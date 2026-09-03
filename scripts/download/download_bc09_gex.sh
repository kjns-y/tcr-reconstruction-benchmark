#!/usr/bin/env bash
# Download BC09 GEX SRA runs with explicit space checks and idempotency.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
THREADS="${1:-8}"
META="$ROOT/data/metadata/GSM3148575_accessions.tsv"
OUTDIR="$ROOT/data/raw/BC09_GEX"
SRA_DIR="$ROOT/data/raw/sra"
TMPDIR="$ROOT/data/raw/fasterq_tmp"
LOG="$ROOT/logs/download_bc09_gex.log"
SPACE_MULTIPLIER="${TCRBENCH_SPACE_MULTIPLIER:-8}"
SAFETY_GB="${TCRBENCH_SAFETY_GB:-20}"

for command_name in prefetch fasterq-dump awk df; do require_command "$command_name"; done
compressor="gzip"
command -v pigz >/dev/null 2>&1 && compressor="pigz"
[[ -s "$META" ]] || { echo "ERROR: missing metadata: $META. Run resolve_accessions.sh first." >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be a positive integer." >&2; exit 2; }
mkdir -p "$OUTDIR" "$SRA_DIR" "$TMPDIR" "$(dirname "$LOG")"

size_mb="$(awk -F '\t' 'NR>1 {s+=$7} END {printf "%.0f", s}' "$META")"
required_bytes="$((size_mb * 1024 * 1024 * SPACE_MULTIPLIER + SAFETY_GB * 1024 * 1024 * 1024))"
available_bytes="$(df -PB1 "$ROOT" | awk 'NR==2 {print $4}')"
{
    echo "[$(timestamp)] NCBI compressed estimate: ${size_mb} MB"
    echo "Required free space estimate: $required_bytes bytes (${SPACE_MULTIPLIER}x + ${SAFETY_GB} GiB safety)"
    echo "Available: $available_bytes bytes"
} | tee -a "$LOG"
if (( available_bytes < required_bytes )); then
    echo "ERROR: insufficient safe free space; download not started." | tee -a "$LOG" >&2
    exit 3
fi

tail -n +2 "$META" | cut -f4 | while IFS= read -r run; do
    [[ "$run" =~ ^SRR[0-9]+$ ]] || { echo "ERROR: invalid run in metadata: $run" >&2; exit 4; }
    r1="$OUTDIR/${run}_1.fastq.gz"
    r2="$OUTDIR/${run}_2.fastq.gz"
    if [[ -s "$r1" && -s "$r2" ]]; then
        echo "[$(timestamp)] SKIP complete: $run" | tee -a "$LOG"
        continue
    fi
    if [[ -e "$r1" || -e "$r2" ]]; then
        echo "ERROR: partial compressed output exists for $run; inspect it before retrying." | tee -a "$LOG" >&2
        exit 5
    fi
    echo "[$(timestamp)] START $run" | tee -a "$LOG"
    prefetch --max-size u --output-directory "$SRA_DIR" "$run" >>"$LOG" 2>&1
    run_tmp="$TMPDIR/$run"
    mkdir -p "$run_tmp"
    fasterq-dump "$SRA_DIR/$run/$run.sra" --split-files --threads "$THREADS" --temp "$run_tmp" --outdir "$OUTDIR" >>"$LOG" 2>&1
    [[ -s "$OUTDIR/${run}_1.fastq" && -s "$OUTDIR/${run}_2.fastq" ]] || {
        echo "ERROR: expected paired FASTQ was not produced for $run" | tee -a "$LOG" >&2
        exit 6
    }
    if [[ "$compressor" == "pigz" ]]; then
        pigz -p "$THREADS" "$OUTDIR/${run}_1.fastq" "$OUTDIR/${run}_2.fastq" >>"$LOG" 2>&1
    else
        gzip "$OUTDIR/${run}_1.fastq" "$OUTDIR/${run}_2.fastq" >>"$LOG" 2>&1
    fi
    rmdir "$run_tmp" 2>/dev/null || true
    echo "[$(timestamp)] DONE $run" | tee -a "$LOG"
done

echo "Download complete. Run scripts/preprocessing/check_fastq.sh next." | tee -a "$LOG"
