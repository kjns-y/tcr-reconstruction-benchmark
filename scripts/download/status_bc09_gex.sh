#!/usr/bin/env bash
# Report the background BC09 download without modifying it.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
PID_FILE="$ROOT/logs/download_bc09_background.pid"
OUTDIR="$ROOT/data/raw/BC09_GEX"

if [[ -s "$PID_FILE" ]]; then
    job_pid="$(tr -d '[:space:]' < "$PID_FILE")"
    if [[ "$job_pid" =~ ^[1-9][0-9]*$ ]] && kill -0 "$job_pid" 2>/dev/null; then
        ps -o pid,ppid,etime,stat,cmd -p "$job_pid"
    else
        echo "Background PID is not active: ${job_pid:-invalid}"
    fi
else
    echo "No background PID file."
fi

printf 'completed_pairs\t'
find "$OUTDIR" -maxdepth 1 -type f -name 'SRR*_1.fastq.gz' 2>/dev/null \
    | sed 's/_1\.fastq\.gz$//' | while IFS= read -r prefix; do
        [[ -s "${prefix}_2.fastq.gz" ]] && basename "$prefix"
    done | wc -l
du -sh "$ROOT/data/raw/sra" "$OUTDIR" "$ROOT/data/raw/fasterq_tmp" 2>/dev/null || true
tail -n 20 "$ROOT/logs/download_bc09_gex.log" 2>/dev/null || true
