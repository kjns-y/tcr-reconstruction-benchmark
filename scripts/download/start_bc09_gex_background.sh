#!/usr/bin/env bash
# Start the reviewed full BC09 download as a server-side resumable background job.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
THREADS="${1:-8}"
PID_FILE="$ROOT/logs/download_bc09_background.pid"
OUT="$ROOT/logs/download_bc09_background.out"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be a positive integer." >&2; exit 2; }
for command_name in prefetch fasterq-dump; do require_command "$command_name"; done

if [[ -s "$PID_FILE" ]]; then
    previous_pid="$(tr -d '[:space:]' < "$PID_FILE")"
    if [[ "$previous_pid" =~ ^[1-9][0-9]*$ ]] && kill -0 "$previous_pid" 2>/dev/null; then
        echo "Download is already active as PID $previous_pid."
        exit 0
    fi
fi

active="$(pgrep -f "$ROOT/scripts/download/download_bc09_gex.sh" || true)"
[[ -z "$active" ]] || {
    echo "ERROR: a BC09 downloader exists but is not represented by the PID file: $active" >&2
    exit 3
}

stale_locks="$(find "$ROOT/data/raw/sra" -type f -name '*.sra.lock' -print 2>/dev/null || true)"
[[ -z "$stale_locks" ]] || {
    echo "ERROR: SRA lock file(s) exist and no downloader is active:" >&2
    printf '%s\n' "$stale_locks" >&2
    echo "Verify that prefetch is absent, then archive the stale lock before resuming." >&2
    exit 4
}

mkdir -p "$ROOT/logs"
nohup "$ROOT/run_pipeline.sh" download "$THREADS" >"$OUT" 2>&1 < /dev/null &
job_pid=$!
printf '%s\n' "$job_pid" > "$PID_FILE"
sleep 2
kill -0 "$job_pid" 2>/dev/null || {
    echo "ERROR: background downloader exited immediately; see $OUT" >&2
    exit 5
}
echo "Started BC09 background download: PID $job_pid"
echo "Status: $ROOT/scripts/download/status_bc09_gex.sh"
