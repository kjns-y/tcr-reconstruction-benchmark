#!/usr/bin/env bash
# Generate four depth tables from the same 500-cell receptor truth.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
export PYTHONPATH="$ROOT/repo/yasim-sctcr${PYTHONPATH:+:$PYTHONPATH}"
export LOG_FILE_NAME="$ROOT/logs/yasim_sctcr_runtime.log"
BARCODES="$ROOT/simulation/truth/barcodes.txt"
[[ -s "$BARCODES" ]] || { echo "ERROR: run 01_generate_tcr_truth.sh first." >&2; exit 2; }
for depth in 2 10 50 100; do
    out="$ROOT/simulation/truth/scTCR.depth_${depth}x.tsv"
    python "$SCRIPT_DIR/yasim_sctcr_compat.py" generate_tcr_depth -b "$BARCODES" -o "$out" -d "$depth"
done
