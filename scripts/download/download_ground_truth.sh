#!/usr/bin/env bash
# Download the small official processed VDJ annotation without overwriting it.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3148nnn/GSM3148580/suppl/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv.gz"
OUT="$ROOT/data/ground_truth/raw/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv.gz"
LOG="$ROOT/logs/download_ground_truth.log"
mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"

if [[ -s "$OUT" ]]; then
    gzip -t "$OUT"
    echo "Already present and valid: $OUT" | tee -a "$LOG"
    exit 0
fi

tmp="${OUT}.partial"
trap 'rm -f "$tmp"' EXIT
{
    echo "[$(timestamp)] GET $URL"
    curl --fail --location --retry 3 --retry-delay 3 --output "$tmp" "$URL"
    gzip -t "$tmp"
    mv "$tmp" "$OUT"
    sha256sum "$OUT"
} 2>&1 | tee -a "$LOG"
trap - EXIT
