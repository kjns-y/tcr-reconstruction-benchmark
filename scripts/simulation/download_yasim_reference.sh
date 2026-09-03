#!/usr/bin/env bash
# Download only the small YASIM-scTCR receptor reference tables (not scRNA data).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
OUT="$ROOT/simulation/reference"
LOG="$ROOT/logs/yasim_reference_download.log"
mkdir -p "$OUT" "$(dirname "$LOG")"

files=(
  tcr_cache.json.xz
  usage_bias.min.json.xz
  cdr3_deletion_table.min.json.xz
  cdr3_insertion_table.min.json.xz
)
for file in "${files[@]}"; do
    compressed="$OUT/$file"
    plain="$OUT/${file%.xz}"
    if [[ -s "$plain" ]]; then
        echo "SKIP existing: $plain" | tee -a "$LOG"
        continue
    fi
    if [[ ! -s "$compressed" ]]; then
        url="https://zenodo.org/api/records/12513698/files/$file/content"
        tmp="${compressed}.partial"
        curl --fail --location --retry 3 --output "$tmp" "$url" 2>> "$LOG"
        mv "$tmp" "$compressed"
    fi
    xz -t "$compressed"
    xz -dc "$compressed" > "${plain}.partial"
    mv "${plain}.partial" "$plain"
    sha256sum "$compressed" "$plain" | tee -a "$LOG"
done

