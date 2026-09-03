#!/usr/bin/env bash
# Usage: run_mixcr_simulation.sh 2|10|50|100
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
DEPTH="${1:-}"
[[ "$DEPTH" =~ ^(2|10|50|100)$ ]] || { echo "Usage: $0 2|10|50|100" >&2; exit 2; }
require_command mixcr
if ! mixcr -v >/dev/null 2>&1 || ! mixcr listPresets 2>/dev/null | grep -Fq '10x-sc-5gex'; then
    echo "WAITING_FOR_MIXCR_LICENSE: usable MiXCR with 10x-sc-5gex is required." >&2
    exit 10
fi
input="$ROOT/simulation/depth_${DEPTH}x"
barcodes="$ROOT/simulation/truth/barcodes.txt"
r1="$input/BC09SIM_${DEPTH}x_1.fastq.gz"
r2="$input/BC09SIM_${DEPTH}x_2.fastq.gz"
[[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: simulated FASTQs missing for ${DEPTH}x." >&2; exit 3; }
[[ -s "$barcodes" ]] || { echo "ERROR: simulation barcode whitelist is missing: $barcodes" >&2; exit 3; }
outdir="$ROOT/results/mixcr/simulation/depth_${DEPTH}x"
mkdir -p "$outdir"
prefix="$outdir/MiXCR_SIM_${DEPTH}x"
# Synthetic cell barcodes are intentionally not members of the 10x 737K
# whitelist. Use the exact shared 500-cell truth list rather than allowing the
# native preset to discard every simulated read during tag refinement.
cmd=(mixcr analyze 10x-sc-5gex --species hsa \
    --set-whitelist "CELL=file:$barcodes" "$r1" "$r2" "$prefix")
printf '%q ' "${cmd[@]}" > "$ROOT/logs/mixcr_simulation_${DEPTH}x_command.txt"; printf '\n' >> "$ROOT/logs/mixcr_simulation_${DEPTH}x_command.txt"
"${cmd[@]}" > "$ROOT/logs/mixcr_simulation_${DEPTH}x.log" 2>&1 || {
    grep -Eqi 'licen[cs]e|activation' "$ROOT/logs/mixcr_simulation_${DEPTH}x.log" && exit 10
    exit 4
}
clns="$(find "$outdir" -maxdepth 1 -type f \( -name '*.clns' -o -name '*.clna' \) | sort | tail -n 1)"
[[ -s "$clns" ]] || { echo "ERROR: no MiXCR clone file at ${DEPTH}x." >&2; exit 5; }
exported="$outdir/MiXCR_SIM_${DEPTH}x_cell_receptors.tsv"
mixcr exportClones --drop-default-fields --split-by-tags Cell -cellId -vGene -jGene \
    -aaFeature CDR3 --chains TRA,TRB "$clns" "$exported" >> "$ROOT/logs/mixcr_simulation_${DEPTH}x.log" 2>&1
python "$SCRIPT_DIR/parse_mixcr.py" --input "$exported" \
    --output "$ROOT/data/processed/simulation_${DEPTH}x_mixcr_prediction.tsv"
