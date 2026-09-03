#!/usr/bin/env bash
# Usage: run_trust4_simulation.sh 2|10|50|100
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
DEPTH="${1:-}"
THREADS="${2:-8}"
[[ "$DEPTH" =~ ^(2|10|50|100)$ ]] || { echo "Usage: $0 2|10|50|100 [threads]" >&2; exit 2; }
require_command run-trust4

find_reference() {
    local candidate
    for candidate in \
        "$ROOT/repo/TRUST4/human_IMGT+C.fa" \
        "${TRUST4_REFERENCE_DIR:-}/human_IMGT+C.fa" \
        "${CONDA_PREFIX:-}/share/trust4/human_IMGT+C.fa" \
        "${CONDA_PREFIX:-}/share/trust4-"*/human_IMGT+C.fa \
        "$(dirname "$(command -v run-trust4)")/human_IMGT+C.fa" \
        "$(dirname "$(command -v run-trust4)")/../share/trust4/"*/human_IMGT+C.fa; do
        [[ -f "$candidate" ]] && { printf '%s\n' "$candidate"; return; }
    done
    return 1
}
imgt="$(find_reference || true)"
[[ -n "$imgt" ]] || { echo "ERROR: human_IMGT+C.fa not found; set TRUST4_REFERENCE_DIR." >&2; exit 3; }
input="$ROOT/simulation/depth_${DEPTH}x"
r1="$input/BC09SIM_${DEPTH}x_1.fastq.gz"
r2="$input/BC09SIM_${DEPTH}x_2.fastq.gz"
[[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: simulated FASTQs missing for ${DEPTH}x." >&2; exit 4; }
outdir="$ROOT/results/trust4/simulation/depth_${DEPTH}x"
mkdir -p "$outdir"
cmd=(run-trust4 -f "$imgt" --ref "$imgt" -u "$r2" --barcode "$r1" --UMI "$r1" \
    --readFormat 'bc:0:15,um:16:25' -t "$THREADS" -o "TRUST4_SIM_${DEPTH}x" --od "$outdir")
printf '%q ' "${cmd[@]}" > "$ROOT/logs/trust4_simulation_${DEPTH}x_command.txt"; printf '\n' >> "$ROOT/logs/trust4_simulation_${DEPTH}x_command.txt"
"${cmd[@]}" > "$ROOT/logs/trust4_simulation_${DEPTH}x.log" 2>&1
report=""
for candidate in "$outdir"/*_barcode_airr.tsv "$outdir"/*_barcode_report.tsv; do
    [[ -s "$candidate" ]] && { report="$candidate"; break; }
done
[[ -n "$report" ]] || { echo "ERROR: TRUST4 barcode result absent at ${DEPTH}x." >&2; exit 5; }
python "$SCRIPT_DIR/parse_trust4.py" --input "$report" \
    --output "$ROOT/data/processed/simulation_${DEPTH}x_trust4_prediction.tsv"
