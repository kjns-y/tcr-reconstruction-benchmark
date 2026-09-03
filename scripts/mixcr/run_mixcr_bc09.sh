#!/usr/bin/env bash
# Run the installed MiXCR version only when its preset and license are usable.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
THREADS="${1:-8}"
FASTQ_DIR="${BC09_FASTQ_DIR:-$ROOT/data/raw/BC09_GEX}"
OUTDIR="$ROOT/results/mixcr/BC09"
INPUT_DIR="$OUTDIR/input_links"
COMMAND_LOG="$ROOT/logs/mixcr_bc09_command.txt"
RUN_LOG="$ROOT/logs/mixcr_bc09.log"
require_command mixcr
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be a positive integer." >&2; exit 2; }

if ! mixcr -v > "$ROOT/logs/mixcr_version.txt" 2>&1; then
    echo "WAITING_FOR_MIXCR_LICENSE: mixcr -v failed; see logs/mixcr_version.txt" >&2
    exit 10
fi
if ! mixcr listPresets > "$ROOT/logs/mixcr_presets.txt" 2>&1; then
    echo "WAITING_FOR_MIXCR_LICENSE: cannot list presets; see logs/mixcr_presets.txt" >&2
    exit 10
fi
grep -Fq '10x-sc-5gex' "$ROOT/logs/mixcr_presets.txt" || {
    echo "ERROR: installed MiXCR does not provide preset 10x-sc-5gex." >&2
    exit 11
}
mapfile -t r1 < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_1.fastq.gz' | sort)
mapfile -t r2 < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_2.fastq.gz' | sort)
(( ${#r1[@]} > 0 && ${#r1[@]} == ${#r2[@]} )) || { echo "ERROR: paired FASTQ sets are absent or unequal." >&2; exit 12; }
mkdir -p "$INPUT_DIR"
for i in "${!r1[@]}"; do
    lane="$(printf '%03d' "$((i + 1))")"
    ln -sfn "${r1[$i]}" "$INPUT_DIR/BC09_L${lane}_R1_001.fastq.gz"
    ln -sfn "${r2[$i]}" "$INPUT_DIR/BC09_L${lane}_R2_001.fastq.gz"
done

prefix="$OUTDIR/BC09_mixcr"
cmd=(mixcr analyze 10x-sc-5gex --species hsa --threads "$THREADS" \
    "$INPUT_DIR/BC09_L{{n}}_R1_001.fastq.gz" \
    "$INPUT_DIR/BC09_L{{n}}_R2_001.fastq.gz" "$prefix")
printf '%q ' "${cmd[@]}" > "$COMMAND_LOG"; printf '\n' >> "$COMMAND_LOG"
"${cmd[@]}" > "$RUN_LOG" 2>&1 || {
    if grep -Eqi 'licen[cs]e|activation' "$RUN_LOG"; then
        echo "WAITING_FOR_MIXCR_LICENSE: analysis blocked by MiXCR licensing; see $RUN_LOG" >&2
        exit 10
    fi
    echo "ERROR: MiXCR analysis failed; see $RUN_LOG" >&2
    exit 13
}

clns="$(find "$OUTDIR" -maxdepth 1 -type f \( -name '*.clns' -o -name '*.clna' \) | sort | tail -n 1)"
[[ -s "$clns" ]] || { echo "ERROR: no MiXCR .clns/.clna result found." >&2; exit 14; }
exported="$OUTDIR/BC09_mixcr_cell_receptors.tsv"
export_cmd=(mixcr exportClones --drop-default-fields --split-by-tags Cell \
    -cellId -vGene -jGene -aaFeature CDR3 --chains TRA,TRB "$clns" "$exported")
printf '%q ' "${export_cmd[@]}" >> "$COMMAND_LOG"; printf '\n' >> "$COMMAND_LOG"
"${export_cmd[@]}" >> "$RUN_LOG" 2>&1
python3 "$SCRIPT_DIR/parse_mixcr.py" --input "$exported" --output "$ROOT/data/processed/BC09_mixcr_prediction.tsv"
