#!/usr/bin/env bash
# Download a bounded BC09 GEX preview to verify the real read structure.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
SPOTS="${1:-10000}"
META="$ROOT/data/metadata/GSM3148575_accessions.tsv"
OUTDIR="$ROOT/data/raw/BC09_GEX_preview"
STATS="$ROOT/data/metadata/BC09_fastq_preview_stats.tsv"
INFERENCE="$ROOT/data/metadata/BC09_fastq_preview_inference.tsv"
KNOWN_BARCODES="$ROOT/data/ground_truth/BC09_ground_truth.tsv"
LOG="$ROOT/logs/preview_bc09_gex.log"

for command_name in fastq-dump seqkit python3 awk gzip; do
    require_command "$command_name"
done
[[ -s "$META" ]] || {
    echo "ERROR: missing metadata: $META. Run resolve_accessions.sh first." >&2
    exit 2
}
[[ "$SPOTS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: spots must be a positive integer." >&2
    exit 2
}
if (( SPOTS > 100000 )); then
    echo "ERROR: preview is capped at 100000 spots; use the reviewed full downloader for larger data." >&2
    exit 2
fi

run="$(awk -F '\t' 'NR == 2 {print $4}' "$META")"
[[ "$run" =~ ^SRR[0-9]+$ ]] || {
    echo "ERROR: metadata did not yield a valid first SRR accession: $run" >&2
    exit 3
}

mkdir -p "$OUTDIR" "$(dirname "$LOG")"
r1="$OUTDIR/${run}_1.fastq.gz"
r2="$OUTDIR/${run}_2.fastq.gz"

if [[ -s "$r1" && -s "$r2" ]]; then
    gzip -t "$r1" "$r2"
    echo "[$(timestamp)] SKIP existing validated preview pair: $r1, $r2" | tee -a "$LOG"
elif [[ -e "$r1" || -e "$r2" ]]; then
    echo "ERROR: partial preview output exists; inspect it before retrying: $OUTDIR" | tee -a "$LOG" >&2
    exit 4
else
    echo "[$(timestamp)] Downloading first $SPOTS spots from $run" | tee -a "$LOG"
    fastq-dump \
        --maxSpotId "$SPOTS" \
        --split-files \
        --gzip \
        --outdir "$OUTDIR" \
        "$run" >>"$LOG" 2>&1
    [[ -s "$r1" && -s "$r2" ]] || {
        echo "ERROR: preview did not produce the expected paired FASTQs." | tee -a "$LOG" >&2
        exit 5
    }
    gzip -t "$r1" "$r2"
fi

seqkit stats --tabular --all "$r1" "$r2" >"$STATS" 2>>"$LOG"
inference_args=(
    --read1 "$r1"
    --read2 "$r2"
    --output "$INFERENCE"
    --sample-reads "$SPOTS"
)
if [[ -s "$KNOWN_BARCODES" ]]; then
    inference_args+=(--known-barcodes "$KNOWN_BARCODES")
fi
python3 "$ROOT/scripts/preprocessing/infer_fastq_structure.py" \
    "${inference_args[@]}" 2>>"$LOG" | tee -a "$LOG"

echo "Saved preview FASTQs: $OUTDIR" | tee -a "$LOG"
echo "Saved stats: $STATS" | tee -a "$LOG"
echo "Saved structure inference: $INFERENCE" | tee -a "$LOG"
