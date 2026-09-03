#!/usr/bin/env bash
# Generate one 500-cell receptor repertoire that every depth reuses.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
TRUTH="$ROOT/simulation/truth"
REF="$ROOT/simulation/reference"
LOG="$ROOT/logs/yasim_truth.log"
PARAMS="$ROOT/logs/yasim_parameters.tsv"
export PYTHONPATH="$ROOT/repo/yasim-sctcr${PYTHONPATH:+:$PYTHONPATH}"
export LOG_FILE_NAME="$ROOT/logs/yasim_sctcr_runtime.log"
mkdir -p "$TRUTH" "$(dirname "$LOG")"
for file in tcr_cache.json usage_bias.min.json cdr3_deletion_table.min.json cdr3_insertion_table.min.json; do
    [[ -s "$REF/$file" ]] || { echo "ERROR: missing $REF/$file; run download_yasim_reference.sh" >&2; exit 2; }
done
python -c 'import yasim_sctcr' >/dev/null 2>&1 || { echo "ERROR: yasim-sctcr is not installed in the active environment." >&2; exit 3; }
require_command seqkit

barcodes="$TRUTH/barcodes.txt"
python "$SCRIPT_DIR/generate_barcodes.py" --cells 500 --output "$barcodes" | tee -a "$LOG"
prefix="$TRUTH/sim_tcr"
if [[ ! -s "$prefix.stats" && ! -s "$prefix.stats.tsv" ]]; then
    python "$SCRIPT_DIR/yasim_sctcr_compat.py" rearrange_tcr \
        --tcr_cache_path "$REF/tcr_cache.json" \
        --cdr3_deletion_table_path "$REF/cdr3_deletion_table.min.json" \
        --cdr3_insertion_table_path "$REF/cdr3_insertion_table.min.json" \
        --usage_bias_json "$REF/usage_bias.min.json" \
        --portion_non_productive 0 -n 500 -o "$prefix" >> "$LOG" 2>&1
fi
stats="$prefix.stats"; [[ -s "$stats" ]] || stats="$prefix.stats.tsv"
cell_fasta="$TRUTH/sim_t_cell.nt.fa"
if [[ ! -s "$cell_fasta" ]]; then
    python "$SCRIPT_DIR/yasim_sctcr_compat.py" generate_tcr_clonal_expansion \
        -b "$barcodes" --src_tcr_stats_tsv "$stats" --dst_nt_fasta "$cell_fasta" --alpha 1 >> "$LOG" 2>&1
fi
barcode_map="$cell_fasta.stats"; [[ -s "$barcode_map" ]] || barcode_map="$cell_fasta.stats.tsv"
[[ -s "$barcode_map" ]] || { echo "ERROR: YASIM barcode map not found beside $cell_fasta" >&2; exit 4; }
python "$SCRIPT_DIR/build_simulation_truth.py" --tcr-stats "$stats" \
    --barcode-map "$barcode_map" --output "$TRUTH/simulation_ground_truth.tsv" | tee -a "$LOG"

rc_fasta="$TRUTH/sim_t_cell.rc.nt.fa"
if [[ ! -s "$rc_fasta" ]]; then
    seqkit seq --seq-type DNA --validate-seq --reverse --complement "$cell_fasta" > "$rc_fasta"
fi
if [[ ! -d "$rc_fasta.d" ]]; then
    python -m labw_utils.bioutils split_fasta "$rc_fasta" >> "$LOG" 2>&1
fi
{
    printf 'parameter\tvalue\n'
    printf 'cells\t500\nread_length\t150\nportion_non_productive\t0\n'
    printf 'repertoire_reused_across_depths\ttrue\n'
    printf 'random_seed\tNA: yasim-sctcr 1.0.1 uses random.SystemRandom and exposes no seed option\n'
    printf 'yasim_sctcr_version\t'; python -c 'import yasim_sctcr; print(yasim_sctcr.__version__)'
    printf 'art_version\t'; (art_illumina 2>&1 || true) | grep -m1 'Version'
} > "$PARAMS"
