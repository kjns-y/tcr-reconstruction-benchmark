#!/usr/bin/env bash
# Simulate fixed-length 5'-oriented reads at each depth, then add 10x R1 tags.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
export PYTHONPATH="$ROOT/repo/yasim-sctcr${PYTHONPATH:+:$PYTHONPATH}"
export LOG_FILE_NAME="$ROOT/logs/yasim_runtime.log"
THREADS="${1:-8}"
ART_BASE_SEED="${TCRBENCH_ART_SEED:-20260829}"
PARAMS="$ROOT/logs/yasim_art_parameters.tsv"
SEED_LOG_RAW="$ROOT/logs/yasim_art_receptor_seeds.raw.tsv"
SEED_LOG="$ROOT/logs/yasim_art_receptor_seeds.tsv"
ART_WRAPPER="$SCRIPT_DIR/art_seeded_wrapper.sh"
FASTA_DIR="$ROOT/simulation/truth/sim_t_cell.rc.nt.fa.d"
[[ -d "$FASTA_DIR" ]] || { echo "ERROR: run 01_generate_tcr_truth.sh first." >&2; exit 2; }
[[ "$ART_BASE_SEED" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: TCRBENCH_ART_SEED must be a positive integer." >&2; exit 2; }
require_command art_illumina
require_command sha256sum
require_command seqkit
[[ -x "$ART_WRAPPER" ]] || { echo "ERROR: seeded ART wrapper is not executable: $ART_WRAPPER" >&2; exit 2; }
: > "$SEED_LOG_RAW"
export TCRBENCH_ART_BASE_SEED="$ART_BASE_SEED"
export TCRBENCH_ART_SEED_LOG="$SEED_LOG_RAW"
export TCRBENCH_ART_REAL_EXECUTABLE="$(command -v art_illumina)"
{
    printf 'parameter\tvalue\n'
    printf 'art_base_seed\t%s\n' "$ART_BASE_SEED"
    printf 'per_receptor_seed\tSHA256(base_seed:receptor_basename) mod 2147483646 + 1\n'
    printf 'sequencer_profile\tHS25\n'
    printf 'read_length\t150\n'
    printf 'layout\tsingle_end_5prime_amplicon\n'
    printf 'threads\t%s\n' "$THREADS"
} > "$PARAMS"
for depth in 2 10 50 100; do
    outdir="$ROOT/simulation/depth_${depth}x"
    prefix="$outdir/yasim_tcr_${depth}x"
    depth_file="$ROOT/simulation/truth/scTCR.depth_${depth}x.tsv"
    log="$ROOT/logs/yasim_art_depth_${depth}x.log"
    command_log="$ROOT/logs/yasim_art_depth_${depth}x_command.txt"
    mkdir -p "$outdir"
    [[ -s "$depth_file" ]] || { echo "ERROR: missing $depth_file; run 02_generate_depths.sh" >&2; exit 3; }
    cmd=(python -m yasim art -F "$FASTA_DIR" -o "$prefix" --sequencer_name HS25 \
        --read_length 150 -d "$depth_file" -e "$ART_WRAPPER" --preserve_intermediate_files \
        -j "$THREADS" --amplicon)
    printf '%q ' "${cmd[@]}" > "$command_log"; printf '\n' >> "$command_log"
    if [[ ! -s "$prefix.fq" ]]; then
        "${cmd[@]}" > "$log" 2>&1
    fi
    [[ -s "$prefix.fq" ]] || { echo "ERROR: ART/YASIM did not create $prefix.fq; see $log" >&2; exit 4; }
    r1="$outdir/BC09SIM_${depth}x_1.fastq.gz"
    r2="$outdir/BC09SIM_${depth}x_2.fastq.gz"
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        python "$SCRIPT_DIR/make_10x_fastqs.py" --input "$prefix.fq" --r1 "$r1" --r2 "$r2" | tee -a "$log"
    fi
done
{
    printf 'receptor_fasta\tart_seed\n'
    sort -u "$SEED_LOG_RAW"
} > "$SEED_LOG"
rm -f "$SEED_LOG_RAW"

mapfile -t simulation_fastqs < <(
    find "$ROOT/simulation" -maxdepth 2 -type f -name 'BC09SIM_*x_[12].fastq.gz' | sort
)
(( ${#simulation_fastqs[@]} == 8 )) || {
    echo "ERROR: expected eight simulated FASTQs, observed ${#simulation_fastqs[@]}." >&2
    exit 5
}
SIMULATION_STATS="$ROOT/data/metadata/simulation_fastq_stats.tsv"
seqkit stats --tabular --all "${simulation_fastqs[@]}" > "$SIMULATION_STATS"
python "$SCRIPT_DIR/validate_simulation_fastq.py" --stats "$SIMULATION_STATS"
