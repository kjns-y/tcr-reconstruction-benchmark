#!/usr/bin/env bash
# Stage dispatcher. With no stage it prints help and performs no downloads.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
stage="${1:-help}"
threads="${2:-8}"

case "$stage" in
    setup)
        "$ROOT/scripts/setup/check_environment.sh"
        ;;
    ground_truth)
        "$ROOT/scripts/download/resolve_accessions.sh"
        "$ROOT/scripts/download/download_ground_truth.sh"
        python "$ROOT/scripts/preprocessing/build_ground_truth.py"
        ;;
    preview)
        "$ROOT/scripts/download/preview_bc09_gex.sh" "${2:-10000}"
        ;;
    download)
        echo "This stage starts the large BC09 GEX SRA download after a space check."
        "$ROOT/scripts/download/download_bc09_gex.sh" "$threads"
        "$ROOT/scripts/preprocessing/check_fastq.sh"
        ;;
    trust4)
        "$ROOT/scripts/trust4/run_trust4_bc09.sh" "$threads"
        ;;
    mixcr)
        "$ROOT/scripts/mixcr/run_mixcr_bc09.sh" "$threads"
        ;;
    evaluate)
        truth="$ROOT/data/ground_truth/BC09_ground_truth.tsv"
        if [[ -s "$ROOT/data/processed/BC09_trust4_prediction.tsv" ]]; then
            python "$ROOT/scripts/evaluation/evaluate.py" --truth "$truth" \
                --prediction "$ROOT/data/processed/BC09_trust4_prediction.tsv" --method TRUST4 \
                --output "$ROOT/results/benchmark/BC09_TRUST4_metrics.tsv"
        fi
        if [[ -s "$ROOT/data/processed/BC09_mixcr_prediction.tsv" ]]; then
            python "$ROOT/scripts/evaluation/evaluate.py" --truth "$truth" \
                --prediction "$ROOT/data/processed/BC09_mixcr_prediction.tsv" --method MiXCR \
                --output "$ROOT/results/benchmark/BC09_MiXCR_metrics.tsv"
        fi
        if [[ -s "$ROOT/results/benchmark/BC09_TRUST4_cell_details.tsv" && -s "$ROOT/results/benchmark/BC09_MiXCR_cell_details.tsv" ]]; then
            python "$ROOT/scripts/evaluation/select_example_cells.py"
        fi
        ;;
    simulate)
        "$ROOT/scripts/simulation/download_yasim_reference.sh"
        "$ROOT/scripts/simulation/inspect_yasim_cli.sh"
        "$ROOT/scripts/simulation/01_generate_tcr_truth.sh"
        "$ROOT/scripts/simulation/02_generate_depths.sh"
        "$ROOT/scripts/simulation/03_simulate_reads.sh" "$threads"
        for depth in 2 10 50 100; do
            [[ -s "$ROOT/data/processed/simulation_${depth}x_trust4_prediction.tsv" ]] || \
                "$ROOT/scripts/trust4/run_trust4_simulation.sh" "$depth" "$threads"
        done
        if command -v mixcr >/dev/null 2>&1 && mixcr -v >/dev/null 2>&1 && \
            mixcr listPresets 2>/dev/null | grep -Fq '10x-sc-5gex'; then
            for depth in 2 10 50 100; do
                [[ -s "$ROOT/data/processed/simulation_${depth}x_mixcr_prediction.tsv" ]] || \
                    "$ROOT/scripts/mixcr/run_mixcr_simulation.sh" "$depth"
            done
        else
            echo "WAITING_FOR_MIXCR_LICENSE: evaluating available TRUST4 simulation results only." >&2
        fi
        python "$ROOT/scripts/evaluation/evaluate_simulation.py"
        ;;
    plot)
        if [[ -s "$ROOT/results/benchmark/BC09_TRUST4_metrics.tsv" || \
              -s "$ROOT/results/benchmark/BC09_MiXCR_metrics.tsv" ]]; then
            python "$ROOT/scripts/plotting/plot_realdata_benchmark.py"
        else
            echo "SKIP real-data figure: no complete BC09 metrics are available yet."
        fi
        [[ -s "$ROOT/results/benchmark/simulation_metrics.tsv" ]] && python "$ROOT/scripts/plotting/plot_depth_benchmark.py"
        ;;
    help|-h|--help)
        cat <<'EOF'
Usage: ./run_pipeline.sh STAGE [threads]

Stages:
  setup          inspect tools only; does not install
  ground_truth   resolve NCBI accessions and build the small official truth
  preview        download at most 100,000 spots to verify real FASTQ structure
  download       LARGE: download BC09 GEX SRA data after a disk-space gate
  trust4         run TRUST4 on real BC09 GEX
  mixcr          run MiXCR if preset/license are usable
  evaluate       compute available real-data metrics and example cells
  simulate       generate, reconstruct, and evaluate the four-depth simulation
  plot           plot whichever metrics are available

With no stage, nothing is downloaded or executed.
EOF
        ;;
    *)
        echo "ERROR: unknown stage: $stage" >&2
        exit 2
        ;;
esac
