#!/usr/bin/env bash
# Capture the installed CLI rather than trusting documentation alone.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
OUT="$ROOT/logs/yasim_cli_help.txt"
subcommands=(rearrange_tcr generate_tcr_clonal_expansion generate_tcr_depth)
export PYTHONPATH="$ROOT/repo/yasim-sctcr${PYTHONPATH:+:$PYTHONPATH}"
export LOG_FILE_NAME="$ROOT/logs/yasim_sctcr_runtime.log"
{
    echo '### Raw installed CLI: python -m yasim_sctcr --help'
    python -m yasim_sctcr --help; echo "exit_code=$?"
    echo
    echo '### Raw installed CLI: python -m yasim_sctcr lscmd'
    python -m yasim_sctcr lscmd; echo "exit_code=$?"
    for subcommand in "${subcommands[@]}"; do
        echo
        echo "### Raw installed CLI: python -m yasim_sctcr $subcommand --help"
        python -m yasim_sctcr "$subcommand" --help; echo "exit_code=$?"
        echo
        echo "### Compatibility launcher: $subcommand --help"
        python "$SCRIPT_DIR/yasim_sctcr_compat.py" "$subcommand" --help
        echo "compat_${subcommand}_exit_code=$?"
    done
    echo
    echo '### python -m yasim art --help (bulk TCR depth table adapter)'
    python -m yasim art --help; echo "exit_code=$?"
    echo
    echo '### python -m yasim_sc art --help (expects a directory of per-cell depth files)'
    python -m yasim_sc art --help; echo "exit_code=$?"
} > "$OUT" 2>&1
if ! grep -q '^compat_rearrange_tcr_exit_code=0$' "$OUT"; then
    echo "ERROR: compatibility launcher did not provide usable rearrange_tcr help; see $OUT" >&2
    exit 1
fi
echo "Saved: $OUT"
