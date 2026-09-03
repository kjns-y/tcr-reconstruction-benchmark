#!/usr/bin/env bash
# Source this file: source scripts/setup/activate_tcrbench.sh
if [[ -n "${BASH_VERSION:-}" && "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "ERROR: source this script instead of executing it." >&2
    exit 2
fi
if [[ -f "$PWD/scripts/setup/activate_tcrbench.sh" ]]; then
    _tcrbench_root="$PWD"
elif [[ -n "${BASH_VERSION:-}" ]]; then
    _tcrbench_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
else
    echo "ERROR: change to the project root before sourcing this script." >&2
    return 2
fi
_tcrbench_mamba="${MAMBA_EXE:-}"
for _candidate in "$_tcrbench_mamba" "$(command -v micromamba 2>/dev/null || true)" \
    "$HOME/.local/bin/micromamba" "$_tcrbench_root/../tcr_projects/environments/bin/micromamba"; do
    [[ -x "$_candidate" ]] && { _tcrbench_mamba="$_candidate"; break; }
done
if [[ -x "$_tcrbench_mamba" ]]; then
    export MAMBA_EXE="$_tcrbench_mamba"
    export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/share/mamba}"
    eval "$("$_tcrbench_mamba" shell hook --shell bash)"
    micromamba activate tcrbench
    unset _candidate _tcrbench_mamba _tcrbench_root
    return 0
fi
_tcrbench_conda="${CONDA_EXE:-}"
for _candidate in "$_tcrbench_conda" "$(command -v conda 2>/dev/null || true)" "$HOME/anaconda3/bin/conda" \
    "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    [[ -x "$_candidate" ]] && { _tcrbench_conda="$_candidate"; break; }
done
[[ -x "$_tcrbench_conda" ]] || { echo "ERROR: conda/micromamba not found." >&2; return 2; }
source "$(dirname "$(dirname "$_tcrbench_conda")")/etc/profile.d/conda.sh"
conda activate tcrbench
unset _candidate _tcrbench_conda _tcrbench_root
