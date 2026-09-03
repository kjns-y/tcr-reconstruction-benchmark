#!/usr/bin/env bash
# Create an isolated tcrbench environment; never modifies base in place.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
YML="$ROOT/envs/environment.yml"

micromamba_bin="${MAMBA_EXE:-}"
for candidate in "$micromamba_bin" "$(command -v micromamba 2>/dev/null || true)" \
    "$HOME/.local/bin/micromamba" "$ROOT/../tcr_projects/environments/bin/micromamba"; do
    [[ -x "$candidate" ]] && { micromamba_bin="$candidate"; break; }
done
if [[ -x "$micromamba_bin" ]]; then
    if "$micromamba_bin" env list | awk 'NR>2 {print $1}' | grep -Fxq tcrbench; then
        echo "Environment tcrbench already exists; no changes made."
        exit 0
    fi
    "$micromamba_bin" create -y -f "$YML"
    echo "Created isolated Micromamba environment: tcrbench"
    echo "Activate with: source scripts/setup/activate_tcrbench.sh"
    exit 0
fi

conda_bin="${CONDA_EXE:-}"
for candidate in "$conda_bin" "$(command -v conda 2>/dev/null || true)" "$HOME/anaconda3/bin/conda" \
    "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    [[ -x "$candidate" ]] && { conda_bin="$candidate"; break; }
done
[[ -x "$conda_bin" ]] || { echo "ERROR: conda/micromamba was not found. Install Miniforge without changing base." >&2; exit 2; }

if "$conda_bin" env list | awk '{print $1}' | grep -Fxq tcrbench; then
    echo "Environment tcrbench already exists; no changes made."
    echo "To update explicitly: $conda_bin env update -n tcrbench -f $YML --prune"
    exit 0
fi
"$conda_bin" env create -f "$YML"
echo "Created isolated environment: tcrbench"
echo "Activate with: source \"$(dirname "$(dirname "$conda_bin")")/etc/profile.d/conda.sh\" && conda activate tcrbench"

