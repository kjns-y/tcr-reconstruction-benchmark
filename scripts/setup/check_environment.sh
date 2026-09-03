#!/usr/bin/env bash
# Check every required executable without stopping when one is missing.
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
OUT="$ROOT/logs/environment_check.txt"
mkdir -p "$(dirname "$OUT")"

version_for() {
    local tool="$1"
    local executable="$2"
    case "$tool" in
        python|python3|git|conda|mamba) "$executable" --version 2>&1 | head -n 1 ;;
        java) "$executable" -version 2>&1 | head -n 1 ;;
        wget|curl|pigz) "$executable" --version 2>&1 | head -n 1 ;;
        seqkit) "$executable" version 2>&1 | head -n 1 ;;
        samtools) "$executable" --version 2>&1 | head -n 1 ;;
        prefetch|fasterq-dump|fastq-dump) "$executable" -V 2>&1 | sed -n '/[^[:space:]]/ {p;q;}' ;;
        run-trust4) "$executable" --version 2>&1 | head -n 1 || "$executable" 2>&1 | head -n 1 ;;
        mixcr) "$executable" -v 2>&1 | head -n 1 ;;
        art_illumina) "$executable" 2>&1 | grep -m1 'Version' || true ;;
        *) "$executable" --version 2>&1 | head -n 1 ;;
    esac
}

path_for() {
    local tool="$1" candidate
    command -v "$tool" 2>/dev/null && return 0
    case "$tool" in
        conda)
            for candidate in "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
                [[ -x "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
            done
            ;;
        mamba)
            for candidate in "$HOME/.local/bin/micromamba" "$ROOT/../tcr_projects/environments/bin/micromamba"; do
                [[ -x "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
            done
            ;;
    esac
    return 1
}

{
    printf 'tool\tstatus\tversion\tpath\n'
    for tool in git python python3 conda mamba java wget curl pigz seqkit samtools prefetch fasterq-dump fastq-dump run-trust4 mixcr art_illumina; do
        if path="$(path_for "$tool")"; then
            version="$(version_for "$tool" "$path" | tr '\t\r\n' '   ' | sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//')"
            [[ -n "$version" ]] || version="version unavailable"
            printf '%s\tOK\t%s\t%s\n' "$tool" "$version" "$path"
        else
            printf '%s\tMISSING\t-\t-\n' "$tool"
        fi
    done
    for module in yasim_sctcr pandas numpy scipy matplotlib yaml tqdm; do
        if python3 -c "import $module" >/dev/null 2>&1; then
            version="$(python3 -c "import $module as m; print(getattr(m, '__version__', 'installed'))" 2>/dev/null | head -n 1)"
            printf 'python:%s\tOK\t%s\t%s\n' "$module" "$version" "$(command -v python3)"
        else
            printf 'python:%s\tMISSING\t-\t-\n' "$module"
        fi
    done
} | tee "$OUT"

echo "Saved: $OUT" >&2
