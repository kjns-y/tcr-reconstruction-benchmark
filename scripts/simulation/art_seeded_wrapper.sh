#!/usr/bin/env bash
# Give each receptor a stable ART seed derived from one recorded base seed.
set -euo pipefail

BASE_SEED="${TCRBENCH_ART_BASE_SEED:-20260829}"
SEED_LOG="${TCRBENCH_ART_SEED_LOG:-}"
REAL_ART="${TCRBENCH_ART_REAL_EXECUTABLE:-$(command -v art_illumina 2>/dev/null || true)}"
[[ -x "$REAL_ART" ]] || { echo "ERROR: art_illumina is not executable." >&2; exit 127; }
[[ "$BASE_SEED" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid TCRBENCH_ART_BASE_SEED." >&2; exit 2; }

source_fasta=""
arguments=("$@")
for ((index = 0; index < ${#arguments[@]}; index++)); do
    if [[ "${arguments[$index]}" == "--in" || "${arguments[$index]}" == "-i" ]]; then
        (( index + 1 < ${#arguments[@]} )) || { echo "ERROR: ART input argument lacks a value." >&2; exit 2; }
        source_fasta="${arguments[$((index + 1))]}"
        break
    fi
done
[[ -n "$source_fasta" ]] || { echo "ERROR: could not find ART --in argument." >&2; exit 2; }

receptor_id="$(basename "$source_fasta")"
digest="$(printf '%s' "${BASE_SEED}:${receptor_id}" | sha256sum | awk '{print $1}')"
seed=$((16#${digest:0:8} % 2147483646 + 1))
if [[ -n "$SEED_LOG" ]]; then
    printf '%s\t%s\n' "$receptor_id" "$seed" >> "$SEED_LOG"
fi

exec "$REAL_ART" "${arguments[@]}" --rndSeed "$seed"
