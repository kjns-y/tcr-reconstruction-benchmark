#!/usr/bin/env bash

# Shared path helpers. Source this file; do not execute it directly.
set -o pipefail

project_root() {
    local source_path
    source_path="${BASH_SOURCE[0]}"
    while [[ -L "$source_path" ]]; do
        source_path="$(readlink "$source_path")"
    done
    cd "$(dirname "$source_path")/../.." >/dev/null 2>&1 && pwd
}

require_command() {
    local command_name="$1"
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "ERROR: required command is missing: $command_name" >&2
        return 127
    fi
}

timestamp() {
    date '+%Y-%m-%dT%H:%M:%S%z'
}

