#!/usr/bin/env bash
# Install a pinned MiXCR distribution into the isolated tcrbench environment.
# This script does not activate or store a license. The user must obtain and
# activate a valid license through the official MiLaboratories process.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"
ROOT="$(project_root)"
VERSION="${1:-4.7.0}"
LOG="$ROOT/logs/install_mixcr.log"

[[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || {
    echo "ERROR: MiXCR version must have the form X.Y.Z." >&2
    exit 2
}
require_command java
java_major="$(java -version 2>&1 | awk -F '[\".]' 'NR == 1 {print $2}')"
[[ "$java_major" =~ ^[0-9]+$ ]] || { echo "ERROR: could not determine Java version." >&2; exit 3; }
(( java_major >= 11 )) || { echo "ERROR: MiXCR requires Java 11 or newer." >&2; exit 3; }

if command -v mixcr >/dev/null 2>&1; then
    installed="$(mixcr -v 2>&1 | sed -n '1p' || true)"
    if grep -Fq "v$VERSION" <<<"$installed"; then
        echo "MiXCR $VERSION is already installed: $(command -v mixcr)" | tee -a "$LOG"
        exit 0
    fi
    echo "ERROR: another MiXCR is already on PATH: $installed" >&2
    echo "Refusing to replace it implicitly." >&2
    exit 4
fi

micromamba_bin="${MAMBA_EXE:-}"
for candidate in "$micromamba_bin" "$(command -v micromamba 2>/dev/null || true)" \
    "$HOME/.local/bin/micromamba" "$ROOT/../tcr_projects/environments/bin/micromamba"; do
    [[ -x "$candidate" ]] && { micromamba_bin="$candidate"; break; }
done
[[ -x "$micromamba_bin" ]] || {
    echo "ERROR: micromamba is required for the separate official-channel installation." >&2
    exit 5
}

{
    echo "[$(timestamp)] Installing MiXCR $VERSION from the official milaboratories channel."
    echo "Target environment: tcrbench"
    echo "Java: $(java -version 2>&1 | sed -n '1p')"
} | tee -a "$LOG"

# System Java already satisfies the declared dependency. Avoid modifying the
# rest of the scientific environment merely to install the distribution JAR.
"$micromamba_bin" install -y -n tcrbench --override-channels \
    -c milaboratories -c conda-forge --no-deps "mixcr=$VERSION" >>"$LOG" 2>&1

hash -r
require_command mixcr
mixcr -v >"$ROOT/logs/mixcr_version.txt" 2>&1 || {
    echo "MiXCR was installed, but execution/license validation failed." >&2
    echo "Activate a valid license with the official 'mixcr activate-license' workflow." >&2
    exit 10
}
echo "Installed MiXCR $VERSION: $(command -v mixcr)" | tee -a "$LOG"
