#!/usr/bin/env bash
# Source with a site config; BDX_RUNTIME_ROOT identifies this private GEMC.
set -eo pipefail
set +u
source "$1"
: "${BDX_RUNTIME_ROOT:?Set BDX_RUNTIME_ROOT before sourcing environment.sh}"
export GEMC="$BDX_RUNTIME_ROOT" GEMC_DATA_DIR="$BDX_RUNTIME_ROOT" BDX_NO_EVIO=1
command -v geant4-config >/dev/null
if [[ -n ${BDX_EXPECTED_GEANT4:-} && $(geant4-config --version) != "$BDX_EXPECTED_GEANT4" ]]; then
    echo "Geant4 version differs from BDX_EXPECTED_GEANT4=$BDX_EXPECTED_GEANT4" >&2
    return 1
fi
set -u
