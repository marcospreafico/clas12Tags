#!/usr/bin/env bash
# Source safely from interactive Bash or strict batch scripts. Never enable
# errexit/nounset/pipefail in the calling shell. Modules may need nounset off.
_bdx_load_environment() {
    local bdx_had_nounset=0 bdx_status=0
    case $- in *u*) bdx_had_nounset=1 ;; esac
    set +u
    if [[ -z ${1:-} || ! -f $1 ]]; then
        echo 'Provide an existing BDX environment config' >&2
        bdx_status=1
    else
        source "$1"
        bdx_status=$?
    fi
    if [[ $bdx_status == 0 && -z ${BDX_RUNTIME_ROOT:-} ]]; then
        echo 'Set BDX_RUNTIME_ROOT before sourcing environment.sh' >&2
        bdx_status=1
    fi
    if [[ $bdx_status == 0 ]]; then
        if BDX_RUNTIME_ROOT=$(cd "$BDX_RUNTIME_ROOT" && pwd -P); then
            export BDX_RUNTIME_ROOT
            export GEMC="$BDX_RUNTIME_ROOT" GEMC_DATA_DIR="$BDX_RUNTIME_ROOT" BDX_NO_EVIO=1
        else
            bdx_status=1
        fi
    fi
    if [[ $bdx_status == 0 ]]; then
        if ! command -v geant4-config >/dev/null; then
            echo 'geant4-config missing after environment setup' >&2
            bdx_status=1
        elif [[ -n ${BDX_EXPECTED_GEANT4:-} && $(geant4-config --version) != "$BDX_EXPECTED_GEANT4" ]]; then
            echo "Geant4 version differs from BDX_EXPECTED_GEANT4=$BDX_EXPECTED_GEANT4" >&2
            bdx_status=1
        fi
    fi
    if [[ $bdx_had_nounset == 1 ]]; then set -u; fi
    return "$bdx_status"
}
_bdx_load_environment "$@"
