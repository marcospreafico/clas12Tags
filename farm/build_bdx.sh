#!/usr/bin/env bash
# Build on ifarm, inside a fresh copy so Mac/Docker objects cannot be reused.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ ${1:-} == --help ]]; then
    echo 'Usage: bash farm/build_bdx.sh [farm/ifarm.env] [compile cores=2]'; exit 0
fi
env_file=${1:-"$root/farm/ifarm.env"}
cores=${2:-2}
[[ $cores =~ ^[1-9][0-9]*$ && -f $env_file ]] || { echo 'Provide an environment file and positive core count' >&2; exit 2; }
release=$(python3 "$root/farm/release.py" prepare "$root" "$env_file")
export BDX_RUNTIME_ROOT="$release/runtime"
source "$root/farm/environment.sh" "$release/runtime/farm/ifarm.env"
(
    cd "$BDX_RUNTIME_ROOT/source"
    scons -j"$cores" OPT=1
    mkdir -p ../build/bdx_optics
    g++ -std=c++17 -O2 -I"$HIPO/include" ../scripts/bdx_hipo.cc \
        -L"$HIPO/lib" "-Wl,-rpath,$HIPO/lib" -lhipo4 -llz4 -lfmt -pthread \
        -o ../build/bdx_optics/bdx_hipo
) 2>&1 | tee "$release/build.log"
python3 "$root/farm/release.py" finish "$release"
echo "Use --release $release/release.json when planning jobs."
