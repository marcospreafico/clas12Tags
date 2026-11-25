#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
image=${BDX_GEMC_IMAGE:-jeffersonlab/gemc:dev-almalinux94@sha256:4dd5fa8c6714a7f9241106c07edb32f6c3e2e645ec55958b508b8f0d5f95ac74}
mkdir -p "$repo/build/bdx_optics"
docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
    -v "$repo:/work/clas12Tags" -w /work/clas12Tags/source \
    "$image" -lc '
        module load ccdb/1.99.6 hipo/4.5.0
        export GEMC=/work/clas12Tags GEMC_DATA_DIR=/work/clas12Tags BDX_NO_EVIO=1
        set -euo pipefail
        scons -j"${1:-4}"
        g++ -std=c++17 -O2 -I"$HIPO/include" ../scripts/bdx_hipo.cc \
            -L"$HIPO/lib" -lhipo4 -llz4 -lfmt -pthread -o ../build/bdx_optics/bdx_hipo
    ' bdx-build "${1:-4}" > "$repo/build/bdx_optics/compile.log" 2>&1 || {
        tail -80 "$repo/build/bdx_optics/compile.log" >&2; exit 1;
    }
echo "GEMC BDX compilato: $repo/source/gemc"
