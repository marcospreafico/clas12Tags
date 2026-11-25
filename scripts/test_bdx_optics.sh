#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
image=${BDX_GEMC_IMAGE:-jeffersonlab/gemc:dev-almalinux94@sha256:4dd5fa8c6714a7f9241106c07edb32f6c3e2e645ec55958b508b8f0d5f95ac74}
python3 -m unittest discover -s "$repo/tests/bdx_optics" -p 'test_*.py' -v
mkdir -p "$repo/build/bdx_optics/transport_tests"
docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
    -v "$repo:/work/clas12Tags" -w /work/clas12Tags "$image" -lc '
        set -euo pipefail
        module load hipo/4.5.0
        b=build/bdx_optics/transport_tests
        read -r -a flags <<< "$(geant4-config --cflags)"
        read -r -a libs <<< "$(geant4-config --libs)"
        g++ "${flags[@]}" -std=c++17 -O2 -Isource/src -I"$HIPO/include" tests/bdx_optics/transport_test.cc \
            source/src/BdxCrystalOptics.cc -o "$b/transport_test" "${libs[@]}" -L"$HIPO/lib" -lhipo4 -llz4 -lfmt -pthread > "$b/compile.log" 2>&1 || {
            tail -60 "$b/compile.log"; exit 1;
        }
        g++ "${flags[@]}" -std=c++17 -O2 -Isource/src -I"$HIPO/include" tests/bdx_optics/relocation_test.cc \
            source/src/BdxCrystalOptics.cc -o "$b/relocation_test" "${libs[@]}" -L"$HIPO/lib" -lhipo4 -llz4 -lfmt -pthread > "$b/compile-relocation.log" 2>&1 || {
            tail -60 "$b/compile-relocation.log"; exit 1;
        }
        g++ -std=c++17 -O2 -Isource/src -I"$HIPO/include" tests/bdx_optics/hipo_roundtrip_test.cc \
            -L"$HIPO/lib" -lhipo4 -llz4 -lfmt -pthread -o "$b/hipo_roundtrip_test"
        "$b/hipo_roundtrip_test" "$b/large_event.hipo"
        build/bdx_optics/bdx_hipo "$b/large_event.hipo" --check "$b/large_event.json"
        "$b/relocation_test" /work/clas12Tags | tee "$b/relocation_checks.txt"
        python3 tests/bdx_optics/run_transport_checks.py
    '

python3 "$repo/tests/bdx_optics/check_postselection.py"
