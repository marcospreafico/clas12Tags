#!/usr/bin/env bash
# Run this with bash, not source: site modules stay in a child shell.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd -P)
if [[ ${1:-} == --help ]]; then
    echo 'Usage: bash farm/test_bdx.sh [run_bdx_optics.py options]'
    echo 'Defaults: types 1 and 8, distances 20 and 200 mm, 2 events/point, 100 photons/MeV.'
    echo 'Select a completed release with BDX_RELEASE_JSON, otherwise the latest is used.'
    exit 0
fi
BDX_RELEASE_JSON=$(python3 - "$root" "${BDX_RELEASE_JSON:-}" <<'PY'
from pathlib import Path
import sys
if sys.argv[2]:
    path=Path(sys.argv[2]).resolve()
    if not path.is_file(): raise SystemExit('BDX_RELEASE_JSON does not exist')
else:
    candidates=sorted((Path(sys.argv[1])/'farm/build').glob('*/release.json'))
    if not candidates: raise SystemExit('No completed release: run bash farm/build_bdx.sh first')
    path=candidates[-1].resolve()
print(path)
PY
)
export BDX_RUNTIME_ROOT="$(dirname "$BDX_RELEASE_JSON")/runtime"
source "$root/farm/environment.sh" "$BDX_RUNTIME_ROOT/farm/ifarm.env"
printf 'Testing native release: %s\n' "$BDX_RELEASE_JSON"
exec python3 "$BDX_RUNTIME_ROOT/scripts/run_bdx_optics.py" --native \
    --types 1 8 --distances-mm 20 200 --events 2 --set light_yield_per_MeV=100 "$@"
