#!/usr/bin/env bash
# These three staged inputs are frozen by the campaign plan.
set -euo pipefail
[[ $# == 2 ]] || exit 2
printf '%s  runtime.tar.gz\n' "$1" | sha256sum -c -
printf '%s  job.json\n' "$2" | sha256sum -c -
tar -xzf runtime.tar.gz
export BDX_RUNTIME_ROOT="$PWD/runtime"
source runtime/farm/environment.sh runtime/farm/ifarm.env
exec python3 runtime/farm/worker.py "$PWD/job.json"
