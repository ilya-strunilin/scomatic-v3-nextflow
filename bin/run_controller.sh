#!/usr/bin/env bash
set -euo pipefail

base=${1:?repository path required}
params=${2:?params JSON required}
resume=${3:-}
run_root=$(python3 - "$params" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))['run_root'])
PY
)
export NXF_HOME="$run_root/.nextflow"
export NXF_WORK="$run_root/controller_work"
mkdir -p "$NXF_HOME" "$NXF_WORK"
cd "$base"
if [[ "$resume" == --resume ]]; then
  exec nextflow run main.nf -profile compute1 -params-file "$params" -resume
else
  exec nextflow run main.nf -profile compute1 -params-file "$params"
fi
