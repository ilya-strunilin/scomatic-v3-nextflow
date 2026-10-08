#!/usr/bin/env bash
# Submit the lightweight Nextflow controller to Compute1 LSF.
set -euo pipefail

usage() { echo "Usage: $0 --params PATH [--resume]" >&2; exit 2; }
params=; resume=
while (($#)); do
  case "$1" in --params) params=$2; shift 2;; --resume) resume='-resume'; shift;; *) usage;; esac
done
[[ -s "$params" ]] || usage
command -v bsub >/dev/null || { echo "bsub is required" >&2; exit 127; }
command -v nextflow >/dev/null || { echo "nextflow is required on the controller host" >&2; exit 127; }

base=$(cd "$(dirname "$0")/.." && pwd)
run_root=$(python3 - "$params" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))['run_root'])
PY
)
[[ -n "$run_root" ]] || { echo "run_root is required in params JSON" >&2; exit 2; }
mkdir -p "$run_root"/{logs,controller_work,submission}

# If Compute1 Docker needs bind mounts, export LSF_DOCKER_VOLUMES before calling
# this script. The paths are site/user-specific and intentionally not tracked.
bsub -J nf-scomatic-v3 -q "${LSF_QUEUE:-general}" -n 2 -R 'rusage[mem=8GB]' -M 8GB \
  -oo "$run_root/logs/nextflow_controller_%J.log" \
  /bin/bash "$base/bin/run_controller.sh" "$base" "$params" "$resume" \
  | tee "$run_root/submission/controller_submission.log"
