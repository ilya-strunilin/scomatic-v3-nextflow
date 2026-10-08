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

base=$(cd "$(dirname "$0")/.." && pwd)
mapfile -t submission_config < <(python3 - "$params" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    params = json.load(handle)

print(params.get("run_root", ""))
print(params.get("lsf_queue") or "general")
print(params.get("lsf_group") or "")
print(params.get("lsf_container") or "strunyadocker/mamba_minimal")
print(params.get("controller_runtime") or "")
print(params.get("lsf_docker_volumes") or "")
PY
)
run_root=${submission_config[0]:-}
params_queue=${submission_config[1]:-general}
params_group=${submission_config[2]:-}
params_container=${submission_config[3]:-strunyadocker/mamba_minimal}
controller_runtime=${submission_config[4]:-}
params_docker_volumes=${submission_config[5]:-}
[[ -n "$run_root" ]] || { echo "run_root is required in params JSON" >&2; exit 2; }
[[ -n "$controller_runtime" && -x "$controller_runtime/bin/nextflow" && -s "$controller_runtime/runtime.env" ]] || {
  echo "controller_runtime must contain executable bin/nextflow and runtime.env" >&2
  exit 2
}
mkdir -p "$run_root"/{logs,controller_work,submission}

# Environment overrides are useful for a one-off site policy, while the JSON
# keeps the normal submission reproducible and visible in the run inputs.
queue=${LSF_QUEUE:-$params_queue}
group=${LSF_GROUP:-$params_group}
container=${LSF_CONTROLLER_CONTAINER:-$params_container}
docker_volumes=${LSF_DOCKER_VOLUMES:-$params_docker_volumes}
[[ -n "$docker_volumes" ]] || { echo "lsf_docker_volumes or LSF_DOCKER_VOLUMES is required" >&2; exit 2; }
export LSF_DOCKER_VOLUMES="$docker_volumes"
export LSF_DOCKER_PRESERVE_ENVIRONMENT=false
submission_args=(-J nf-scomatic-v3 -q "$queue" -a "docker(${container})")
[[ -n "$group" ]] && submission_args+=(-G "$group")

# If Compute1 Docker needs bind mounts, export LSF_DOCKER_VOLUMES before calling
# this script. The paths are site/user-specific and intentionally not tracked.
bsub "${submission_args[@]}" -n 2 -R 'rusage[mem=8GB]' -M 8GB \
  -oo "$run_root/logs/nextflow_controller_%J.log" \
  /bin/bash "$base/bin/run_controller.sh" "$base" "$params" "$resume" \
  | tee "$run_root/submission/controller_submission.log"
