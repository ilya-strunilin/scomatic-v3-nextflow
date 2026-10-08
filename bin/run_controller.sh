#!/usr/bin/env bash
set -euo pipefail

base=${1:?repository path required}
params=${2:?params JSON required}
resume=${3:-}
mapfile -t controller_config < <(python3 - "$params" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    params = json.load(handle)

print(params.get('run_root', ''))
print(params.get('controller_runtime', ''))
print(params.get('lsf_docker_volumes', ''))
PY
)
run_root=${controller_config[0]:-}
controller_runtime=${controller_config[1]:-}
params_docker_volumes=${controller_config[2]:-}
[[ -n "$run_root" ]] || { echo "run_root is required in params JSON" >&2; exit 2; }
[[ -x "$controller_runtime/bin/nextflow" && -s "$controller_runtime/runtime.env" ]] || {
  echo "controller_runtime is unavailable inside the controller container" >&2
  exit 2
}

# The Compute1 Docker image does not provide a Nextflow or LSF client runtime.
# The user-created runtime is mounted into the container and contains both.
source "$controller_runtime/runtime.env"
export LSF_DOCKER_VOLUMES="${LSF_DOCKER_VOLUMES:-$params_docker_volumes}"
[[ -n "$LSF_DOCKER_VOLUMES" ]] || { echo 'lsf_docker_volumes is required inside the controller' >&2; exit 3; }
export LSF_DOCKER_PRESERVE_ENVIRONMENT=false
export LSF_BINDIR="$controller_runtime/lsf_client/bin"
export LSF_ENVDIR="$controller_runtime/lsf_client/conf"
export LSF_LIBDIR="$controller_runtime/lsf_client/lib"
export LSF_SERVERDIR="$controller_runtime/lsf_client/etc"
export EGO_CONFDIR="$LSF_ENVDIR/ego/compute1-lsf/kernel"
export EGO_ESRVDIR="$LSF_ENVDIR/ego/compute1-lsf/eservice"
export EGO_SEC_CONF="$EGO_CONFDIR"
export PATH="$LSF_BINDIR:$PATH"
for lsf_command in bsub bjobs bkill; do
  [[ -x "$LSF_BINDIR/$lsf_command" ]] || {
    echo "Controller runtime lacks LSF command: $LSF_BINDIR/$lsf_command" >&2
    exit 3
  }
done
export NXF_HOME="$run_root/.nextflow"
export NXF_WORK="$run_root/controller_work"
mkdir -p "$NXF_HOME" "$NXF_WORK"
cd "$base"
if [[ "$resume" == --resume ]]; then
  exec "$controller_runtime/bin/nextflow" run main.nf -profile compute1 -params-file "$params" -resume
else
  exec "$controller_runtime/bin/nextflow" run main.nf -profile compute1 -params-file "$params"
fi
