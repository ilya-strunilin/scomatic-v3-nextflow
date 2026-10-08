#!/usr/bin/env bash
# Create a user-scoped Nextflow/Java and LSF-client runtime for a Compute1
# Docker controller. Run this once from a Compute1 login host, not in a task.
set -euo pipefail
umask 022

usage() {
  echo "Usage: $0 --runtime-root /scratch1/.../nextflow-runtime [--nextflow-version VERSION]" >&2
  exit 2
}

runtime_root=
nextflow_version=${NEXTFLOW_VERSION:-26.04.6}
while (($#)); do
  case "$1" in
    --runtime-root) runtime_root=${2:-}; shift 2 ;;
    --nextflow-version) nextflow_version=${2:-}; shift 2 ;;
    *) usage ;;
  esac
done
[[ "$runtime_root" == /scratch1/* ]] || usage
[[ "$nextflow_version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || usage

mamba_bin=${MAMBA_BIN:-/opt/conda/bin/mamba}
lsf_root=${LSF_ROOT:-/opt/ibm/lsfsuite/lsf/10.1/linux2.6-glibc2.3-x86_64}
env_prefix="$runtime_root/envs/nextflow-$nextflow_version"
mamba_root="$runtime_root/mamba_root"
mamba_pkgs="$runtime_root/mamba_pkgs"
nxf_home="$runtime_root/nxf_home"
bin_dir="$runtime_root/bin"
lsf_client="$runtime_root/lsf_client"

[[ -x "$mamba_bin" ]] || { echo "Mamba is unavailable: $mamba_bin" >&2; exit 3; }
for directory in bin conf etc lib; do
  [[ -d "$lsf_root/$directory" ]] || { echo "LSF directory is unavailable: $lsf_root/$directory" >&2; exit 3; }
  [[ ! -e "$lsf_client/$directory" ]] || {
    echo "Refusing to overwrite existing LSF client directory: $lsf_client/$directory" >&2
    exit 3
  }
done

mkdir -p "$runtime_root"/{envs,mamba_root,mamba_pkgs,nxf_home,bin,home,xdg_cache}
mkdir -p "$lsf_client"
for directory in bin conf etc lib; do cp -a "$lsf_root/$directory" "$lsf_client/$directory"; done

export MAMBA_ROOT_PREFIX="$mamba_root"
export MAMBA_PKGS_DIRS="$mamba_pkgs"
export CONDA_PKGS_DIRS="$mamba_pkgs"
export HOME="$runtime_root/home"
export XDG_CACHE_HOME="$runtime_root/xdg_cache"
export NXF_HOME="$nxf_home"
export NXF_VER="$nextflow_version"

if [[ ! -x "$env_prefix/bin/java" ]]; then
  "$mamba_bin" create --yes --prefix "$env_prefix" --strict-channel-priority \
    --channel conda-forge openjdk=21
fi
find "$env_prefix/bin" -maxdepth 1 -type f -exec chmod u+x {} +
find "$env_prefix/lib/jvm" -type f -name jspawnhelper -exec chmod u+x {} +

export JAVA_HOME="$env_prefix"
export PATH="$env_prefix/bin:$bin_dir:$PATH"
if [[ ! -x "$bin_dir/nextflow" ]]; then
  (cd "$bin_dir" && curl -fsSL https://get.nextflow.io | bash)
fi

cat > "$runtime_root/runtime.env" <<EOF
export JAVA_HOME="$env_prefix"
export PATH="$env_prefix/bin:$bin_dir:\$PATH"
export NXF_HOME="$nxf_home"
export NXF_VER="$nextflow_version"
EOF

"$bin_dir/nextflow" -version
LSF_BINDIR="$lsf_client/bin" \
LSF_ENVDIR="$lsf_client/conf" \
LSF_LIBDIR="$lsf_client/lib" \
LSF_SERVERDIR="$lsf_client/etc" \
PATH="$lsf_client/bin:$PATH" \
"$lsf_client/bin/bsub" -V

touch "$runtime_root"
find "$env_prefix" "$mamba_pkgs" "$nxf_home" "$bin_dir" -exec touch -h {} +
printf 'COMPUTE1_CONTROLLER_RUNTIME_OK\n' > "$runtime_root/COMPUTE1_CONTROLLER_RUNTIME_OK"
