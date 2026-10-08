#!/usr/bin/env bash
# Copy private lightweight metadata and SComatic reference resources into a
# durable resource tree. Supply source paths explicitly; nothing is downloaded.
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage: provision_compute1_resources.sh --resource-root PATH --metadata PATH \
  --editing-2020 PATH --pon-2020 PATH --bed-2020 PATH \
  --editing-3 PATH --pon-3 PATH --bed-3 PATH
USAGE
  exit 2
}

root= metadata= editing2020= pon2020= bed2020= editing3= pon3= bed3=
while (($#)); do
  case "$1" in
    --resource-root) root=$2; shift 2;; --metadata) metadata=$2; shift 2;;
    --editing-2020) editing2020=$2; shift 2;; --pon-2020) pon2020=$2; shift 2;; --bed-2020) bed2020=$2; shift 2;;
    --editing-3) editing3=$2; shift 2;; --pon-3) pon3=$2; shift 2;; --bed-3) bed3=$2; shift 2;;
    *) usage;;
  esac
done
for path in "$root" "$metadata" "$editing2020" "$pon2020" "$bed2020" "$editing3" "$pon3" "$bed3"; do [[ -n "$path" ]] || usage; done
for path in "$metadata" "$editing2020" "$pon2020" "$bed2020" "$editing3" "$pon3" "$bed3"; do [[ -s "$path" ]] || { echo "missing source: $path" >&2; exit 1; }; done

copy_once() {
  local source=$1 destination=$2 temporary
  [[ -s "$destination" ]] && return 0
  temporary="${destination}.partial.$$"
  rm -f "$temporary"
  cp "$source" "$temporary"
  mv "$temporary" "$destination"
}

mkdir -p "$root"/metadata "$root"/references/GRCh38_2020_A "$root"/references/GRCh38_3_0_0 "$root"/configs
copy_once "$metadata" "$root/metadata/$(basename "$metadata")"
copy_once "$editing2020" "$root/references/GRCh38_2020_A/$(basename "$editing2020")"
copy_once "$pon2020" "$root/references/GRCh38_2020_A/$(basename "$pon2020")"
copy_once "$bed2020" "$root/references/GRCh38_2020_A/$(basename "$bed2020")"
copy_once "$editing3" "$root/references/GRCh38_3_0_0/$(basename "$editing3")"
copy_once "$pon3" "$root/references/GRCh38_3_0_0/$(basename "$pon3")"
copy_once "$bed3" "$root/references/GRCh38_3_0_0/$(basename "$bed3")"
(cd "$root" && find metadata references -type f -print0 | sort -z | xargs -0 sha256sum) > "$root/configs/resource_sha256sums.txt"
