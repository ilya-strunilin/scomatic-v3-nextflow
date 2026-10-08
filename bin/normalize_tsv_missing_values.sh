#!/usr/bin/env bash
set -euo pipefail
(( $# > 0 )) || { echo "usage: $0 TABLE [TABLE ...]" >&2; exit 2; }
for table in "$@"; do
  [[ -s "$table" ]] || { echo "missing or empty table: $table" >&2; exit 2; }
  directory=$(dirname "$table"); name=$(basename "$table"); temporary=$(mktemp "$directory/.${name}.XXXXXX")
  awk -F '\t' -v OFS='\t' '{ for (i = 1; i <= NF; i++) if ($i == "" || $i == ".") $i = "NA"; print }' "$table" > "$temporary"
  mv "$temporary" "$table"
done
