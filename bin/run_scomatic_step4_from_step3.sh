#!/usr/bin/env bash
set -euo pipefail

usage() { echo "Usage: $0 --step3-file PATH --output-dir PATH --tag TAG --repo PATH --env-prefix PATH --reference PATH --bed PATH --editing PATH --pon PATH --max-cell-types N" >&2; exit 2; }
STEP3_FILE= OUTPUT_DIR= TAG= REPO= ENV_PREFIX= REFERENCE= BED= EDITING= PON= MAX_CELL_TYPES=
while (($#)); do
  case "$1" in
    --step3-file) STEP3_FILE=$2; shift 2;; --output-dir) OUTPUT_DIR=$2; shift 2;; --tag) TAG=$2; shift 2;;
    --repo) REPO=$2; shift 2;; --env-prefix) ENV_PREFIX=$2; shift 2;; --reference) REFERENCE=$2; shift 2;;
    --bed) BED=$2; shift 2;; --editing) EDITING=$2; shift 2;; --pon) PON=$2; shift 2;; --max-cell-types) MAX_CELL_TYPES=$2; shift 2;;
    *) usage;;
  esac
done
for required in "$STEP3_FILE" "$OUTPUT_DIR" "$TAG" "$REPO" "$ENV_PREFIX" "$REFERENCE" "$BED" "$EDITING" "$PON" "$MAX_CELL_TYPES"; do [[ -n "$required" && -e "$required" ]] || { [[ "$required" == "$OUTPUT_DIR" || "$required" == "$TAG" || "$required" == "$MAX_CELL_TYPES" ]] || { echo "Missing required path: $required" >&2; exit 1; }; }; done
export HOME="$OUTPUT_DIR/../tmp/home_step4" XDG_CACHE_HOME="$OUTPUT_DIR/../tmp/xdg_step4" TMPDIR="$OUTPUT_DIR/../tmp/tmp_step4" PATH="$ENV_PREFIX/bin:$PATH"
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$TMPDIR" "$OUTPUT_DIR"
"$ENV_PREFIX/bin/python" "$REPO/scripts/BaseCellCalling/BaseCellCalling.step1.py" --infile "$STEP3_FILE" --outfile "$OUTPUT_DIR/$TAG" --ref "$REFERENCE" --max_cell_types "$MAX_CELL_TYPES"
"$ENV_PREFIX/bin/python" "$REPO/scripts/BaseCellCalling/BaseCellCalling.step2.py" --infile "$OUTPUT_DIR/$TAG.calling.step1.tsv" --outfile "$OUTPUT_DIR/$TAG" --editing "$EDITING" --pon "$PON"
"$ENV_PREFIX/bin/bedtools" intersect -header -a "$OUTPUT_DIR/$TAG.calling.step2.tsv" -b "$BED" | awk '$1 ~ /^#/ || $6 == "PASS"' > "$OUTPUT_DIR/$TAG.calling.step2.pass.tsv"
