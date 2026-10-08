#!/usr/bin/env bash
# Internal task wrapper for the paired CSF/PBMC v3 graph.
set -euo pipefail

mode=${1:?mode}; shift
RUNROOT=${RUNROOT:?set RUNROOT}
ENV=${SCOMATIC_ENV:-}
REPO=${SCOMATIC_REPO:-}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
IFS=',' read -r -a TYPES <<< "${CELL_TYPES:-B,CSF_DC,CSF_mac,CSF_mono,NK,PBMC_DC,PBMC_mac,PBMC_mono,Plasma,T,other_myeloid,progenitors}"

die() { echo "ERROR: $*" >&2; exit 1; }

runtime() {
  local root=$1 name=$2
  export HOME="$root/tmp/home_$name"
  export XDG_CACHE_HOME="$root/tmp/xdg_$name"
  export TMPDIR="$root/tmp/$name"
  export PATH="$ENV/bin:$PATH"
  mkdir -p "$HOME" "$XDG_CACHE_HOME" "$TMPDIR"
}

check_environment() {
  [[ -n "$ENV" && -n "$REPO" ]] || die "set SCOMATIC_ENV and SCOMATIC_REPO for SComatic stages"
  [[ -x "$ENV/bin/python" && -x "$ENV/bin/samtools" && -x "$ENV/bin/bedtools" ]] || die "SComatic tools missing under $ENV"
  [[ -s "$REPO/scripts/SplitBam/SplitBamCellTypes.py" ]] || die "SComatic source tree is incomplete: $REPO"
}

load_reference_config() {
  local root=$1
  IFS=$'\t' read -r REF_KEY REF GTF BED EDITING PON < <(tail -n +2 "$root/input/reference_config.tsv")
  [[ -n ${REF:-} && -s "$REF.fai" && -s ${GTF:-} && -s ${BED:-} && -s ${EDITING:-} && -s ${PON:-} ]] || die "invalid reference configuration"
}

make_metadata() {
  local meta=$1 source_id=$2 meta_prefix=$3 tissue=$4 out=$5
  # The public input contract has barcode, tissue, and v3 cell type in columns 1-3.
  # `other` is intentionally excluded while `other_myeloid` is retained.
  if [[ "$meta" == *.gz ]]; then decoder=(gzip -cd -- "$meta"); else decoder=(cat -- "$meta"); fi
  "${decoder[@]}" | awk -F '\t' -v id="$source_id" -v prefix="$meta_prefix" -v tissue="$tissue" '
    BEGIN { OFS="\t"; print "Index", "Cell_type" }
    NR > 1 && index($1, prefix "_") == 1 && $2 == tissue && $3 != "" && $3 != "NA" && $3 != "other" {
      barcode=$1
      sub("^" prefix "_", "", barcode)
      sub(/-1$/, "", barcode)
      print id "_" barcode, $3
      print barcode, $3
    }
  ' > "$out"
}

case "$mode" in
  prepare)
    donor=$1; donor_id=$2; ref_key=$3; manifest=$4; meta=$5; ref_config=$6
    root="$RUNROOT/$donor"
    runtime "$root" prepare; check_environment
    mkdir -p "$root/input" "$root/markers"
    [[ -s "$manifest" && -s "$meta" && -s "$ref_config" ]] || die "missing manifest, metadata, or reference configuration"
    awk -F '\t' -v donor="$donor" 'NR == 1 || $2 == donor' "$manifest" > "$root/input/libraries.tsv"
    [[ $(wc -l < "$root/input/libraries.tsv") -eq 3 ]] || die "each donor must have exactly paired CSF/PBMC rows"
    awk -F '\t' -v key="$ref_key" 'NR == 1 || $1 == key' "$ref_config" > "$root/input/reference_config.tsv"
    [[ $(wc -l < "$root/input/reference_config.tsv") -eq 2 ]] || die "reference_key $ref_key is missing or duplicated"
    csf=$(awk -F '\t' 'NR > 1 && $4 == "CSF" {n++} END {print n+0}' "$root/input/libraries.tsv")
    pbmc=$(awk -F '\t' 'NR > 1 && $4 == "PBMC" {n++} END {print n+0}' "$root/input/libraries.tsv")
    [[ "$csf" -eq 1 && "$pbmc" -eq 1 ]] || die "each donor needs one CSF and one PBMC manifest row"
    while IFS=$'\t' read -r library_donor_id library_donor tissue source_id meta_prefix bam library_ref; do
      [[ "$library_donor" == "$donor" && "$library_donor_id" == "$donor_id" && "$library_ref" == "$ref_key" ]] || die "inconsistent donor manifest"
      make_metadata "$meta" "$source_id" "$meta_prefix" "$tissue" "$root/input/$source_id.pre.tsv"
      [[ $(wc -l < "$root/input/$source_id.pre.tsv") -gt 1 ]] || die "no usable metadata cells for $source_id"
      awk -F '\t' 'NR > 1 && $1 !~ /_/ { print $1 }' "$root/input/$source_id.pre.tsv" | sort -u > "$root/input/$source_id.barcodes.txt"
    done < <(tail -n +2 "$root/input/libraries.tsv")
    cat "$root"/input/*.barcodes.txt | sort | uniq -d > "$root/input/overlapping_plain_barcodes_removed.txt"
    while IFS=$'\t' read -r _ library_donor tissue source_id meta_prefix bam library_ref; do
      awk -F '\t' -v badfile="$root/input/overlapping_plain_barcodes_removed.txt" '
        BEGIN { OFS="\t"; while ((getline barcode < badfile) > 0) bad[barcode] = 1 }
        NR == 1 { print; next }
        { barcode=$1; sub(/^.*_/, "", barcode); if (!(barcode in bad)) print }
      ' "$root/input/$source_id.pre.tsv" > "$root/input/$source_id.cell_types.tsv"
    done < <(tail -n +2 "$root/input/libraries.tsv")
    touch "$donor.prepare.ok"
    ;;

  step1)
    donor=$1; tissue=$2; source_id=$3; bam=$4
    root="$RUNROOT/$donor"; runtime "$root" "step1_$source_id"; check_environment; load_reference_config "$root"
    meta="$root/input/$source_id.cell_types.tsv"; out="$root/results/Step1_BamCellTypes/$source_id"
    [[ -s "$bam" && -s "$bam.bai" && -s "$meta" ]] || die "missing Step 1 BAM, index, or metadata for $source_id"
    mkdir -p "$out" "$root/markers"
    if [[ -f "$root/markers/step1_$source_id.ok" ]] && find "$out" -maxdepth 1 -name '*.bam' -size +0c -print -quit | grep -q .; then
      echo "RECOVERY_REUSE_STEP1 $source_id"; touch "$donor.$source_id.step1.ok"; exit 0
    fi
    "$ENV/bin/samtools" quickcheck "$bam"
    "$ENV/bin/python" "$REPO/scripts/SplitBam/SplitBamCellTypes.py" --bam "$bam" --meta "$meta" --id "$source_id" --n_trim 5 --max_nM 5 --max_NH 1 --outdir "$out"
    while IFS= read -r -d '' f; do "$ENV/bin/samtools" quickcheck "$f"; [[ -s "$f.bai" ]] || "$ENV/bin/samtools" index -@ 1 "$f"; done < <(find "$out" -maxdepth 1 -name '*.bam' -type f -print0)
    find "$out" -maxdepth 1 -name '*.bam' -size +0c -print -quit | grep -q . || die "no Step 1 BAMs for $source_id"
    touch "$root/markers/step1_$source_id.ok" "$donor.$source_id.step1.ok"
    ;;

  joint)
    donor=$1; root="$RUNROOT/$donor"; runtime "$root" joint; check_environment; load_reference_config "$root"
    out="$root/results/Step1JointBams"; mkdir -p "$out" "$root/markers"
    if [[ -f "$root/markers/joint.ok" ]] && find "$out" -maxdepth 1 -name '*.bam' -size +0c -print -quit | grep -q .; then
      echo "RECOVERY_REUSE_JOINT $donor"; touch "$donor.joint.ok"; exit 0
    fi
    printf 'Cell_type\tsource_id\ttissue\tjoint_BAM\n' > "$root/input/joint_celltype_sources.tsv"
    for type in "${TYPES[@]}"; do
      inputs=()
      while IFS=$'\t' read -r _ _ tissue source_id _ _ _; do
        case "$type" in CSF_*) [[ "$tissue" == CSF ]] || continue ;; PBMC_*) [[ "$tissue" == PBMC ]] || continue ;; esac
        candidate="$root/results/Step1_BamCellTypes/$source_id/$source_id.$type.bam"
        if [[ -s "$candidate" ]]; then inputs+=("$candidate"); printf '%s\t%s\t%s\t%s\n' "$type" "$source_id" "$tissue" "$out/$type.bam" >> "$root/input/joint_celltype_sources.tsv"; fi
      done < <(tail -n +2 "$root/input/libraries.tsv")
      ((${#inputs[@]})) || continue
      if ((${#inputs[@]} == 1)); then cp "${inputs[0]}" "$out/$type.bam"; else "$ENV/bin/samtools" merge -f -@ 2 "$out/$type.bam" "${inputs[@]}"; fi
      "$ENV/bin/samtools" index -@ 1 "$out/$type.bam"; "$ENV/bin/samtools" quickcheck "$out/$type.bam"
    done
    find "$out" -maxdepth 1 -name '*.bam' -size +0c -print -quit | grep -q . || die "no joint BAMs"
    touch "$root/markers/joint.ok" "$donor.joint.ok"
    ;;

  step2)
    donor=$1; type=$2; root="$RUNROOT/$donor"; runtime "$root" "step2_$type"; check_environment; load_reference_config "$root"
    bam="$root/results/Step1JointBams/$type.bam"; out="$root/results/Step2_BaseCellCounts"; outfile="$out/$type.tsv"; mkdir -p "$out" "$root/markers"
    if [[ -f "$root/markers/step2_$type.ok" ]] && { [[ ! -s "$bam" ]] || [[ -s "$outfile" ]]; }; then echo "RECOVERY_REUSE_STEP2 $type"; touch "$donor.$type.step2.ok"; exit 0; fi
    if [[ ! -s "$bam" ]]; then printf 'SKIPPED\n' > "$root/markers/step2_$type.ok"; touch "$donor.$type.step2.ok"; exit 0; fi
    "$ENV/bin/python" "$REPO/scripts/BaseCellCounter/BaseCellCounter.py" --bam "$bam" --ref "$REF" --chrom all --out_folder "$out" --min_bq "${MIN_BASE_QUALITY:-30}" --tmp_dir "$TMPDIR" --nprocs "${NPROCS:-10}"
    touch "$root/markers/step2_$type.ok" "$donor.$type.step2.ok"
    ;;

  merge)
    donor=$1; root="$RUNROOT/$donor"; runtime "$root" merge; check_environment; load_reference_config "$root"
    out="$root/results/Step3_BaseCellCountsMerged"; merged="$out/$donor.BaseCellCounts.AllCellTypes.tsv"; mkdir -p "$out" "$root/markers"
    if [[ -f "$root/markers/step3.ok" && -s "$merged" ]]; then echo "RECOVERY_REUSE_STEP3 $donor"; touch "$donor.merge.ok"; exit 0; fi
    "$ENV/bin/python" "$REPO/scripts/MergeCounts/MergeBaseCellCounts.py" --tsv_folder "$root/results/Step2_BaseCellCounts" --outfile "$merged"
    [[ -s "$merged" ]] || die "Step 3 output absent"
    touch "$root/markers/step3.ok" "$donor.merge.ok"
    ;;

  step4)
    donor=$1; root="$RUNROOT/$donor"; runtime "$root" step4; check_environment; load_reference_config "$root"
    step3="$root/results/Step3_BaseCellCountsMerged/$donor.BaseCellCounts.AllCellTypes.tsv"; out="$root/results/Step4_VariantCalling_maxCellTypes${MAX_CELL_TYPES:-15}"
    if [[ -f "$root/markers/postprocess.ok" ]] && find "$out" -maxdepth 1 -name '*.calling.step2.pass.ref_support.annotated.tsv' -size +0c -print -quit | grep -q .; then echo "RECOVERY_REUSE_STEP4 $donor"; touch "$donor.step4.ok"; exit 0; fi
    "$SCRIPT_DIR/run_scomatic_step4_from_step3.sh" --step3-file "$step3" --output-dir "$out" --tag "${donor}_maxCellTypes${MAX_CELL_TYPES:-15}" --repo "$REPO" --env-prefix "$ENV" --reference "$REF" --bed "$BED" --editing "$EDITING" --pon "$PON" --max-cell-types "${MAX_CELL_TYPES:-15}"
    "$SCRIPT_DIR/postprocess_step4.sh" "$root" "$out" "$root/markers/postprocess.ok" "$GTF" "$BED" "$ENV"
    touch "$donor.step4.ok"
    ;;

  gnomad_popmax)
    manifest=$1; [[ -s "$manifest" ]] || die "missing library manifest"
    inputs=(); declare -A seen=()
    while IFS=$'\t' read -r _ donor _ _ _ _ _; do
      [[ ${seen[$donor]:-} ]] && continue; seen[$donor]=1
      table=$(find "$RUNROOT/$donor/results" -type f -name '*.calling.step2.pass.ref_support.annotated.tsv' -print -quit)
      [[ -s "$table" ]] || die "missing annotated PASS table for $donor"; inputs+=("$donor=$table")
    done < <(tail -n +2 "$manifest")
    "$SCRIPT_DIR/gnomad_remote_vcf_filter.sh" --cutoff "${GNOMAD_AF_CUTOFF:-0.001}" "$RUNROOT/gnomad_v4_${GNOMAD_RELEASE:-4.1.1}_remote_vcf_popmax" "${inputs[@]}"
    ;;

  *) die "unknown mode: $mode" ;;
esac
