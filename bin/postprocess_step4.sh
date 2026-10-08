#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?run root required}; STEP4=${2:?Step 4 directory required}; MARKER=${3:?marker path required}; GTF=${4:?reference GTF required}; BED=${5:?high-quality BED required}; ENV_PREFIX=${6:?SComatic environment required}
export HOME="$ROOT/tmp/home_postprocess" XDG_CACHE_HOME="$ROOT/tmp/cache_postprocess" TMPDIR="$ROOT/tmp/postprocess" PATH="$ENV_PREFIX/bin:$PATH"
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$TMPDIR"
[[ -x "$ENV_PREFIX/bin/bedtools" && -s "$GTF" && -s "$BED" ]] || { echo "Missing reference or tool" >&2; exit 1; }
pass=$(find "$STEP4" -maxdepth 1 -type f -name '*.calling.step2.pass.tsv' -print -quit)
[[ -n "$pass" ]] || { echo "No Step 4 PASS table found" >&2; exit 1; }
base=${pass%.calling.step2.pass.tsv}; clean=${base}.calling.step2.pass.annotated.tsv; support=${base}.calling.step2.pass.ref_support.tsv; final=${base}.calling.step2.pass.ref_support.annotated.tsv
tmp=$(mktemp -d "$TMPDIR/work.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
awk 'BEGIN{OFS="\t"} /^#CHROM/{if(!h++){sub(/^#/,""); print} next} $0 !~ /^#/ && $6=="PASS"{print}' "$pass" > "$tmp/clean.tsv"
awk 'BEGIN{OFS="\t"} NR==1{next} {print $1,$2-1,$3,NR-1}' "$tmp/clean.tsv" > "$tmp/mutations.bed"
awk 'BEGIN{OFS="\t";FS="\t"} function attr(s,key,n,a,i){n=split(s,a,";");for(i=1;i<=n;i++){gsub(/^[ \t]+|[ \t]+$/, "", a[i]);if(a[i] ~ ("^" key "[ \t]+")){sub("^" key "[ \t]+", "", a[i]);gsub(/^\"|\"$/, "", a[i]);return a[i]}}return "NA"} $3=="gene"{print $1,$4-1,$5,attr($9,"gene_id"),attr($9,"gene_name")}' "$GTF" > "$tmp/genes.bed"
sort -k1,1 -k2,2n -k3,3n "$tmp/mutations.bed" > "$tmp/mutations.sorted.bed"; sort -k1,1 -k2,2n -k3,3n "$tmp/genes.bed" > "$tmp/genes.sorted.bed"
"$ENV_PREFIX/bin/bedtools" closest -a "$tmp/mutations.sorted.bed" -b "$tmp/genes.sorted.bed" -d -t first > "$tmp/closest.tsv"
awk 'BEGIN{OFS="\t"}{print $4,$8,$9,$10}' "$tmp/closest.tsv" | sort -k1,1n > "$tmp/annotations.tsv"
awk 'BEGIN{OFS="\t"} FNR==NR{gene[$1]=$2;symbol[$1]=$3;distance[$1]=$4;next} FNR==1{print $0,"NearestGene","NearestGeneName","NearestGeneDistance";next}{id=FNR-1;print $0,gene[id],symbol[id],distance[id]}' "$tmp/annotations.tsv" "$tmp/clean.tsv" > "$clean"

"$ENV_PREFIX/bin/python" - "$clean" "$final" <<'PY'
import csv, sys
src, dst = sys.argv[1:]
new = ["N_REF_Cell_types","REF_Cell_types","N_REF_VAF_1_Cell_types","REF_VAF_1_Cell_types","N_no_ALT_support_Cell_types","no_ALT_support_Cell_types"]
alleles = {"A": 1, "C": 2, "T": 3, "G": 4}
def integer(x):
    try: return int(x)
    except (TypeError, ValueError): return 0
def parsed(x):
    if x == "NA": return None
    p = x.split("|")
    return (integer(p[0]), integer(p[1]), p[3].split(":"), p[4].split(":")) if len(p) >= 5 else None
with open(src, newline="") as f, open(dst, "w", newline="") as g:
    reader = csv.reader(f, delimiter="\t"); writer = csv.writer(g, delimiter="\t", lineterminator="\n")
    header = next(reader); names = [x.lstrip("#") for x in header]; info, gene = names.index("INFO"), names.index("NearestGene")
    cells = list(range(info + 1, gene)); writer.writerow(header[:gene] + new + header[gene:])
    for row in reader:
        ref, alt = row[names.index("REF")].upper(), row[names.index("ALT")].split(",", 1)[0].upper(); rn, an = alleles.get(ref, 0), alleles.get(alt, 0)
        ref_types, ref_vaf1, no_alt = [], [], []
        for i in cells:
            q = parsed(row[i] if i < len(row) else "")
            if q is None: continue
            dp, nc, cc, bc = q; rbc = integer(bc[rn-1]) if rn and len(bc) >= rn else 0; rcc = integer(cc[rn-1]) if rn and len(cc) >= rn else 0; abc = integer(bc[an-1]) if an and len(bc) >= an else 0; acc = integer(cc[an-1]) if an and len(cc) >= an else 0; ct = header[i].lstrip("#")
            if dp >= 5 and nc >= 5 and rbc >= 3 and rcc >= 2:
                ref_types.append(ct)
                if sum(integer(v) for j, v in enumerate(bc) if j != rn - 1) == 0: ref_vaf1.append(ct)
            if dp >= 5 and nc >= 5 and not (abc >= 3 and acc >= 2): no_alt.append(ct)
        writer.writerow(row[:gene] + [str(len(ref_types)), ",".join(ref_types), str(len(ref_vaf1)), ",".join(ref_vaf1), str(len(no_alt)), ",".join(no_alt)] + row[gene:])
PY
awk 'BEGIN{OFS="\t"} FNR==1{for(i=1;i<=NF;i++) if($i=="NearestGene"){g=i;break}; for(i=1;i<g;i++) printf "%s%s",$i,(i==g-1?ORS:OFS);next}{for(i=1;i<g;i++) printf "%s%s",$i,(i==g-1?ORS:OFS)}' "$final" > "$support"
printf 'SCOMATIC_JOINT_POSTPROCESS_OK\n' > "$MARKER"
