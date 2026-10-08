#!/usr/bin/env bash
# Query the public indexed gnomAD v4 VCFs with bcftools/HTSlib. Candidate
# coordinates are requested from Google Cloud; no gnomAD VCF is downloaded.
set -euo pipefail

cutoff=0.001
if [[ ${1:-} == --cutoff ]]; then cutoff=${2:?cutoff required}; shift 2; fi
out_dir=${1:?output directory required}; shift
(( $# > 0 )) || { echo "at least one DONOR=TABLE input is required" >&2; exit 2; }
command -v bcftools >/dev/null || { echo "bcftools is required" >&2; exit 127; }
[[ "$cutoff" =~ ^0\.[0-9]+$ ]] || { echo "cutoff must be a decimal fraction" >&2; exit 2; }

mkdir -p "$out_dir"/{regions,raw,results}
candidates="$out_dir/results/gnomad_v4_1_1_popmax_candidates.tsv"
awk -F '\t' '
  BEGIN { OFS="\t" }
  FNR == 1 { delete c; for (i=1;i<=NF;i++) { n=$i; sub(/^#/, "", n); c[n]=i }; if (!("CHROM" in c && "Start" in c && "REF" in c && "ALT" in c)) exit 2; next }
  { chrom=$c["CHROM"]; if (chrom !~ /^chr/) chrom=(chrom=="M" || chrom=="MT") ? "chrM" : "chr" chrom; n=split($c["ALT"],alts,","); for(i=1;i<=n;i++) if(alts[i]!="" && alts[i]!=".") print chrom,$c["Start"],toupper($c["REF"]),toupper(alts[i]) }
' "${@#*=}" | sort -k1,1 -k2,2n -k3,3 -k4,4 -u > "$candidates"
[[ -s "$candidates" ]] || { echo "no candidate alleles found" >&2; exit 1; }

queryable="$out_dir/regions/queryable_loci.tsv"
awk -F '\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1,$2,$2}' OFS='\t' "$candidates" | sort -u > "$queryable"
printf 'source\tchrom\turl\trecords\n' > "$out_dir/results/query_manifest.tsv"
while IFS=$'\t' read -r chrom _ _; do echo "$chrom"; done < "$queryable" | sort -u | while read -r chrom; do
  region="$out_dir/regions/$chrom.tsv"; awk -F '\t' -v c="$chrom" '$1==c {print $1,$2,$3}' OFS='\t' "$queryable" > "$region"
  for source in exomes genomes; do
    url="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/${source}/gnomad.${source}.v4.1.1.sites.${chrom}.vcf.bgz"
    raw="$out_dir/raw/${source}.${chrom}.tsv"
    bcftools query -R "$region" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/AF_grpmax\t%INFO/grpmax\t%INFO/AF\t%INFO\n' "$url" > "$raw"
    printf '%s\t%s\t%s\t%s\n' "$source" "$chrom" "$url" "$(wc -l < "$raw")" >> "$out_dir/results/query_manifest.tsv"
  done
done

python3 - "$candidates" "$out_dir" "$cutoff" "$@" <<'PY'
import csv, glob, os, sys
from collections import defaultdict

candidate_path, out_dir, cutoff_text, *specs = sys.argv[1:]
cutoff = float(cutoff_text)
na = 'NA'

def present(x): return x not in ('', '.', 'NA', None)
def clean(x): return x if present(x) else na
def num(x): return float(x) if present(x) else None
def key(chrom, pos, ref, alt):
    chrom = chrom if chrom.startswith('chr') else ('chrM' if chrom in ('M', 'MT') else 'chr' + chrom)
    return chrom, str(pos), ref.upper(), alt.upper()

candidates = set()
with open(candidate_path) as handle:
    for chrom, pos, ref, alt in csv.reader(handle, delimiter='\t'):
        candidates.add(key(chrom, pos, ref, alt))

records = {source: {} for source in ('exomes', 'genomes')}
for source in records:
    for path in glob.glob(os.path.join(out_dir, 'raw', f'{source}.chr*.tsv')):
        with open(path) as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) < 9: continue
                chrom, pos, ref, alts, ident, grpafs, groups, afs, info = row
                alts, grpafs, groups, afs = alts.split(','), grpafs.split(','), groups.split(','), afs.split(',')
                for i, alt in enumerate(alts):
                    k = key(chrom, pos, ref, alt)
                    if k in candidates:
                        records[source][k] = {'id': ident, 'grpaf': grpafs[i] if i < len(grpafs) else na, 'group': groups[i] if i < len(groups) else na, 'af': afs[i] if i < len(afs) else na, 'info': info}

header = ['CHROM','Start','REF','ALT','gnomAD_v4_1_1_variant_id','gnomAD_v4_1_1_rsid','gnomAD_v4_1_1_exomes_AF_grpmax','gnomAD_v4_1_1_exomes_grpmax_ancestry','gnomAD_v4_1_1_exomes_AF','gnomAD_v4_1_1_genomes_AF_grpmax','gnomAD_v4_1_1_genomes_grpmax_ancestry','gnomAD_v4_1_1_genomes_AF','gnomAD_v4_1_1_popmax_AF','gnomAD_v4_1_1_popmax_ancestry','gnomAD_v4_1_1_generic_AF','gnomAD_v4_1_1_generic_AF_source','gnomAD_v4_1_1_match_status','gnomAD_common_variant_filter','gnomAD_v4_1_1_exomes_INFO','gnomAD_v4_1_1_genomes_INFO']

def annotation(k):
    ex, gn = records['exomes'].get(k), records['genomes'].get(k)
    if ex and gn: status = 'IN_GNOMAD_EXOMES_AND_GENOMES'
    elif ex: status = 'IN_GNOMAD_EXOMES'
    elif gn: status = 'IN_GNOMAD_GENOMES'
    else: status = 'NOT_IN_GNOMAD'
    ex = ex or {}; gn = gn or {}
    pop_entries = [(num(ex.get('grpaf')), ex.get('group'), 'exomes'), (num(gn.get('grpaf')), gn.get('group'), 'genomes')]
    generic_entries = [(num(ex.get('af')), 'exomes'), (num(gn.get('af')), 'genomes')]
    pop_entries = [x for x in pop_entries if x[0] is not None]; generic_entries = [x for x in generic_entries if x[0] is not None]
    popmax, ancestry = (max(pop_entries, key=lambda x:x[0])[0:2] if pop_entries else (None, na))
    generic, generic_source = (max(generic_entries, key=lambda x:x[0]) if generic_entries else (None, na))
    if popmax is not None: decision = 'FILTER_POPMAX_AF_GE_' + cutoff_text if popmax >= cutoff else 'RETAIN_POPMAX_AF_LT_' + cutoff_text
    elif generic is not None: decision = 'FILTER_GENERIC_AF_GE_' + cutoff_text if generic >= cutoff else 'RETAIN_GENERIC_AF_LT_' + cutoff_text
    elif status == 'NOT_IN_GNOMAD': decision = 'RETAIN_NOT_IN_GNOMAD'
    else: decision = 'RETAIN_GNOMAD_NO_AF_OR_GRPMAX'
    variant_id = na if status == 'NOT_IN_GNOMAD' else f'{k[0][3:]}-{k[1]}-{k[2]}-{k[3]}'
    rsid = ex.get('id') or gn.get('id') or na
    return [variant_id, clean(rsid), clean(ex.get('grpaf')), clean(ex.get('group')), clean(ex.get('af')), clean(gn.get('grpaf')), clean(gn.get('group')), clean(gn.get('af')), clean(popmax), clean(ancestry), clean(generic), clean(generic_source), status, decision, clean(ex.get('info')), clean(gn.get('info'))]

cache_path = os.path.join(out_dir, 'results', 'gnomad_v4_1_1_popmax_cache.tsv')
annotations = {k: annotation(k) for k in candidates}
with open(cache_path, 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n'); writer.writerow(header)
    for k in sorted(candidates): writer.writerow(list(k) + annotations[k])

summary_path = os.path.join(out_dir, 'results', f'gnomad_v4_1_1_filtered_{cutoff_text}_summary.tsv')
with open(summary_path, 'w', newline='') as summary:
    sw = csv.writer(summary, delimiter='\t', lineterminator='\n'); sw.writerow(['donor','input_pass_variants','retained_total','filtered_popmax','filtered_generic','retained_popmax','retained_generic','retained_unknown'])
    for spec in specs:
        donor, path = spec.split('=', 1)
        annotated_path = os.path.join(out_dir, 'results', f'{donor}.calling.step2.pass.ref_support.annotated.gnomad_popmax.tsv')
        filtered_path = os.path.join(out_dir, 'results', f'{donor}.calling.step2.pass.ref_support.annotated.gnomad_filtered.tsv')
        counts = defaultdict(int); input_n = retained_n = 0
        with open(path, newline='') as source, open(annotated_path, 'w', newline='') as annotated, open(filtered_path, 'w', newline='') as filtered:
            reader = csv.DictReader(source, delimiter='\t'); fieldnames = reader.fieldnames or []
            for required in ('CHROM','Start','REF','ALT'):
                if required not in fieldnames: raise RuntimeError(f'{required} missing from {path}')
            names = fieldnames + header[4:]
            aw = csv.DictWriter(annotated, fieldnames=names, delimiter='\t', lineterminator='\n'); fw = csv.DictWriter(filtered, fieldnames=names, delimiter='\t', lineterminator='\n'); aw.writeheader(); fw.writeheader()
            for row in reader:
                input_n += 1; ordered_alts = []
                for alt in row['ALT'].split(','):
                    if alt and alt not in ordered_alts: ordered_alts.append(alt)
                anns = [annotations[key(row['CHROM'], row['Start'], row['REF'], alt)] for alt in ordered_alts]
                appended = [','.join(values[i] for values in anns) if len(anns) > 1 else anns[0][i] for i in range(len(header)-4)]
                record = dict(row); record.update(dict(zip(header[4:], appended))); aw.writerow(record)
                decisions = [values[13] for values in anns]
                common = any(x.startswith('FILTER_') for x in decisions)
                for x in decisions: counts[x] += 1
                if not common: fw.writerow(record); retained_n += 1
        sw.writerow([donor, input_n, retained_n, counts['FILTER_POPMAX_AF_GE_'+cutoff_text], counts['FILTER_GENERIC_AF_GE_'+cutoff_text], counts['RETAIN_POPMAX_AF_LT_'+cutoff_text], counts['RETAIN_GENERIC_AF_LT_'+cutoff_text], counts['RETAIN_NOT_IN_GNOMAD'] + counts['RETAIN_GNOMAD_NO_AF_OR_GRPMAX']])
PY

"$(dirname "$0")/normalize_tsv_missing_values.sh" "$out_dir"/results/*.gnomad_popmax.tsv "$out_dir"/results/*.gnomad_filtered.tsv
touch "$out_dir/GNOMAD_V4_1_1_REMOTE_VCF_POPMAX_GENERIC_FALLBACK_FULL_INFO_OK"
