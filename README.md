# SComatic v3 Nextflow

This is a parameterized, paired-CSF/PBMC SComatic workflow for Compute1. It
splits each source BAM by cell type, creates tissue-aware joint BAMs, runs
SComatic Steps 2–4, adds nearest-gene and reference-support annotations, and
filters common germline variants using public gnomAD v4.1.1 VCF queries through
HTSlib/bcftools.

The published defaults reproduce the v3 configuration used here:

- paired CSF and PBMC source BAMs are split separately, then joint BAMs are
  made per tissue-aware cell type;
- `CSF_mac`, `CSF_mono`, `CSF_DC`, `PBMC_mac`, `PBMC_mono`, and `PBMC_DC` stay
  separate;
- only the exact cell label `other` is excluded; `other_myeloid` is retained;
- BaseCellCounter uses minimum base quality 30;
- Step 4 uses `max_cell_types = 15`;
- gnomAD excludes variants with `AF_grpmax >= 0.001` (0.1%). When gnomAD has
  no `AF_grpmax` for an allele, the larger cohort-wide `AF` from exomes/genomes
  is used as a fallback.

Project-specific clinical tables, metadata, manifests, BAMs, run outputs,
credentials, and LSF client files are deliberately not in this repository.

## Install SComatic and gnomAD query dependencies

Install the SComatic checkout in durable storage and its environment in a
project-scoped Scratch1 path. The environment must be available to every LSF
task through the Docker bind mounts described below. Do not install packages at
task runtime.

SComatic's published requirements use Python 3.7, R 3.6.1 or later,
`samtools`, `bedtools`, and `datamash`. This pipeline additionally requires
`bcftools` (and the HTSlib it brings) for gnomAD queries.

```bash
# Choose paths owned by your Compute1 account/project.
export PROJECT_CODE=/storage3/fs1/<group>/Active/<user>/Projects/<project>/scomatic-v3-nextflow
export SCOMATIC_REPO=/storage3/fs1/<group>/Active/<user>/Projects/<project>/SComatic
export PROJECT_SCRATCH=/scratch1/fs1/<group>/<user>/codex/scomatic-v3-runtime
export SCOMATIC_ENV=$PROJECT_SCRATCH/envs/scomatic
export MAMBA_ROOT_PREFIX=$PROJECT_SCRATCH/mamba_root
export MAMBA_PKGS_DIRS=$PROJECT_SCRATCH/mamba_pkgs
export SCOMATIC_REV=5fd34e6d6b52524a0e392eeee6cb2f86d497863e

mkdir -p "$PROJECT_SCRATCH"/{envs,mamba_root,mamba_pkgs}
git clone https://github.com/cortes-ciriano-lab/SComatic.git "$SCOMATIC_REPO"
git -C "$SCOMATIC_REPO" checkout "$SCOMATIC_REV"

# This is the exact upstream SComatic commit used by the v3 analyses that this
# workflow configuration reproduces. Record it with every analysis.
git -C "$SCOMATIC_REPO" rev-parse HEAD

mamba create --yes --prefix "$SCOMATIC_ENV" --strict-channel-priority \
  --channel conda-forge --channel bioconda \
  python=3.7 r-base=3.6.1 samtools bedtools datamash bcftools
"$SCOMATIC_ENV/bin/python" -m pip install --no-cache-dir \
  -r "$SCOMATIC_REPO/requirements.txt"
"$SCOMATIC_ENV/bin/Rscript" "$SCOMATIC_REPO/r_requirements_install.R"
```

Validate the runtime before launching a cohort:

```bash
"$SCOMATIC_ENV/bin/python" - <<'PY'
import numpy, numpy_groupies, pandas, pybedtools, pysam, rpy2, scipy
print("SComatic Python dependencies: OK")
PY
"$SCOMATIC_ENV/bin/samtools" --version | head -1
"$SCOMATIC_ENV/bin/bedtools" --version
"$SCOMATIC_ENV/bin/bcftools" --version | head -1
```

`bcftools` performs the gnomAD part of the workflow. No local gnomAD VCF or
Hail installation is needed: `bin/gnomad_remote_vcf_filter.sh` queries the
public, indexed gnomAD v4.1.1 exome and genome VCFs over HTTPS, keeps a local
allele cache in the run output, and filters on the configured AF threshold.
Therefore, the task container needs outbound HTTPS access to public Google
Cloud Storage. The first real gnomAD task is also the practical connectivity
check. The workflow does not upload BAMs or metadata to gnomAD.

SComatic's upstream instructions and license are available at
https://github.com/cortes-ciriano-lab/SComatic. Keep the SComatic checkout,
the installed environment, and the reference resource versions recorded in
your run documentation.

### SComatic version used for the published v3 defaults

The v3 analyses were run against upstream SComatic commit
`5fd34e6d6b52524a0e392eeee6cb2f86d497863e` (`5fd34e6` in abbreviated form).
That checkout was detached at this commit and had no tracked source edits. Its
only local extra file was a coordinate-compatible BED reference resource, which
is an input to this workflow rather than a modification of SComatic code.

## Input files

Keep manifests, metadata, and lightweight SComatic resources in permanent
`/storage3` storage. Large source BAMs and managed reference FASTAs/GTFs may
remain in their existing locations. `docs/resources.md` gives a recommended
layout; the schemas in `assets/` define the required column order.

### Library manifest: `library_manifest.tsv`

The manifest has one row per source library and uses the exact header in
`assets/library_manifest.schema.tsv`:

| Column | Requirement and use |
| --- | --- |
| `donor_id` | Original donor identifier. It is retained for provenance and must agree across that donor's two rows. |
| `donor_safe` | Stable filesystem-safe donor name used in result paths and report names. It must agree across both rows and be unique across donors. |
| `reference_key` | Key that joins the donor to one row of `reference_config.tsv`. Both libraries for a donor must use the same key. |
| `tissue` | Exactly `CSF` or `PBMC`. Each donor must have exactly one of each. |
| `source_id` | Unique source-library identifier. It prefixes barcodes during splitting and distinguishes CSF/PBMC libraries with overlapping raw barcodes. |
| `meta_prefix` | Prefix before the first underscore in metadata barcodes for this library, for example `GSM123` in `GSM123_AAAC...`. |
| `bam` | Absolute path to a coordinate-sorted, indexed BAM. The workflow requires both this file and `<bam>.bai`. The BAM must carry 10x cell barcodes in the `CB` tag; `nM` and `NH` tags are required by this workflow's Step 1 read filters. |

Do not include multiple replicates for one tissue in this initial paired-input
implementation. Merge/select them before creating the manifest, or extend the
workflow deliberately with matching provenance checks.

### Cell metadata: `cell_metadata.tsv[.gz]`

The workflow reads a plain or gzip-compressed TSV. Its first three columns, in
this order, are:

| Position | Meaning |
| --- | --- |
| 1 | Barcode, prefixed as `<meta_prefix>_<10x barcode>`. |
| 2 | Tissue label, exactly `CSF` or `PBMC`. |
| 3 | v3 cell-type label. |

Only manifest-matched tissue rows are used. The workflow recognizes both
prefixed and unprefixed barcode forms when passing metadata to SComatic, and
removes raw barcode collisions between the two source libraries. It excludes
only the exact label `other`; `other_myeloid` is retained. Exclude doublets and
unassigned cells upstream rather than relying on this workflow to identify
them. The default tissue-aware labels are
`B, CSF_DC, CSF_mac, CSF_mono, NK, PBMC_DC, PBMC_mac, PBMC_mono, Plasma, T,
other_myeloid, progenitors`.

### Reference configuration: `reference_config.tsv`

The header and column order are defined by `assets/reference_config.schema.tsv`.
Every `reference_key` used in the manifest must match exactly one row.

| Column | Requirement |
| --- | --- |
| `ref` | Reference FASTA; the matching `<ref>.fai` index must exist. |
| `gtf` | GTF used for nearest-gene annotation. |
| `high_quality_bed` | Coordinate-compatible high-quality/callable-region BED supplied to SComatic Step 4. |
| `rna_editing` | RNA-editing site table supplied to SComatic Step 4.2. |
| `pon` | SComatic panel of normals supplied to Step 4.2. |

All resources must use the same genome build and contig naming convention as
the BAMs. `bin/provision_compute1_resources.sh` can copy lightweight metadata,
BED, RNA-editing, and PoN files into durable storage; it does not download
references or infer compatible resources.

### Parameters JSON

Start from `assets/compute1.params.example.json`. The required path fields are
`manifest`, `cell_metadata`, `reference_config`, `scomatic_env`,
`scomatic_repo`, and `run_root`. The first five must already exist; use a new,
writable Scratch1 directory for `run_root` because it stores large intermediate
BAMs, logs, markers, and Nextflow work. Keep the params JSON outside the Git
checkout if it contains project-specific locations.

The published scientific defaults are `max_cell_types = 15`,
`min_base_quality = 30`, and `gnomad_af_cutoff = 0.001`. The latter removes an
allele when gnomAD's conservative exome/genome `AF_grpmax` reaches 0.1%; it
uses generic cohort AF only when grpmax is unavailable.

## How to run on Compute1

1. Clone this repository in durable project storage and complete the
   installation above.
2. Prepare and audit the three input files. Confirm every manifest BAM and
   `<bam>.bai` exists, metadata contains cells for each library/tissue, and
   every `reference_key` has exactly one compatible reference-config row.
3. Copy `assets/compute1.params.example.json` to a durable project-specific
   location and fill in the absolute paths:

```bash
cp assets/compute1.params.example.json \
  /storage3/path/to/params/my_v3_run.json
```

4. If your Compute1 Docker integration requires explicit mounts, export
   `LSF_DOCKER_VOLUMES` so the repository, inputs, reference paths, SComatic
   checkout/environment, and selected Scratch1 `run_root` are visible inside
   every task container. The queue defaults to `general`; override `lsf_queue`
   or `lsf_group` only when your account requires it. The submission helper
   applies these fields to both the controller and its child tasks.
5. Validate configuration without submitting scientific tasks:

```bash
bash -n bin/*.sh
nextflow run main.nf -profile compute1 \
  -params-file /storage3/path/to/params/my_v3_run.json -preview
```

6. Submit the small Nextflow controller, which in turn submits LSF tasks:

```bash
bin/submit_compute1.sh --params /storage3/path/to/params/my_v3_run.json
```

7. Monitor the controller log in `<run_root>/logs/` and LSF child jobs. A task
   is complete only after its LSF status, combined log, success marker, and
   expected outputs agree. `docs/compute1.md` describes the log layout and
   restart behavior.

For an interrupted run with unchanged inputs, parameters, and `run_root`, use:

```bash
bin/submit_compute1.sh --params /storage3/path/to/params/my_v3_run.json --resume
```

Do not use `--resume` after changing the manifest, metadata, reference
configuration, scientific parameters, or SComatic version. Use a new run root
instead.

## Compute1 example: GSE133028_6

`examples/GSE133028_6/` is a two-library, paired CSF/PBMC example using public
GSE133028 BAM paths already available on
Compute1. It includes the exact library manifest and a params template that
references the durable `storage3` metadata and reference configuration. No BAM,
metadata, or result data is included in this repository. See the
`examples/GSE133028_6/README.md` run guide before launching it.

## Outputs

Each donor is written under `<run_root>/<donor_safe>/results/`. The terminal
gnomAD directory includes the full allele cache, one complete annotated PASS
table per donor, one `gnomad_filtered` table per donor, and a retention summary.
The exact `gnomAD_common_variant_filter` value remains in both tables so that
downstream users can audit each retention/removal decision.

## Pipeline-added annotation columns

The following columns are added after standard SComatic Step 4. They are
present in the `*.pass.ref_support.annotated.tsv` table; the gnomAD fields are
then appended to the `*.gnomad_popmax.tsv` and `*.gnomad_filtered.tsv` tables.
`NA` denotes an unavailable value, not a zero measurement.

### Nearest gene

| Column | Meaning |
| --- | --- |
| `NearestGene` | Gene identifier for the closest GTF `gene` feature. |
| `NearestGeneName` | Corresponding gene symbol/name from the GTF. |
| `NearestGeneDistance` | Distance in base pairs to that gene body; `0` means the variant overlaps the gene body. Ties use the first genomic feature deterministically. |

### Cell-type reference-support summaries

These summaries are reconstructed from the per-cell-type SComatic INFO fields.
A cell type is considered sufficiently covered for these calls when it has at
least 5 total bases and at least 5 covered cells at the locus.

| Columns | Meaning |
| --- | --- |
| `N_REF_Cell_types`, `REF_Cell_types` | Number and comma-delimited names of cell types with reference support: at least 3 reference-base observations across at least 2 cells, in addition to sufficient coverage. |
| `N_REF_VAF_1_Cell_types`, `REF_VAF_1_Cell_types` | Subset of reference-supporting cell types with no observed non-reference base among the counted high-quality bases—i.e., reference VAF of 1 at this site under the SComatic count representation. |
| `N_no_ALT_support_Cell_types`, `no_ALT_support_Cell_types` | Number and comma-delimited names of sufficiently covered cell types that do **not** meet the ALT-support rule (at least 3 ALT-base observations across at least 2 cells). This is broader than `REF_VAF_1`: it may include cell types with low-level ALT evidence. |

### gnomAD v4.1.1 annotation and common-variant decision

The pipeline queries the public indexed gnomAD v4.1.1 exome and genome VCFs
with bcftools. Annotation is allele-specific: repeated SComatic ALT entries
arising from multiple supporting cell types are collapsed before lookup. At
multi-allelic sites, one comma-delimited value is emitted per distinct ALT in
first-observed order.

| Columns | Meaning |
| --- | --- |
| `gnomAD_v4_1_1_variant_id` | Normalized `chrom-pos-ref-alt` identifier constructed for a matched gnomAD allele. |
| `gnomAD_v4_1_1_rsid` | dbSNP rsID reported in the gnomAD VCF `ID` field, if present. |
| `gnomAD_v4_1_1_exomes_AF_grpmax`, `gnomAD_v4_1_1_exomes_grpmax_ancestry`, `gnomAD_v4_1_1_exomes_AF` | gnomAD exome group-maximum AF, the ancestry group achieving it, and generic cohort-wide AF. |
| `gnomAD_v4_1_1_genomes_AF_grpmax`, `gnomAD_v4_1_1_genomes_grpmax_ancestry`, `gnomAD_v4_1_1_genomes_AF` | Equivalent values from gnomAD genomes. |
| `gnomAD_v4_1_1_popmax_AF`, `gnomAD_v4_1_1_popmax_ancestry` | Conservative maximum `AF_grpmax` across the exome and genome records and the associated ancestry group. This is not a combined-cohort AF. |
| `gnomAD_v4_1_1_generic_AF`, `gnomAD_v4_1_1_generic_AF_source` | Maximum generic `AF` across exomes/genomes and the source dataset. It is used only when `popmax_AF` is unavailable. |
| `gnomAD_v4_1_1_match_status` | Whether the exact allele matched neither, exomes only, genomes only, or both gnomAD datasets. |
| `gnomAD_common_variant_filter` | Audit-ready retention/removal decision. `FILTER_POPMAX_AF_GE_0.001` removes alleles with `popmax_AF >= 0.1%`; `FILTER_GENERIC_AF_GE_0.001` applies the generic-AF fallback only if grpmax is absent. `RETAIN_*` values identify the corresponding retained condition. |
| `gnomAD_v4_1_1_exomes_INFO`, `gnomAD_v4_1_1_genomes_INFO` | Complete original VCF INFO strings retained in the allele cache/output for future annotation without repeating the remote lookup. |

The default `gnomad_filtered` table removes only `FILTER_*` variants. Alleles
not present in gnomAD, or matched alleles without an available AF/grpmax value,
remain and are explicitly labelled rather than treated as frequency zero.

## License and citation

This workflow wrapper and its original documentation are distributed under the
MIT License; see `LICENSE`. The repository invokes SComatic as an external
dependency and does not relicense SComatic, gnomAD, or any referenced input
resources. Users must comply with each upstream resource's terms, including
SComatic's academic-use licensing conditions.

Please cite this workflow using `CITATION.cff`, as well as the SComatic and
gnomAD publications/resources used in an analysis.
