# Inputs

## Library manifest

Use one row for CSF and one row for PBMC per donor. The exact column order is
shown in `assets/library_manifest.schema.tsv`.

- `donor_id`: original source donor identifier, retained for provenance.
- `donor_safe`: stable, filesystem-safe name used in results and reports.
- `reference_key`: joins to `reference_config.tsv`.
- `tissue`: exactly `CSF` or `PBMC`.
- `source_id`: source library name used to prefix cell barcodes.
- `meta_prefix`: prefix in the first metadata column before `_`.
- `bam`: absolute path to an indexed source BAM. `<bam>.bai` must exist.

## Cell metadata

This v3 implementation expects a tabular file (plain or gzip-compressed) whose
first three columns are, in order: barcode, tissue, and v3 cell type. Barcodes
must begin with `<meta_prefix>_`; both prefixed and unprefixed forms are made
available to SComatic to handle source BAM barcode conventions. The workflow
uses only rows matching the manifest tissue, excludes exactly `other`, and
preserves `other_myeloid`.

Do not commit clinical or cell metadata tables. Store them under durable
project storage and record their paths in your params JSON.

## Reference configuration

`assets/reference_config.schema.tsv` maps each `reference_key` to absolute
paths for FASTA (with `.fai`), GTF, high-quality BED, RNA-editing table, and
SComatic panel of normals. These resources must use coordinates compatible with
the source BAMs.
