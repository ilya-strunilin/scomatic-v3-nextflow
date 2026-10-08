# Persistent resource layout

Keep durable, project-specific inputs outside Git. A suggested Compute1 layout
is:

```text
/storage3/fs1/gfwu/Active/<user>/Projects/<project>/scomatic-v3-resources/
  metadata/cell_metadata_v3.tsv.gz
  manifests/library_manifest_<cohort>.tsv
  references/reference_config.tsv
  references/<small SComatic resource files>
  params/<run>.json
```

Large reference FASTAs/GTFs and source BAMs can remain at their existing
managed locations; their absolute paths belong in `reference_config.tsv` and
the library manifest. Copy lightweight metadata, manifests, panels of normals,
editing tables, and BED resources into this durable tree with checksums before
starting a cohort. Keep `run_root` in scratch because it contains large
intermediate BAMs and Nextflow work.

This split makes the GitHub repository portable while keeping project data and
HPC paths private.
