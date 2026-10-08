# GSE133028_6 paired CSF/PBMC example on Compute1

This example runs one paired donor from public dataset GSE133028 through the
published v3 workflow. It uses these indexed source BAMs on Compute1:

| Tissue | Source library | BAM |
| --- | --- | --- |
| CSF | `GSM5264639` | `/storage3/fs1/gfwu/Active/public_datasets/GSE133028/full_output/GSM5264639/outs/possorted_genome_bam.bam` |
| PBMC | `GSM5264640` | `/storage3/fs1/gfwu/Active/public_datasets/GSE133028/full_output/GSM5264640/outs/possorted_genome_bam.bam` |

`library_manifest.tsv` records those paths, the verified barcode prefixes, and
the `GRCh38_2020_A` reference key. The metadata and reference-config paths in
`params.json` point to the durable private resources provisioned for this
Compute1 deployment. They are referenced only; neither resource is committed
to this repository.

## Before running

1. Place this repository at the durable location encoded in `params.json`, or
   update `manifest` to the absolute path of your clone's
   `examples/GSE133028_6/library_manifest.tsv`.
2. Set `scomatic_env` to a validated SComatic Python environment and
   `scomatic_repo` to the corresponding SComatic checkout. Both paths must be
   readable from LSF Docker tasks.
3. Choose a new, writable scratch `run_root`. Do not point at a previous
   production GSE133028_6 run.
4. Confirm that both `<bam>.bai` files exist and that the `GRCh38_2020_A` row
   is present in the referenced `reference_config.tsv`.
5. Export account-appropriate Docker mounts if required by Compute1. For this
   setup, the repository, the `/storage3` inputs/resources, the selected
   `/scratch1` run root, and the SComatic environment must be visible inside
   the task container.

## Run

Copy the template outside the Git checkout so local edits are never committed:

```bash
cp examples/GSE133028_6/params.json \
  /storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF/Analysis/codex/testing/scomatic-v3-nextflow-config/GSE133028_6.params.json
```

Edit the copied JSON as described above, then submit the controller from the
repository root:

```bash
bin/submit_compute1.sh \
  --params /storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF/Analysis/codex/testing/scomatic-v3-nextflow-config/GSE133028_6.params.json
```

For an interrupted run with unchanged paths and parameters, use:

```bash
bin/submit_compute1.sh \
  --params /storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF/Analysis/codex/testing/scomatic-v3-nextflow-config/GSE133028_6.params.json \
  --resume
```

The default configuration retains the v3 settings: separate tissue-aware
myeloid labels, `max_cell_types = 15`, base quality 30, and gnomAD removal at
`AF_grpmax >= 0.001` (with generic-AF fallback only when grpmax is absent).

## Expected outputs

The workflow creates `GSE133028_6/results/` below the selected `run_root`.
Completion requires the Step 1, joint-BAM, Step 2, Step 3, Step 4, and gnomAD
success markers, plus the terminal `gnomad_filtered` table and retention
summary. See `README.md`, "Outputs," for the expected table classes and
annotation columns.
