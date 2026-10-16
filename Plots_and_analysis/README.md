# Plots and analysis

`Plots_and_analysis/csf_analyze_SComatic_codex.R` plots completed, gnomAD-filtered
v3 results. It writes four derived variant tables and five PDFs per donor:
cell-type-specific ALT/REF support, putative clonal ALT/REF support, and three
clonal heatmaps with the additional annotations (support, ALT VAF, and ALT+
cell fraction). The quantitative heatmaps follow the annotated support plot's
row and column order.

## How to run

Use an R environment containing `dplyr`, `tidyr`, `stringr`, `purrr`, `tibble`,
`data.table`, `RColorBrewer`, `ComplexHeatmap`, `fastcluster`, `circlize`,
`R.utils`, and `png`. `grid` is included with R. Packages must already be
installed. The script does not require `every.script.v1.R`.

From the repository root, in that R environment:

```bash
SCOMATIC_DONOR_ID='GSE133028:6' \
SCOMATIC_OUTPUT_ROOT=/absolute/path/to/new-plot-results \
Rscript --vanilla Plots_and_analysis/csf_analyze_SComatic_codex.R
```

On Compute1, run this command inside an LSF job with the input directories
mounted. Testing used R 4.5.3 and ComplexHeatmap 2.26.1 in:

```text
/storage1/fs1/martyomov/Active/snakemake_caches/i.strunilin/.conda/envs/550_base1
```

The default donor is `GSE133028:6`. Default input paths are encoded near the
top of the script and point to existing Storage3 files. The plotting project
root is:

```text
/storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF/Analysis/codex/testing/scomatic_v3_plotting
```

Under that root, `config/v3_plot_inputs.tsv` maps original donor IDs to their
preserved `gnomad_filtered` tables. `inputs/` contains the original v0.3
clinical table and C4's three Step 1 metadata files for its replicate-aware
cell counts. The grouped cell metadata is the v0.5 table in
`Analysis/codex/testing/scomatic-v3-nextflow/resources/metadata/` under the
same CSF project. These project input files remain on Compute1.

Override paths with `SCOMATIC_INPUT_MANIFEST`, `SCOMATIC_CLINICAL`,
`SCOMATIC_CELL_METADATA`, or `SCOMATIC_C4_INPUT_ROOT`. An alternative manifest
must have `donor_id` and `variant_table` columns. Donor output folders use the
original ID with `:` replaced by `_`.

## Current generated plots

The tested outputs for all 13 completed donors are on permanent Storage3:

```text
/storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF/Analysis/codex/testing/scomatic_v3_plotting/runs/20261016_r1/results/
```

For example, GSE133028:6 is in the `GSE133028_6/` subdirectory. Local copies
of these remote outputs are in:

```text
/Users/istr/Library/CloudStorage/Box-Box/Artyomov_Lab/Projects/CSF/codex/csf_install_SComatic/results/v3_plot_cleanup_20261016/remote/results/
```

The run produced 65 PDFs, and all 52 derived tables matched the existing
results. Logs, donor success markers, and `validation_summary.tsv` are retained
with the remote run. A representative plot check confirmed the layout and
row order; the remote R environment orders two pure-REF legend entries
differently. Choose a new output root when running another analysis.
