### csf_analyze_SComatic_codex
###
### desc: Analyze preserved SComatic v3 results and plot donor heatmaps.

options(digits = 3, scipen = 3)
set.seed(1)

## Load packages

my0packages <- c(
  'dplyr',
  'tidyr',
  'stringr',
  'purrr',
  'tibble',
  'data.table',
  'RColorBrewer',
  'ComplexHeatmap',
  'fastcluster',
  'grid'
)

for (pkg_name_tmp in my0packages) {
  library(package = pkg_name_tmp, character.only = T, warn.conflicts = FALSE)
}


### Functions -----

dfread <- function(l.path, ...) {
  data.table::fread(file = l.path, data.table = FALSE, ...)
}

nafwrite <- function(l.table, l.path, ...) {
  data.table::fwrite(x = l.table, file = l.path, na = 'NA',
                    quote = FALSE, sep = '\t', ...)
}

### -----


# Set SCOMATIC_DONOR_ID for a different donor; the default is the public example.
donor_id = Sys.getenv('SCOMATIC_DONOR_ID', unset = 'GSE133028:6')
donor_folder_id = str_replace_all(donor_id, ':', '_')

# All scientific inputs are read from permanent Storage3.
project_root = '/storage3/fs1/gfwu/Active/i.strunilin/Projects/CSF'
plot_root = file.path(project_root, 'Analysis/codex/testing/scomatic_v3_plotting')
input_manifest = Sys.getenv(
  'SCOMATIC_INPUT_MANIFEST',
  unset = file.path(plot_root, 'config/v3_plot_inputs.tsv')
)
clinical_path = Sys.getenv(
  'SCOMATIC_CLINICAL',
  unset = file.path(plot_root, 'inputs/CSF_atlas_clinical_v0.3_2026-08-12.tsv')
)
cell_metadata_path = Sys.getenv(
  'SCOMATIC_CELL_METADATA',
  unset = file.path(
    project_root, 'Analysis/codex/testing/scomatic-v3-nextflow/resources/metadata',
    'cell_meta_grouped_csf_SComatic_prepare_anno2_clin_0.5.tsv.gz'
  )
)
c4_input_root = Sys.getenv(
  'SCOMATIC_C4_INPUT_ROOT', unset = file.path(plot_root, 'inputs/C4')
)
output_root = Sys.getenv(
  'SCOMATIC_OUTPUT_ROOT', unset = file.path(plot_root, 'results')
)
output_dir <- file.path(output_root, donor_folder_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)


### Load files ----
clinical = dfread(clinical_path)
cell_meta2 = dfread(cell_metadata_path)
### ----


#### Analyze joint CSF/PBMC SComatic v3 results.
####
#### v3 preserves source-specific myeloid labels (CSF_mac, CSF_mono, CSF_DC,
#### PBMC_mac, PBMC_mono, PBMC_DC) and uses Step 4.1 --max_cell_types 15.
#### Inputs are preserved gnomAD-filtered final Step 4.2 PASS tables with
#### reference-support and nearest-gene annotations. ----

# v3 input: final Step 4.2 PASS table with reference-support and gene
# annotations. v3 uses --max_cell_types 15 and does not add donor/joint
# prefixes to cell-type labels.
v3_donor_input = dfread(input_manifest) %>%
  filter(.data$donor_id == .env$donor_id)
CP_joint_bam_results_path = v3_donor_input$variant_table

# Source-specific labels are intentional: all CSF and PBMC myeloid cells
# remain separate after the joint BAM analysis.
v3_cell_type_order = c(
  'CSF_mac', 'CSF_mono', 'CSF_DC',
  'PBMC_mac', 'PBMC_mono', 'PBMC_DC',
  'B', 'T', 'NK', 'Plasma', 'progenitors', 'other_myeloid', 'other'
)
# Derive each donor's CSF/PBMC library IDs from the clinical table.  The
# barcode prefixes in cell_meta2 omit the study prefix (for example,
# GSE133028:GSM5264639 becomes GSM5264639).
v3_sample_ids = clinical %>%
  filter(donor == donor_id, tissue %in% c('CSF', 'PBMC')) %>%
  arrange(match(tissue, c('CSF', 'PBMC')), sample) %>%
  pull(sample) %>%
  str_remove('^[^:]+:')

# Cell counts are taken from the v3 source-specific metadata used to build
# the joint analysis. They are displayed as a column annotation below Cell.
v3_cell_counts = if (donor_id == 'GSE243905:C4') {
  # C4 has one CSF and two PBMC libraries.  Use the exact Step 1 metadata
  # from this v3 run, after inter-library barcode-overlap removal, and sum
  # the two PBMC replicate counts before constructing the joint annotation.
  c4_input = tibble(
    path = file.path(
      c4_input_root,
      c('C4_CSF.cell_types.tsv',
        'C4_PBMC_REP1.cell_types.tsv',
        'C4_PBMC_REP2.cell_types.tsv')
    ),
    tissue = c('CSF', 'PBMC', 'PBMC')
  )

  map2_dfr(c4_input$path, c4_input$tissue, function(path, tissue) {
    dfread(path) %>%
      count(Cell_type, name = 'n_cells') %>%
      transmute(
        alt_major = case_when(
          tissue == 'CSF' & Cell_type %in% c('mac', 'mono', 'DC') ~ paste0('CSF_', Cell_type),
          tissue == 'PBMC' & Cell_type %in% c('mac', 'mono', 'DC') ~ paste0('PBMC_', Cell_type),
          TRUE ~ Cell_type
        ),
        n_cells
      )
  }) %>%
    group_by(alt_major) %>%
    summarise(n_cells = sum(n_cells), .groups = 'drop') %>%
    right_join(tibble(alt_major = v3_cell_type_order), by = 'alt_major') %>%
    mutate(n_cells = coalesce(n_cells, 0L)) %>%
    arrange(match(alt_major, v3_cell_type_order))
} else {
  cell_meta2 %>%
    filter(str_detect(V1, paste0('^(', paste(v3_sample_ids, collapse = '|'), ')_'))) %>%
    count(alt_major, name = 'n_cells') %>%
    right_join(tibble(alt_major = v3_cell_type_order), by = 'alt_major') %>%
    mutate(n_cells = coalesce(n_cells, 0L)) %>%
    arrange(match(alt_major, v3_cell_type_order))
}

stopifnot(identical(v3_cell_counts$alt_major, v3_cell_type_order))

v3_cell_type_annotation = tibble(
  cell_type = v3_cell_type_order,
  Tissue = case_when(
    str_starts(cell_type, 'CSF_') ~ 'CSF',
    str_starts(cell_type, 'PBMC_') ~ 'PBMC',
    TRUE ~ 'Both'
  ),
  Cell = str_remove(cell_type, '^(CSF|PBMC)_')
) %>%
  mutate(
    Tissue = factor(Tissue, levels = c('CSF', 'PBMC', 'Both')),
    Cell = factor(
      Cell,
      levels = c('mac', 'mono', 'DC', 'B', 'T', 'NK', 'Plasma',
                 'progenitors', 'other_myeloid', 'other')
    )
  )

# Qualitative ColorBrewer Set2 colors, used consistently across v3 plots.
v3_tissue_colors = c(CSF = '#66C2A5', PBMC = '#FC8D62', Both = '#8DA0CB')
v3_broad_type_colors = c(
  mac = '#8DA0CB', mono = '#E78AC3', DC = '#FC8D62',
  B = '#66C2A5', T = '#7570B3', NK = '#A6D854', Plasma = '#FFD92F',
  progenitors = '#E66101', other_myeloid = '#B3B3B3', other = '#BDBDBD'
)

# Count comma-delimited ALT-supporting cell types. Cell_types has no
# donor/joint prefix in v3, so these values are directly comparable across
# the v3 donors.
v3_count_cell_types = function(x) {
  if_else(is.na(x) | x == '', 0L, str_count(x, ',') + 1L)
}

# Match one exact comma-delimited cell-type label.  Wrapping with commas
# prevents `other` from matching labels such as `other_myeloid`.
v3_has_cell_type = function(x, cell_type) {
  str_detect(paste0(',', coalesce(x, ''), ','), paste0(',', cell_type, ','))
}

# A strict single-label check for REF-VAF-1 support.
v3_only_cell_type = function(x, cell_type) {
  !is.na(x) & x == cell_type
}

CP_joint_bam_full_results = dfread(CP_joint_bam_results_path) %>%
  mutate(
    ALT = str_sub(ALT, 1, 1),
    variant_id = paste0(CHROM, ':', Start, ':', REF, '>', ALT),
    safe_variant_id = paste(CHROM, Start, End, REF, ALT, sep = '_'),
    N_ALT_Cell_types = v3_count_cell_types(Cell_types)
  ) %>%
  relocate(variant_id, safe_variant_id, .before = 1) %>% 
  # remove variants with pure REF support only in other or other_myeloid
  filter(!(REF_VAF_1_Cell_types %in% c('other', 'other_myeloid')))

message(paste(nrow(CP_joint_bam_full_results), 'v3 final PASS variants'))
nafwrite(CP_joint_bam_full_results, 'v3_full_results.tsv')

CP_joint_bam_variants_with_ref_celltype_df = CP_joint_bam_full_results %>%
  filter(N_REF_VAF_1_Cell_types > 0)
CP_joint_bam_specific_variants_df = CP_joint_bam_variants_with_ref_celltype_df %>%
  filter(N_ALT_Cell_types == 1)
CP_joint_bam_clonal_variants_df = CP_joint_bam_variants_with_ref_celltype_df %>%
  filter(N_ALT_Cell_types > 1) %>%
  # Do not classify a site as clonal when its only second ALT+ label is the
  # catch-all `other` compartment.
  filter(!(N_ALT_Cell_types == 2L & v3_has_cell_type(Cell_types, 'other'))) %>%
  # Keep the clonal rule explicit even though this condition is also
  # excluded above for the complete v3 result set.
  filter(!(N_REF_VAF_1_Cell_types == 1L &
           v3_only_cell_type(REF_VAF_1_Cell_types, 'other')))

message(paste(nrow(CP_joint_bam_variants_with_ref_celltype_df), 'variants with N_REF_VAF_1_Cell_types > 0'))

message(paste(nrow(CP_joint_bam_specific_variants_df), 'variants with N_ALT_Cell_types == 1, N_REF_VAF_1_Cell_types > 0'))

message(paste(nrow(CP_joint_bam_clonal_variants_df), 'variants with N_ALT_Cell_types > 1, N_REF_VAF_1_Cell_types > 0\n'))


CP_joint_bam_variants_with_multi_alt_celltype_df = CP_joint_bam_full_results %>%
  filter(N_ALT_Cell_types > 1)

message(paste(nrow(CP_joint_bam_variants_with_multi_alt_celltype_df), 'variants with N_ALT_Cell_types > 1\n'))


nafwrite(
  CP_joint_bam_variants_with_ref_celltype_df,
  'v3_variants_with_REF_VAF_1_in_gt1_cell_type.tsv'
)
nafwrite(CP_joint_bam_specific_variants_df, 'v3_cell_type_specific_variants.tsv')
nafwrite(CP_joint_bam_clonal_variants_df, 'v3_putative_clonal_variants.tsv')

# Convert ALT and strict REF support to the three-state heatmap:
# black = ALT, grey = REF VAF 1, white = no coverage. ALT takes
# precedence if a label appears in both source columns.
v3_make_alt_ref_matrix = function(results_df) {
  v3_support_df = results_df %>%
    select(variant_id, Cell_types, REF_VAF_1_Cell_types) %>%
    pivot_longer(
      cols = c(Cell_types, REF_VAF_1_Cell_types),
      names_to = 'support_type', values_to = 'cell_type'
    ) %>%
    filter(!is.na(cell_type), cell_type != '') %>%
    separate_rows(cell_type, sep = ',') %>%
    filter(cell_type %in% v3_cell_type_order) %>%
    mutate(value = if_else(support_type == 'Cell_types', 1L, -1L)) %>%
    group_by(variant_id, cell_type) %>%
    summarise(value = max(value), .groups = 'drop')

  v3_mat = matrix(
    0L,
    nrow = nrow(results_df),
    ncol = length(v3_cell_type_order),
    dimnames = list(results_df$variant_id, v3_cell_type_order)
  )
  if (nrow(v3_support_df) > 0) {
    v3_mat[cbind(
      match(v3_support_df$variant_id, rownames(v3_mat)),
      match(v3_support_df$cell_type, colnames(v3_mat))
    )] = v3_support_df$value
  }
  v3_mat
}

# Fixed body geometry keeps the three PDFs aligned for a donor even when
# the legend text or dendrogram widths differ slightly.
v3_clonal_heatmap_body_width = unit(90, 'mm')
v3_clonal_heatmap_row_height = unit(4, 'mm')
v3_clonal_target_row_order = NULL
v3_clonal_target_row_dendrogram = NULL
v3_clonal_support_hclust = function(results_df) {
  if (nrow(results_df) <= 1L) return(FALSE)
  hclust(dist(v3_make_alt_ref_matrix(results_df)), method = 'ward.D2')
}

v3_save_support_heatmap = function(results_df, filename, width, height,
                                   show_myeloid_alt_annotation = TRUE) {
  # Always replace an existing output.  In particular, a rerun that leaves
  # no eligible clonal sites must not retain a stale heatmap from an earlier
  # analysis.
  if (nrow(results_df) == 0) {
    pdf(filename, width = width, height = height, useDingbats = FALSE,
        bg = 'transparent')
    grid.newpage()
    grid.text(
      'No eligible variants after clonal-site filtering',
      gp = gpar(fontsize = 18)
    )
    dev.off()
    return(invisible(NULL))
  }
  
  heatmap_font_size = 12
  
  ht_opt(
    simple_anno_size = unit(heatmap_font_size+4, 'pt'),
    fast_hclust = TRUE,
    legend_title_gp = gpar(fontsize = heatmap_font_size, fontface = 'bold'),
    legend_labels_gp = gpar(fontsize = heatmap_font_size),
    legend_grid_width = unit(5, 'mm'),
    legend_grid_height = unit(5, 'mm'),
    heatmap_row_names_gp = gpar(fontsize = heatmap_font_size-4),
    heatmap_column_names_gp = gpar(fontsize = heatmap_font_size+4),
    HEATMAP_LEGEND_PADDING =  unit(0.2, 'cm')
  )
  

  v3_mat = v3_make_alt_ref_matrix(results_df)
  v3_n_myeloid_alt = results_df$Cell_types %>%
    str_split(',') %>%
    map_int(function(cell_types) {
      cell_types %>%
        str_remove('^(CSF|PBMC)_') %>%
        intersect(c('mac', 'mono', 'DC')) %>%
        length()
    })
  v3_myeloid_alt_annotation = case_when(
    v3_n_myeloid_alt >= 3 ~ '>=3 myeloid types',
    v3_n_myeloid_alt == 2 ~ '2 myeloid types',
    TRUE ~ '<2 myeloid types'
  )
  
  v3_n_ref_vaf1 = case_when(results_df$N_REF_VAF_1_Cell_types >= 2 ~ '>=2',
                            results_df$N_REF_VAF_1_Cell_types == 1 ~ '1',
                            TRUE ~ '0'
                            )
  
  v3_myeloid_alt_ha = rowAnnotation(
    `ALT in mac/mono/DC` = v3_myeloid_alt_annotation,
    `pure REF support` = v3_n_ref_vaf1,
    col = list(
      `ALT in mac/mono/DC` = c(
        '<2 myeloid types' = 'white',
        '2 myeloid types' = '#A1D99B',
        '>=3 myeloid types' = '#006D2C'
      ),
      `pure REF support` = c('0' = 'white', '1' = '#fee391','>=2' = '#fe9929')
    ),
    simple_anno_size = unit(4, 'mm'),
    border = TRUE
  )
  
  v3_cell_count_breaks = pretty(c(0, max(v3_cell_counts$n_cells)), n = 3)
  v3_cell_count_breaks = v3_cell_count_breaks[
    v3_cell_count_breaks >= 0 & v3_cell_count_breaks <= max(v3_cell_counts$n_cells)
  ]
  v3_column_df = v3_cell_type_annotation %>%
    column_to_rownames('cell_type') %>%
    as.data.frame()
  
  v3_column_ha = HeatmapAnnotation(
    Tissue = v3_column_df$Tissue,
    Cell = v3_column_df$Cell,
    col = list(Tissue = v3_tissue_colors, Cell = v3_broad_type_colors),
    `sqrt(N cells)` = anno_barplot(
      sqrt(v3_cell_counts$n_cells),
      baseline = 0,
      border = FALSE,
      gp = gpar(fill = '#5B7FA3', col = NA),
      axis = TRUE,
      axis_param = list(
        side = 'left',
        at = sqrt(v3_cell_count_breaks),
        labels = v3_cell_count_breaks,
        labels_rot = 0
      ),
      numbers_gp = gpar(fontsize = heatmap_font_size-2),
      height = unit(14, 'mm')
    ),
    annotation_name_gp = gpar(fontsize = heatmap_font_size+4, fontface = 'bold'),
    gap = unit(c(1, 1), 'mm')
  )
  v3_column_spacer_ha = HeatmapAnnotation(
    spacer = anno_empty(height = unit(1.5, 'mm'), border = FALSE),
    show_annotation_name = FALSE
  )

  v3_ht = Heatmap(
    v3_mat,
    name = 'Support',
    col = c('-1' = 'grey70', '0' = 'white', '1' = 'black'),
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      labels = c('pure REF', 'No coverage', 'ALT+')
    ),
    top_annotation = c(v3_column_ha, v3_column_spacer_ha),
    left_annotation = if (show_myeloid_alt_annotation) v3_myeloid_alt_ha else NULL,
    cluster_rows = v3_clonal_support_hclust(results_df),
    row_dend_reorder = FALSE,
    cluster_columns = FALSE,
    show_row_names = nrow(results_df) <= 70,
    row_names_side = 'left',
    row_dend_side = 'right',
    row_dend_width = unit(3, 'mm'),
    row_names_max_width = unit(5, 'cm'),
    use_raster = TRUE,
    raster_quality = 3,
    rect_gp = gpar(col = 'grey85')
  )
  
  pdf(filename, width = width, height = height, useDingbats = FALSE, bg = 'transparent')
  v3_drawn_ht = ComplexHeatmap::draw(
    v3_ht,
    heatmap_legend_side = 'right',
    annotation_legend_side = 'right',
    merge_legends = TRUE,
    padding = unit(c(5, 15, 5, 5), 'mm'),
    background = 'transparent'
  )
  dev.off()
}

v3_save_support_heatmap(
  CP_joint_bam_clonal_variants_df,
  'v3_putative_clonal_ALT_REF_heatmap.pdf',
  width = 8.5, height = 7
)

v3_save_support_heatmap(
  CP_joint_bam_specific_variants_df,
  'v3_cell_type_specific_ALT_REF_heatmap.pdf',
  width = 9, height = 8,
  show_myeloid_alt_annotation = FALSE
)

####
#### CSF-mac co-ALT annotation heatmap
# This companion block shares variant filtering with the primary v3 block,
# but constructs its own cell annotations and plotting function.
v3_new_anno_cell_type_order = c(
  'CSF_mac', 'PBMC_mono', 'CSF_mono', 'CSF_DC',
  'PBMC_mac', 'PBMC_DC',
  'B', 'T', 'NK', 'Plasma', 'progenitors', 'other_myeloid'
)

# Infer the broad hematopoietic ancestor from ALT-supporting cell types.
# `other` is not plotted in this heatmap and is excluded here; the explicit
# `other_myeloid` label remains part of the myeloid compartment.
v3_infer_ancestor = function(alt_cell_types) {
  broad_types = alt_cell_types %>%
    coalesce('') %>%
    str_split(',') %>%
    pluck(1) %>%
    str_remove('^(CSF|PBMC)_') %>%
    discard(~ .x %in% c('', 'other')) %>%
    unique()
  lymphoid = c('B', 'T', 'NK', 'Plasma')
  myeloid = c('mac', 'mono', 'DC', 'other_myeloid')
  has_lymphoid = any(broad_types %in% lymphoid)
  has_myeloid = any(broad_types %in% myeloid)

  if (length(broad_types) == 0L) return('Unresolved')
  if ('progenitors' %in% broad_types || (has_lymphoid && has_myeloid)) return('HSPC')
  if (setequal(broad_types, c('T', 'NK'))) return('TNKP')
  if (setequal(broad_types, c('B', 'Plasma'))) return('B')
  if (length(broad_types) == 1L) return('Single')
  if (all(broad_types %in% lymphoid)) return('CLP')
  if (setequal(broad_types, c('mono', 'mac'))) return('MONO')
  if (all(broad_types %in% myeloid)) return('CMP')
  'Unresolved'
}
v3_ancestor_colors = c(
  HSPC = '#1F78B4',
  TNKP = '#D95F02',
  B = '#7570B3',
  CLP = '#E7298A',
  MONO = '#A6761D',
  CMP = '#33A02C',
  Single = '#E6AB02',
  Unresolved = '#BDBDBD'
)

# SComatic stores, for each cell type, DP and NC followed by base-indexed
# cell counts (CC) and read counts (BC). The base order is A,C,T,G,I,D.
# These helpers only parse those recorded summaries; they do not revisit BAM
# files or recompute allele evidence.
v3_scomatic_base_order = c('A', 'C', 'T', 'G', 'I', 'D')
v3_scomatic_alt_metric = function(summary_string, alt_base, metric) {
  if (is.na(summary_string) || summary_string %in% c('', 'NA', '.')) {
    return(NA_real_)
  }
  fields = str_split(summary_string, '\\|', simplify = FALSE)[[1]]
  base_index = match(alt_base, v3_scomatic_base_order)
  if (length(fields) < 4L || is.na(base_index)) return(NA_real_)

  denominator = suppressWarnings(as.numeric(fields[[if (metric == 'vaf') 1L else 2L]]))
  numerator_vector = str_split(fields[[if (metric == 'vaf') 4L else 3L]], ':', simplify = FALSE)[[1]]
  numerator = suppressWarnings(as.numeric(numerator_vector[[base_index]]))
  if (is.na(denominator) || denominator <= 0 || is.na(numerator)) return(NA_real_)
  numerator / denominator
}

v3_make_alt_metric_matrix = function(results_df, metric) {
  stopifnot(metric %in% c('vaf', 'cell_fraction'))
  alt_base = str_sub(results_df$ALT, 1, 1)
  mat = matrix(
    NA_real_, nrow = nrow(results_df), ncol = length(v3_cell_type_order),
    dimnames = list(results_df$variant_id, v3_cell_type_order)
  )
  for (cell_type in v3_cell_type_order) {
    # SComatic omits a cell-type column entirely when the donor has no
    # cells of that type. Retain the matrix column and leave it as no
    # coverage so the layout stays consistent across donors.
    if (!(cell_type %in% names(results_df))) next
    mat[, cell_type] = map2_dbl(
      results_df[[cell_type]], alt_base,
      ~ v3_scomatic_alt_metric(.x, .y, metric)
    )
  }
  mat
}

v3_new_anno_cell_counts = v3_cell_counts %>%
  filter(alt_major != 'other') %>%
  arrange(match(alt_major, v3_new_anno_cell_type_order))
v3_new_anno_cell_type_annotation = tibble(
  cell_type = v3_new_anno_cell_type_order,
  Tissue = case_when(
    str_starts(cell_type, 'CSF_') ~ 'CSF',
    str_starts(cell_type, 'PBMC_') ~ 'PBMC',
    TRUE ~ 'Both'
  ),
  Cell = str_remove(cell_type, '^(CSF|PBMC)_')
) %>%
  mutate(
    Tissue = factor(Tissue, levels = c('CSF', 'PBMC', 'Both')),
    Cell = factor(
      Cell, levels = c('mac', 'mono', 'DC', 'B', 'T', 'NK', 'Plasma',
                       'progenitors', 'other_myeloid')
    )
  )

v3_cell_type_order = v3_new_anno_cell_type_order
v3_cell_counts = v3_new_anno_cell_counts
v3_cell_type_annotation = v3_new_anno_cell_type_annotation
stopifnot(
  identical(v3_cell_counts$alt_major, v3_cell_type_order),
  identical(as.character(v3_cell_type_annotation$cell_type), v3_cell_type_order)
)

v3_new_anno_row_annotation = function(results_df) {
  v3_csf_mac_annotation = map2_chr(
    results_df$Cell_types,
    results_df$REF_VAF_1_Cell_types,
    function(alt_cell_types, ref_cell_types) {
      alt_types = str_split(coalesce(alt_cell_types, ''), ',', simplify = TRUE)
      alt_types = alt_types[alt_types != '']
      if (v3_has_cell_type(alt_cell_types, 'CSF_mac')) {
        co_alt_types = setdiff(alt_types, 'CSF_mac')
        if ('PBMC_mono' %in% co_alt_types) return('CSF mac + PBMC mono')
        if (any(co_alt_types %in% c('B', 'NK', 'T'))) return('CSF mac + B/NK/T')
        if ('PBMC_DC' %in% co_alt_types) return('CSF mac + PBMC DC')
        if (length(co_alt_types) > 0) return('CSF mac + other')
        return('CSF mac')
      }
      'CSF mac REF/no coverage'
    }
  )
  v3_n_ref_vaf1 = case_when(
    results_df$N_REF_VAF_1_Cell_types >= 2 ~ '>=2',
    results_df$N_REF_VAF_1_Cell_types == 1 ~ '1',
    TRUE ~ '0'
  )
  v3_ancestor_annotation = vapply(results_df$Cell_types, v3_infer_ancestor, character(1))

  rowAnnotation(
    Ancestor = v3_ancestor_annotation,
    `ALT+ in CSF mac` = v3_csf_mac_annotation,
    `pure REF support` = v3_n_ref_vaf1,
    col = list(
      Ancestor = v3_ancestor_colors,
      `ALT+ in CSF mac` = c(
        'CSF mac' = '#008837',
        'CSF mac + PBMC mono' = '#4DAF4A',
        'CSF mac + B/NK/T' = '#E41A1C',
        'CSF mac + PBMC DC' = '#377EB8',
        'CSF mac + other' = '#984EA3',
        'CSF mac REF/no coverage' = 'white'
      ),
      `pure REF support` = c('0' = 'white', '1' = '#fee391', '>=2' = '#fe9929')
    ),
    annotation_legend_param = list(
      Ancestor = list(ncol = 1),
      `ALT+ in CSF mac` = list(ncol = 1),
      `pure REF support` = list(ncol = 1)
    ),
    simple_anno_size = unit(4, 'mm'),
    border = TRUE,
    gap = unit(1, 'mm')
  )
}

v3_new_anno_top_annotation = function(heatmap_font_size, sparse = FALSE) {
  v3_cell_count_breaks = pretty(c(0, max(v3_cell_counts$n_cells)), n = 3)
  v3_cell_count_breaks = v3_cell_count_breaks[
    v3_cell_count_breaks >= 0 & v3_cell_count_breaks <= max(v3_cell_counts$n_cells)
  ]
  v3_column_df = v3_cell_type_annotation %>%
    column_to_rownames('cell_type') %>%
    as.data.frame()
  v3_sqrt_cell_count_anno = anno_barplot(
      sqrt(v3_cell_counts$n_cells), baseline = 0, border = FALSE,
      gp = gpar(fill = '#5B7FA3', col = NA), axis = TRUE,
      axis_param = list(
        side = 'left', at = sqrt(v3_cell_count_breaks),
        labels = v3_cell_count_breaks, labels_rot = 0
      ),
      numbers_gp = gpar(fontsize = heatmap_font_size - 2),
      height = unit(14, 'mm')
    )
  if (sparse) {
    v3_column_ha = HeatmapAnnotation(
      Tissue = v3_column_df$Tissue,
      Cell = v3_column_df$Cell,
      sqrt_title = anno_empty(height = unit(4, 'mm'), border = FALSE),
      `sqrt(N cells)` = v3_sqrt_cell_count_anno,
      col = list(Tissue = v3_tissue_colors, Cell = v3_broad_type_colors),
      annotation_legend_param = list(Tissue = list(ncol = 1), Cell = list(ncol = 1)),
      show_annotation_name = FALSE,
      gap = unit(c(1, 1, 1), 'mm')
    )
  } else {
    v3_column_ha = HeatmapAnnotation(
      Tissue = v3_column_df$Tissue,
      Cell = v3_column_df$Cell,
      `sqrt(N cells)` = v3_sqrt_cell_count_anno,
      col = list(Tissue = v3_tissue_colors, Cell = v3_broad_type_colors),
      annotation_name_gp = gpar(fontsize = heatmap_font_size + 4, fontface = 'bold'),
      annotation_name_side = 'right',
      annotation_legend_param = list(Tissue = list(ncol = 1), Cell = list(ncol = 1)),
      show_annotation_name = c(TRUE, TRUE, TRUE),
      gap = unit(c(1, 1), 'mm')
    )
  }
  v3_column_spacer_ha = HeatmapAnnotation(
    spacer = anno_empty(height = unit(1.5, 'mm'), border = FALSE),
    show_annotation_name = FALSE
  )
  c(v3_column_ha, v3_column_spacer_ha)
}


v3_save_new_anno_heatmap = function(results_df, filename, width, height,
                                   show_myeloid_alt_annotation = TRUE) {
  # Always replace an existing output.  In particular, a rerun that leaves
  # no eligible clonal sites must not retain a stale heatmap from an earlier
  # analysis.
  if (nrow(results_df) == 0) {
    pdf(filename, width = width, height = height, useDingbats = FALSE,
        bg = 'transparent')
    grid.newpage()
    grid.text(
      'No eligible variants after clonal-site filtering',
      gp = gpar(fontsize = 18)
    )
    dev.off()
    return(invisible(NULL))
  }
  
  heatmap_font_size = 12
  
  ht_opt(
    simple_anno_size = unit(heatmap_font_size+4, 'pt'),
    fast_hclust = TRUE,
    legend_title_gp = gpar(fontsize = heatmap_font_size, fontface = 'bold'),
    legend_labels_gp = gpar(fontsize = heatmap_font_size),
    legend_grid_width = unit(5, 'mm'),
    legend_grid_height = unit(5, 'mm'),
    heatmap_row_names_gp = gpar(fontsize = heatmap_font_size-4),
    heatmap_column_names_gp = gpar(fontsize = heatmap_font_size+4),
    HEATMAP_LEGEND_PADDING = unit(0.2, 'cm'),
    ANNOTATION_LEGEND_PADDING = unit(0.2, 'cm')
  )
  

  v3_mat = v3_make_alt_ref_matrix(results_df)
  v3_myeloid_alt_ha = v3_new_anno_row_annotation(results_df)
  
  v3_ht = Heatmap(
    v3_mat,
    name = 'Support',
    width = v3_clonal_heatmap_body_width,
    height = min(v3_clonal_heatmap_row_height * nrow(results_df),unit(150,'mm')),
    col = c('-1' = 'grey85', '0' = 'white', '1' = 'black'),
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      labels = c('pure REF', 'No coverage', 'ALT+')
    ),
    top_annotation = v3_new_anno_top_annotation(
      heatmap_font_size, sparse = nrow(results_df) <= 2
    ),
    left_annotation = if (show_myeloid_alt_annotation) v3_myeloid_alt_ha else NULL,
    cluster_rows = TRUE,
    clustering_distance_rows = 'euclidean',
    clustering_method_rows = 'ward.D2',
    cluster_columns = FALSE,
    show_row_names = nrow(results_df) <= 70,
    row_names_side = 'left',
    row_dend_side = 'right',
    row_dend_width = unit(3, 'mm'),
    row_names_max_width = unit(5, 'cm'),
    use_raster = TRUE,
    raster_quality = 3,
    rect_gp = gpar(col = 'grey85')
  )
  
  pdf(filename, width = width, height = height, useDingbats = FALSE, bg = 'transparent')
  v3_drawn_ht = ComplexHeatmap::draw(
    v3_ht,
    heatmap_legend_side = 'right',
    annotation_legend_side = 'right',
    merge_legends = TRUE,
    legend_gap = if (nrow(results_df) <= 2) unit(7, 'mm') else unit(3.5, 'mm'),
    padding = unit(c(5, 15, 5, 5), 'mm'),
    background = 'transparent'
  )
  if (nrow(results_df) <= 2) {
    decorate_annotation('sqrt_title', {
      grid.text(
        'sqrt(N cells)', x = unit(0.5, 'npc'), y = unit(0.5, 'npc'),
        gp = gpar(fontsize = heatmap_font_size + 1, fontface = 'bold')
      )
    })
  }
  # Capture the displayed order from this canonical annotated support plot.
  # The metric plots reuse its dendrogram so leaf orientation stays exact.
  v3_clonal_target_row_order <<- ComplexHeatmap::row_order(v3_drawn_ht)
  v3_clonal_target_row_dendrogram <<- if (nrow(results_df) > 1L) {
    ComplexHeatmap::row_dend(v3_drawn_ht)
  } else {
    FALSE
  }
  dev.off()
}

v3_save_new_anno_heatmap(
  CP_joint_bam_clonal_variants_df,
  'v3_putative_clonal_ALT_REF_heatmap_new_anno.pdf',
  width = 10.5, height = 8.75
)

v3_save_alt_metric_heatmap = function(results_df, filename, metric, width, height) {
  
  if (nrow(results_df) == 0) {
    pdf(filename, width = width, height = height, useDingbats = FALSE, bg = 'transparent')
    grid.newpage()
    grid.text('No eligible variants after clonal-site filtering', gp = gpar(fontsize = 18))
    dev.off()
    return(invisible(NULL))
  }
  heatmap_font_size = 12
  metric_label = if (metric == 'vaf') 'ALT VAF' else 'ALT+ cell fraction'
  metric_mat = v3_make_alt_metric_matrix(results_df, metric)
  metric_max = max(metric_mat, na.rm = TRUE)
  if (!is.finite(metric_max)) metric_max = 0
  color_max = if (metric_max > 0) metric_max else 1
  metric_color_breaks = seq(0, color_max, length.out = 9)
  metric_colors = RColorBrewer::brewer.pal(9, 'PuBuGn')

  ht_opt(
    simple_anno_size = unit(heatmap_font_size + 4, 'pt'), fast_hclust = TRUE,
    legend_title_gp = gpar(fontsize = heatmap_font_size, fontface = 'bold'),
    legend_labels_gp = gpar(fontsize = heatmap_font_size),
    legend_grid_width = unit(5, 'mm'), legend_grid_height = unit(5, 'mm'),
    heatmap_row_names_gp = gpar(fontsize = heatmap_font_size - 4),
    heatmap_column_names_gp = gpar(fontsize = heatmap_font_size + 4),
    HEATMAP_LEGEND_PADDING = unit(0.2, 'cm'),
    ANNOTATION_LEGEND_PADDING = unit(0.2, 'cm')
  )
  metric_ht = Heatmap(
    metric_mat, name = metric_label,
    width = v3_clonal_heatmap_body_width,
    height = v3_clonal_heatmap_row_height * nrow(results_df),
    col = circlize::colorRamp2(metric_color_breaks, metric_colors),
    na_col = 'white',
    heatmap_legend_param = list(
      at = c(0, metric_max), labels = c('0', formatC(metric_max, format = 'f', digits = 3))
    ),
    top_annotation = v3_new_anno_top_annotation(
      heatmap_font_size, sparse = nrow(results_df) <= 2
    ),
    left_annotation = v3_new_anno_row_annotation(results_df),
    cluster_rows = v3_clonal_target_row_dendrogram,
    row_dend_reorder = FALSE,
    cluster_columns = FALSE, show_row_names = nrow(results_df) <= 70,
    row_names_side = 'left', row_dend_side = 'right',
    row_dend_width = unit(3, 'mm'), row_names_max_width = unit(5, 'cm'),
    use_raster = TRUE, raster_quality = 3, rect_gp = gpar(col = 'grey85')
  )
  pdf(filename, width = width, height = height, useDingbats = FALSE, bg = 'transparent')
  metric_drawn_ht = ComplexHeatmap::draw(
    metric_ht, heatmap_legend_side = 'right', annotation_legend_side = 'right',
    merge_legends = TRUE,
    legend_gap = if (nrow(results_df) <= 2) unit(7, 'mm') else unit(3.5, 'mm'),
    padding = unit(c(5, 15, 5, 5), 'mm'), background = 'transparent'
  )
  if (nrow(results_df) <= 2) {
    decorate_annotation('sqrt_title', {
      grid.text(
        'sqrt(N cells)', x = unit(0.5, 'npc'), y = unit(0.5, 'npc'),
        gp = gpar(fontsize = heatmap_font_size + 1, fontface = 'bold')
      )
    })
  }
  if (!identical(
    as.integer(unlist(ComplexHeatmap::row_order(metric_drawn_ht))),
    as.integer(unlist(v3_clonal_target_row_order))
  )) {
    stop('Metric heatmap row order does not match the new-annotation support heatmap.')
  }
  dev.off()
}

v3_save_alt_metric_heatmap(
  CP_joint_bam_clonal_variants_df,
  'v3_putative_clonal_ALT_VAF_heatmap_new_anno.pdf',
  metric = 'vaf', width = 10.5, height = 8.75
)
v3_save_alt_metric_heatmap(
  CP_joint_bam_clonal_variants_df,
  'v3_putative_clonal_ALT_cell_fraction_heatmap_new_anno.pdf',
  metric = 'cell_fraction', width = 10.5, height = 8.75
)

