nextflow.enable.dsl = 2

process PREPARE_METADATA {
  tag "${donor}"
  input:
  tuple val(donor), val(donor_id), val(reference_key)
  output:
  tuple val(donor), val(donor_id), val(reference_key), path("${donor}.prepare.ok"), emit: ready
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' \\
  run_v3_task.sh prepare '${donor}' '${donor_id}' '${reference_key}' '${params.manifest}' '${params.cell_metadata}' '${params.reference_config}'
  """
}

process STEP1_SPLIT {
  tag "${donor}:${source_id}"
  input:
  tuple val(donor), val(donor_id), val(reference_key), val(tissue), val(source_id), val(meta_prefix), val(bam)
  output:
  tuple val(donor), val(reference_key), val(source_id), path("${donor}.${source_id}.step1.ok"), emit: ready
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' \\
  run_v3_task.sh step1 '${donor}' '${tissue}' '${source_id}' '${bam}'
  """
}

process BUILD_JOINT_BAMS {
  tag "${donor}"
  input:
  tuple val(donor), val(reference_keys), val(source_ids), path(markers), val(recovery_revision)
  output:
  tuple val(donor), path("${donor}.joint.ok"), emit: ready
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' CELL_TYPES='${params.cell_types}' \\
  run_v3_task.sh joint '${donor}'
  """
}

process STEP2_COUNTER {
  tag "${donor}:${cell_type}"
  input:
  tuple val(donor), val(cell_type), val(recovery_revision)
  output:
  tuple val(donor), val(cell_type), path("${donor}.${cell_type}.step2.ok"), emit: ready
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' NPROCS='${task.cpus}' MIN_BASE_QUALITY='${params.min_base_quality}' \\
  run_v3_task.sh step2 '${donor}' '${cell_type}'
  """
}

process MERGE_STEP3 {
  tag "${donor}"
  input:
  tuple val(donor), val(cell_types), val(recovery_revisions), path(markers)
  output:
  tuple val(donor), val(recovery_revisions), path("${donor}.merge.ok"), emit: ready
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' \\
  run_v3_task.sh merge '${donor}'
  """
}

process STEP4_MAX15_POSTPROCESS {
  tag "${donor}"
  input:
  tuple val(donor), val(recovery_revisions), path(marker)
  output:
  path("${donor}.step4.ok")
  script:
  """
  RUNROOT='${params.run_root}' SCOMATIC_ENV='${params.scomatic_env}' SCOMATIC_REPO='${params.scomatic_repo}' MAX_CELL_TYPES='${params.max_cell_types}' \\
  run_v3_task.sh step4 '${donor}'
  """
}

process GNOMAD_POPMAX_FILTER {
  tag 'all donors'
  input:
  path(markers)
  output:
  path('gnomad_popmax.ok')
  script:
  """
  RUNROOT='${params.run_root}' GNOMAD_AF_CUTOFF='${params.gnomad_af_cutoff}' GNOMAD_RELEASE='${params.gnomad_release}' \\
  run_v3_task.sh gnomad_popmax '${params.manifest}'
  touch gnomad_popmax.ok
  """
}

workflow {
  if (!params.manifest) error 'Missing required parameter: --manifest'
  if (!params.cell_metadata) error 'Missing required parameter: --cell_metadata'
  if (!params.reference_config) error 'Missing required parameter: --reference_config'
  if (!params.scomatic_env) error 'Missing required parameter: --scomatic_env'
  if (!params.scomatic_repo) error 'Missing required parameter: --scomatic_repo'
  file(params.manifest, checkIfExists: true)
  file(params.cell_metadata, checkIfExists: true)
  file(params.reference_config, checkIfExists: true)
  file(params.scomatic_env, checkIfExists: true)
  file(params.scomatic_repo, checkIfExists: true)
  if (!params.run_root) error 'Missing required parameter: --run_root'

  libraries = Channel
    .fromPath(params.manifest)
    .splitCsv(header: true, sep: '\t')
    .map { r -> tuple(r.donor_safe, r.tissue, r.source_id, r.meta_prefix, r.bam) }

  donors = Channel
    .fromPath(params.manifest)
    .splitCsv(header: true, sep: '\t')
    .map { r -> tuple(r.donor_safe, r.donor_id, r.reference_key) }
    .unique()

  prep = PREPARE_METADATA(donors)
  step1In = prep.ready
    .map { donor, donorId, referenceKey, marker -> tuple(donor, donorId, referenceKey) }
    .combine(libraries, by: 0)
  step1 = STEP1_SPLIT(step1In)

  // One CSF and one PBMC library per donor. Joint work starts as soon as that
  // donor's two Step 1 jobs finish; it does not wait for other donors.
  jointIn = step1.ready
    .groupTuple(size: 2)
    .map { donor, referenceKeys, sourceIds, markers -> tuple(donor, referenceKeys, sourceIds, markers, params.recovery_revision) }
  joint = BUILD_JOINT_BAMS(jointIn)

  cellTypes = params.cell_types.tokenize(',')*.trim()
  step2In = joint.ready.flatMap { donor, marker -> cellTypes.collect { cellType -> [donor, cellType, params.recovery_revision] } }
  step2 = STEP2_COUNTER(step2In)

  // Fixed per-donor grouping prevents one donor's Step 3/4 path from being
  // blocked by Step 2 work for another donor.
  mergeIn = step2.ready
    .groupTuple(size: cellTypes.size())
    .map { donor, groupedCellTypes, markers -> tuple(donor, groupedCellTypes, params.recovery_revision, markers) }
  merged = MERGE_STEP3(mergeIn)
  step4 = STEP4_MAX15_POSTPROCESS(merged.ready)
  GNOMAD_POPMAX_FILTER(step4.collect())
}
