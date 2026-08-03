write_annotation_division_fixture <- function(root) {
  leaf <- file.path(root, "String", "Louvain", "C1")
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  nodes <- data.frame(node = sprintf("N%02d", 1:10), module = 1L)
  edges <- data.frame(
    node1 = sprintf("N%02d", 1:9),
    node2 = sprintf("N%02d", 2:10),
    module = 1L
  )
  utils::write.table(
    nodes,
    file.path(leaf, "node_Module_String_Louvain.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    edges,
    file.path(leaf, "edges_String_Louvain.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

annotation_test_gene_sets <- function() {
  data.frame(
    database = rep("TEST", 20),
    category = rep("TEST_CATEGORY", 20),
    term_id = rep(c("MODULE_TERM", "BACKGROUND_TERM"), each = 10),
    description = rep(c("Module genes", "Background genes"), each = 10),
    gene = c(sprintf("N%02d", 1:10), sprintf("N%02d", 11:20)),
    stringsAsFactors = FALSE
  )
}

test_that("step 5 annotates exactly the modules listed by step 4", {
  division <- tempfile("annotation-division-")
  selection <- tempfile("annotation-selection-")
  output <- tempfile("annotation-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_annotation_division_fixture(division)

  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  annotation <- functional_annotation(
    base_input_path = selection,
    base_output_path = output,
    network_method = "String",
    module_method = "Louvain",
    gene_sets = annotation_test_gene_sets(),
    minGSSize = 2,
    maxGSSize = 20,
    strict = TRUE
  )

  expect_equal(unique(annotation$module_uid), manifest$module_uid)
  expect_identical(annotation$term_id, "MODULE_TERM")
  expect_identical(annotation$gene_count, 10L)
  expect_identical(annotation$analysis_status, "completed")
  expect_true(file.exists(attr(annotation, "annotation_file")))
  expect_true(file.exists(attr(annotation, "top_file")))
  expect_true(file.exists(attr(annotation, "qc_file")))
  expect_true(file.exists(file.path(
    output,
    "modules",
    manifest$module_uid,
    "module_annotation.tsv"
  )))
})

test_that("step 5 records a status row when no term is significant", {
  division <- tempfile("annotation-division-")
  selection <- tempfile("annotation-selection-")
  output <- tempfile("annotation-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_annotation_division_fixture(division)
  module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )

  annotation <- functional_annotation(
    module_manifest_file = file.path(selection, "standardized", "module_manifest.tsv"),
    base_output_path = output,
    network_method = "String",
    module_method = "Louvain",
    gene_sets = annotation_test_gene_sets(),
    minGSSize = 2,
    maxGSSize = 20,
    pvalueCutoff = 1e-12,
    pAdjustCutoff = 1e-12
  )

  expect_equal(nrow(annotation), 1L)
  expect_identical(annotation$analysis_status, "completed_no_significant_terms")
  expect_identical(annotation$database, "NONE")
  top <- utils::read.delim(attr(annotation, "top_file"), check.names = FALSE)
  expect_equal(nrow(top), 0L)
})

test_that("step 5 accepts a custom tabular gene-set file", {
  division <- tempfile("annotation-division-")
  selection <- tempfile("annotation-selection-")
  output <- tempfile("annotation-output-")
  gene_set_file <- tempfile(fileext = ".tsv")
  dir.create(division)
  on.exit(unlink(c(division, selection, output, gene_set_file), recursive = TRUE, force = TRUE), add = TRUE)
  write_annotation_division_fixture(division)
  module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  utils::write.table(
    annotation_test_gene_sets(),
    gene_set_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  annotation <- functional_annotation(
    base_input_path = selection,
    base_output_path = output,
    network_method = "String",
    module_method = "Louvain",
    gene_set_file = gene_set_file,
    minGSSize = 2,
    maxGSSize = 20
  )
  expect_identical(annotation$term_id, "MODULE_TERM")
})

test_that("pipeline runner chains step 4 output directly into step 5", {
  division <- tempfile("annotation-division-")
  selection <- tempfile("annotation-selection-")
  output <- tempfile("annotation-output-")
  gene_set_file <- tempfile(fileext = ".tsv")
  config_file <- tempfile(fileext = ".yml")
  dir.create(division)
  on.exit(
    unlink(c(division, selection, output, gene_set_file, config_file), recursive = TRUE, force = TRUE),
    add = TRUE
  )
  write_annotation_division_fixture(division)
  utils::write.table(
    annotation_test_gene_sets(),
    gene_set_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  config <- yaml::read_yaml(system.file("config", "default.yml", package = "MLSnpDR"))
  config$project$target_subtypes <- "C1"
  config$pipeline$network_methods <- "String"
  config$pipeline$module_methods <- "Louvain"
  config$paths$module_division_dir <- division
  config$paths$module_selection_dir <- selection
  config$paths$module_annotation_dir <- output
  config$paths$annotation_gene_set_file <- gene_set_file
  config$annotation$databases <- "TEST"
  config$annotation$min_gene_set_size <- 2L
  config$annotation$max_gene_set_size <- 20L
  yaml::write_yaml(config, config_file)

  run <- run_ml_snpdr_pipeline(
    config_file,
    stages = c("module_selection", "module_annotation"),
    dry_run = FALSE
  )

  expect_identical(run$status, "completed")
  expect_equal(nrow(run$outputs$module_selection), 1L)
  expect_identical(run$outputs$module_annotation$term_id, "MODULE_TERM")
  expect_true(file.exists(attr(run$outputs$module_annotation, "annotation_file")))
})
