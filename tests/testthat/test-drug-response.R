write_drug_division_fixture <- function(root) {
  leaf <- file.path(root, "String", "Louvain", "C1")
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(node = sprintf("N%02d", 1:10), module = 1L),
    file.path(leaf, "node_Module_String_Louvain.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      node1 = sprintf("N%02d", 1:9),
      node2 = sprintf("N%02d", 2:10),
      module = 1L
    ),
    file.path(leaf, "edges_String_Louvain.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_drug_response_fixture <- function(root, panel = "PRISM") {
  module_dir <- file.path(root, "C1", "Louvain", "string", "M1")
  prediction_dir <- file.path(module_dir, "calcPhenotype_Output")
  dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- "louvain_string_C1_M1"
  utils::write.table(
    data.frame(
      protein = c("N01", "N02"),
      drug = c("DrugA", "DrugB"),
      pvalue = c(0.001, 0.02),
      label = "sig"
    ),
    file.path(module_dir, paste0("DRN_", prefix, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node = c("N01", "N02", "DrugA", "DrugB"), type = c("protein", "protein", "drug", "drug")),
    file.path(module_dir, paste0("DRN_info_", prefix, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(label = c("non-sig", "sig"), count = c(1L, 1L)),
    file.path(module_dir, paste0("drug_response_level_", prefix, ".csv")),
    quote = FALSE,
    row.names = FALSE
  )
  predictions <- data.frame(Sample = c("S1", "S2"), DrugA = c(0.1, 0.2), DrugB = c(0.3, 0.4))
  utils::write.csv(
    predictions,
    file.path(prediction_dir, "DrugPredictions.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  invisible(panel)
}

test_that("step 6 standardizes every manifest module for every panel", {
  division <- tempfile("drug-division-")
  selection <- tempfile("drug-selection-")
  prism <- tempfile("drug-prism-")
  gdsc1 <- tempfile("drug-gdsc1-")
  output <- tempfile("drug-output-")
  dir.create(division)
  dir.create(prism)
  dir.create(gdsc1)
  on.exit(unlink(c(division, selection, prism, gdsc1, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_drug_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  write_drug_response_fixture(prism)
  write_drug_response_fixture(gdsc1)

  summary <- drug_response_analysis(
    module_manifest_file = attr(manifest, "manifest_file"),
    drug_response_path = output,
    drug_response_roots = c(PRISM = prism, GDSC1 = gdsc1),
    strict = TRUE
  )

  expect_equal(nrow(summary), 2L)
  expect_identical(summary$drug_panel, c("GDSC1", "PRISM"))
  expect_true(all(summary$module_uid == manifest$module_uid))
  expect_true(all(summary$drug_number == 1L))
  expect_true(all(summary$tested_drug_number == 2L))
  expect_true(all(summary$drn_edge_number == 2L))
  expect_true(all(summary$analysis_status == "completed"))
  expect_true(file.exists(attr(summary, "summary_file")))
  expect_true(file.exists(attr(summary, "hits_file")))
  expect_true(file.exists(attr(summary, "coverage_file")))
  expect_true(all(file.exists(file.path(output, summary$drn_file))))
  hits <- utils::read.delim(attr(summary, "hits_file"), check.names = FALSE)
  expect_equal(nrow(hits), 4L)
  expect_true(all(hits$evidence_type == "protein_drug_drn"))
})

test_that("step 6 accepts an explicit relative-path input index", {
  division <- tempfile("drug-division-")
  selection <- tempfile("drug-selection-")
  inputs <- tempfile("drug-inputs-")
  output <- tempfile("drug-output-")
  index_file <- tempfile(tmpdir = inputs, fileext = ".tsv")
  dir.create(division)
  dir.create(inputs)
  on.exit(unlink(c(division, selection, inputs, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_drug_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  write_drug_response_fixture(inputs)
  module_dir <- file.path("C1", "Louvain", "string", "M1")
  prefix <- "louvain_string_C1_M1"
  index <- data.frame(
    module_uid = manifest$module_uid,
    drug_panel = "PRISM",
    drn_file = file.path(module_dir, paste0("DRN_", prefix, ".txt")),
    drn_info_file = file.path(module_dir, paste0("DRN_info_", prefix, ".txt")),
    drug_level_file = file.path(module_dir, paste0("drug_response_level_", prefix, ".csv")),
    prediction_file = file.path(module_dir, "calcPhenotype_Output", "DrugPredictions.csv"),
    stringsAsFactors = FALSE
  )
  utils::write.table(index, index_file, sep = "\t", quote = FALSE, row.names = FALSE)

  summary <- drug_response_analysis(
    module_manifest_file = attr(manifest, "manifest_file"),
    drug_response_path = output,
    drug_response_index_file = index_file
  )
  expect_equal(nrow(summary), 1L)
  expect_identical(summary$drug_panel, "PRISM")
})

test_that("step 6 strict mode rejects incomplete panel inputs", {
  division <- tempfile("drug-division-")
  selection <- tempfile("drug-selection-")
  prism <- tempfile("drug-prism-")
  empty_panel <- tempfile("drug-empty-")
  output <- tempfile("drug-output-")
  dir.create(division)
  dir.create(prism)
  dir.create(empty_panel)
  on.exit(unlink(c(division, selection, prism, empty_panel, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_drug_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  write_drug_response_fixture(prism)

  expect_error(
    drug_response_analysis(
      module_manifest_file = attr(manifest, "manifest_file"),
      drug_response_path = output,
      drug_response_roots = c(PRISM = prism, GDSC1 = empty_panel),
      strict = TRUE
    ),
    "Missing drug-response input file"
  )
  expect_false(dir.exists(output))
})

test_that("pipeline runner chains step 4 output directly into step 6", {
  division <- tempfile("drug-division-")
  selection <- tempfile("drug-selection-")
  prism <- tempfile("drug-prism-")
  output <- tempfile("drug-output-")
  config_file <- tempfile(fileext = ".yml")
  dir.create(division)
  dir.create(prism)
  on.exit(
    unlink(c(division, selection, prism, output, config_file), recursive = TRUE, force = TRUE),
    add = TRUE
  )
  write_drug_division_fixture(division)
  write_drug_response_fixture(prism)

  config <- yaml::read_yaml(system.file("config", "default.yml", package = "MLSnpDR"))
  config$project$target_subtypes <- "C1"
  config$pipeline$network_methods <- "String"
  config$pipeline$module_methods <- "Louvain"
  config$paths$module_division_dir <- division
  config$paths$module_selection_dir <- selection
  config$paths$drug_response_dir <- output
  config$drug_response$panels <- "PRISM"
  config$drug_response$input_roots <- list(PRISM = prism)
  yaml::write_yaml(config, config_file)

  run <- run_ml_snpdr_pipeline(
    config_file,
    stages = c("module_selection", "drug_response"),
    dry_run = FALSE
  )

  expect_identical(run$status, "completed")
  expect_equal(nrow(run$outputs$module_selection), 1L)
  expect_equal(nrow(run$outputs$drug_response), 1L)
  expect_identical(run$outputs$drug_response$module_uid, run$outputs$module_selection$module_uid)
})
