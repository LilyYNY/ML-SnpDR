write_feature_division_fixture <- function(root) {
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

feature_fixture_table <- function(module_id = "String|Louvain|C1|M1", module_size = 10) {
  result <- data.frame(
    module_id = module_id,
    module_label = "M1_C1",
    Subtype = "C1",
    Network = "String",
    Method = "Louvain",
    module_size = module_size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (index in 2:34) result[[sprintf("feature_%02d", index)]] <- index / 10
  result
}

test_that("step 6A creates a fixed all-module machine-learning matrix", {
  division <- tempfile("feature-division-")
  selection <- tempfile("feature-selection-")
  feature_file <- tempfile(fileext = ".tsv")
  output <- tempfile("feature-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, feature_file, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_feature_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  utils::write.table(
    feature_fixture_table(),
    feature_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  features <- prepare_module_features(
    module_manifest_file = attr(manifest, "manifest_file"),
    feature_file = feature_file,
    output_dir = output
  )

  expect_equal(nrow(features), 1L)
  expect_identical(features$module_uid, manifest$module_uid)
  expect_identical(features$module_size, 10)
  expect_identical(features$feature_missing_count, 0L)
  expect_true(features$feature_qc_pass)
  expect_true(file.exists(attr(features, "feature_file")))
  expect_true(file.exists(attr(features, "schema_file")))
  expect_true(file.exists(attr(features, "qc_file")))
  schema <- jsonlite::read_json(attr(features, "schema_file"), simplifyVector = TRUE)
  expect_identical(schema$feature_count, 34L)
  expect_identical(schema$feature_columns[[1L]], "module_size")
})

test_that("step 6A rejects missing manifest modules and duplicate identities", {
  division <- tempfile("feature-division-")
  selection <- tempfile("feature-selection-")
  feature_file <- tempfile(fileext = ".tsv")
  output <- tempfile("feature-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, feature_file, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_feature_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  duplicate <- rbind(feature_fixture_table(), feature_fixture_table())
  utils::write.table(duplicate, feature_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_module_features(
      module_manifest_file = attr(manifest, "manifest_file"),
      feature_file = feature_file,
      output_dir = output
    ),
    "duplicate module identities"
  )
})

test_that("step 6A strict QC rejects non-finite values and module-size mismatch", {
  division <- tempfile("feature-division-")
  selection <- tempfile("feature-selection-")
  feature_file <- tempfile(fileext = ".tsv")
  output <- tempfile("feature-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, feature_file, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_feature_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  bad <- feature_fixture_table(module_size = 11)
  bad$feature_02 <- NA_real_
  utils::write.table(bad, feature_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_module_features(
      module_manifest_file = attr(manifest, "manifest_file"),
      feature_file = feature_file,
      output_dir = output
    ),
    "Feature QC failed"
  )
  expect_false(dir.exists(output))
})
