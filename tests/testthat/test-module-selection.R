write_division_fixture <- function(
  root,
  network = "String",
  method = "Louvain",
  subtype = "C1"
) {
  leaf <- file.path(root, network, method, subtype)
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

  nodes <- data.frame(
    node = c(sprintf("N%02d", 1:10), "X01", "X02"),
    module = c(rep(1L, 10), 2L, 2L)
  )
  edges <- data.frame(
    node1 = c(sprintf("N%02d", 1:9), "N10", "N01", "X01"),
    node2 = c(sprintf("N%02d", 2:10), "N10", "N10", "X02"),
    module = c(rep(1L, 10), 0L, 2L)
  )
  utils::write.table(
    nodes,
    file.path(leaf, paste0("node_Module_", network, "_", method, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    edges,
    file.path(leaf, paste0("edges_", network, "_", method, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

test_that("step 4 preserves subnetDR outputs and returns the canonical manifest", {
  input <- tempfile("module-division-")
  output <- tempfile("module-selection-")
  dir.create(input)
  on.exit(unlink(c(input, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_division_fixture(input)

  manifest <- module_selection(
    base_input_path = input,
    base_output_path = output,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )

  legacy_dir <- file.path(output, "String", "Louvain", "C1")
  expect_true(file.exists(file.path(legacy_dir, "Module_Number_String_C1.csv")))
  expect_true(file.exists(file.path(legacy_dir, "Module_select_String_C1.csv")))
  expect_true(file.exists(file.path(legacy_dir, "node_Module_select_String_C1.txt")))
  expect_true(file.exists(file.path(legacy_dir, "edges_select_String_C1.txt")))

  expect_equal(nrow(manifest), 1L)
  expect_identical(manifest$module_uid, "string__louvain__C1__M1")
  expect_identical(manifest$module_size, 10L)
  expect_identical(manifest$edge_count, 9L)
  expect_identical(manifest$self_loops_excluded, 1L)
  expect_true(file.exists(attr(manifest, "manifest_file")))
  expect_true(file.exists(attr(manifest, "qc_file")))

  selected_edges <- utils::read.delim(file.path(legacy_dir, "edges_select_String_C1.txt"))
  expect_equal(nrow(selected_edges), 11L)
  expect_true(any(selected_edges$module == "M0"))
})

test_that("pipeline runner can execute implemented step 4 with custom paths", {
  input <- tempfile("module-division-")
  output <- tempfile("module-selection-")
  config_file <- tempfile(fileext = ".yml")
  dir.create(input)
  on.exit(unlink(c(input, output, config_file), recursive = TRUE, force = TRUE), add = TRUE)
  write_division_fixture(input)

  config <- yaml::read_yaml(system.file("config", "default.yml", package = "MLSnpDR"))
  config$project$target_subtypes <- "C1"
  config$pipeline$network_methods <- "String"
  config$pipeline$module_methods <- "Louvain"
  config$paths$module_division_dir <- input
  config$paths$module_selection_dir <- output
  yaml::write_yaml(config, config_file)

  run <- run_ml_snpdr_pipeline(
    config_file,
    stages = "module_selection",
    dry_run = FALSE
  )

  expect_identical(run$status, "completed")
  expect_equal(nrow(run$outputs$module_selection), 1L)
  expect_true(file.exists(attr(run$outputs$module_selection, "manifest_file")))
})

test_that("step 4 refuses to overwrite an existing output directory", {
  input <- tempfile("module-division-")
  output <- tempfile("module-selection-")
  dir.create(input)
  dir.create(output)
  on.exit(unlink(c(input, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_division_fixture(input)

  expect_error(
    module_selection(
      base_input_path = input,
      base_output_path = output,
      network_method = "String",
      module_method = "Louvain",
      subtypes = "C1"
    ),
    "already exists"
  )
})
