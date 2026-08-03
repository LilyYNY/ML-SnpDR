write_manifest_fixture <- function(
  root,
  network = "String",
  method = "Louvain",
  subtype = "C1",
  declared_count = 3L,
  empty = FALSE
) {
  leaf <- file.path(root, network, method, subtype)
  dir.create(leaf, recursive = TRUE, showWarnings = FALSE)

  if (empty) {
    modules <- data.frame(module = character(), count = integer())
    nodes <- data.frame(node = character(), module = character())
    edges <- data.frame(node1 = character(), node2 = character(), module = character())
  } else {
    modules <- data.frame(module = "M1", count = declared_count)
    nodes <- data.frame(node = c("C", "A", "B"), module = "M1")
    edges <- data.frame(
      node1 = c("B", "A", "C", "A"),
      node2 = c("C", "B", "C", "outside"),
      module = c("M1", "M1", "M1", "M0")
    )
  }

  utils::write.csv(
    modules,
    file.path(leaf, paste0("Module_select_", network, "_", subtype, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    nodes,
    file.path(leaf, paste0("node_Module_select_", network, "_", subtype, ".txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    edges,
    file.path(leaf, paste0("edges_select_", network, "_", subtype, ".txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

test_that("ModuleSelection adapter writes a canonical manifest and module files", {
  input <- tempfile("module-selection-")
  output <- tempfile("module-manifest-")
  dir.create(input)
  on.exit(unlink(c(input, output), recursive = TRUE, force = TRUE), add = TRUE)

  write_manifest_fixture(input)
  write_manifest_fixture(input, network = "chengF", method = "WF", subtype = "C2", empty = TRUE)

  manifest <- build_module_manifest(input, output, min_size_exclusive = 1)

  expect_s3_class(manifest, "data.frame")
  expect_equal(nrow(manifest), 1L)
  expect_identical(manifest$module_uid, "string__louvain__C1__M1")
  expect_identical(manifest$module_size, 3L)
  expect_identical(manifest$edge_count, 2L)
  expect_identical(manifest$self_loops_excluded, 1L)
  expect_true(manifest$prefilter_pass)
  expect_match(manifest$node_sha256, "^[0-9a-f]{64}$")
  expect_match(manifest$edge_sha256, "^[0-9a-f]{64}$")
  expect_false(grepl("^[A-Za-z]:", manifest$node_file))

  manifest_file <- file.path(output, "module_manifest.tsv")
  qc_file <- file.path(output, "module_manifest_qc.tsv")
  expect_true(file.exists(manifest_file))
  expect_true(file.exists(qc_file))

  nodes <- utils::read.delim(file.path(output, manifest$node_file))
  edges <- utils::read.delim(file.path(output, manifest$edge_file))
  qc <- utils::read.delim(qc_file)
  expect_identical(nodes$node, c("A", "B", "C"))
  expect_equal(nrow(edges), 2L)
  expect_false(any(edges$module == "M0"))
  expect_equal(nrow(qc), 2L)
  expect_identical(qc$n_modules[qc$network == "chengF"], 0L)
  expect_identical(qc$excluded_self_loops[qc$network == "String"], 1L)
  expect_identical(qc$excluded_m0_edges[qc$network == "String"], 1L)

  reloaded <- read_module_manifest(manifest_file, verify_hashes = TRUE)
  expect_equal(nrow(reloaded), 1L)
  expect_true(file.exists(reloaded$node_file_abs))
  expect_true(file.exists(reloaded$edge_file_abs))
})

test_that("ModuleSelection adapter rejects declared-size mismatches", {
  input <- tempfile("module-selection-")
  output <- tempfile("module-manifest-")
  dir.create(input)
  on.exit(unlink(c(input, output), recursive = TRUE, force = TRUE), add = TRUE)

  write_manifest_fixture(input, declared_count = 4L)

  expect_error(
    build_module_manifest(input, output, min_size_exclusive = 1),
    "Declared size does not match node count"
  )
  expect_false(dir.exists(output))
})

test_that("ModuleSelection adapter never overwrites an output directory", {
  input <- tempfile("module-selection-")
  output <- tempfile("module-manifest-")
  dir.create(input)
  dir.create(output)
  on.exit(unlink(c(input, output), recursive = TRUE, force = TRUE), add = TRUE)

  write_manifest_fixture(input)

  expect_error(
    build_module_manifest(input, output, min_size_exclusive = 1),
    "output_dir already exists"
  )
})

test_that("installed synthetic ModuleSelection example runs end to end", {
  input <- system.file("extdata", "example", "ModuleSelection", package = "MLSnpDR")
  output <- tempfile("module-manifest-example-")
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)

  expect_true(nzchar(input))
  manifest <- build_module_manifest(input, output)

  expect_equal(nrow(manifest), 1L)
  expect_identical(manifest$module_size, 10L)
  expect_identical(manifest$edge_count, 9L)
  expect_identical(manifest$self_loops_excluded, 1L)
})
