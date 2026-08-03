write_upstream_fixture <- function(root) {
  samples <- c("C1_A", "C1_B", "C1_C", "C2_A", "C2_B", "C2_C")
  expression <- data.frame(
    gene = paste0("G", 1:10),
    C1_A = c(rep(10, 5), rep(1, 5)),
    C1_B = c(rep(11, 5), rep(1, 5)),
    C1_C = c(rep(9, 5), rep(1, 5)),
    C2_A = c(rep(1, 5), rep(10, 5)),
    C2_B = c(rep(1, 5), rep(11, 5)),
    C2_C = c(rep(1, 5), rep(9, 5)),
    check.names = FALSE
  )
  phenotype <- data.frame(
    Sample = samples,
    Subtype = rep(c("C1", "C2"), each = 3),
    stringsAsFactors = FALSE
  )
  clique <- function(nodes) {
    pairs <- utils::combn(nodes, 2)
    data.frame(node1 = pairs[1, ], node2 = pairs[2, ], stringsAsFactors = FALSE)
  }
  ppi <- rbind(clique(paste0("G", 1:5)), clique(paste0("G", 6:10)))
  expression_file <- file.path(root, "expression.tsv")
  phenotype_file <- file.path(root, "subtype.tsv")
  ppi_file <- file.path(root, "string.tsv")
  index_file <- file.path(root, "ppi_index.tsv")
  utils::write.table(expression, expression_file, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(phenotype, phenotype_file, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(ppi, ppi_file, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(
    data.frame(network = "String", edge_file = basename(ppi_file)),
    index_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(
    expression_file = expression_file,
    phenotype_file = phenotype_file,
    ppi_file = ppi_file,
    index_file = index_file
  )
}

test_that("steps 1-3 produce directly chainable canonical outputs", {
  root <- tempfile("upstream-chain-")
  dir.create(root)
  files <- write_upstream_fixture(root)
  dep_dir <- file.path(root, "01_dep")
  network_dir <- file.path(root, "02_network")
  division_dir <- file.path(root, "03_division")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  dep <- run_diff_expr_analysis(
    expression_file = files$expression_file,
    phenotype_file = files$phenotype_file,
    output_dir = dep_dir,
    gene_column = "gene",
    p_threshold = 1,
    write_legacy_excel = FALSE,
    subtypes = c("C1", "C2")
  )
  expect_equal(nrow(dep), 20L)
  expect_equal(sum(dep$label == "up"), 10L)
  expect_true(file.exists(attr(dep, "result_file")))

  networks <- run_network_construction(
    diff_file = attr(dep, "result_file"),
    ppi_index_file = files$index_file,
    output_dir = network_dir,
    ppi_method = "String",
    subtypes = c("C1", "C2")
  )
  expect_equal(nrow(networks), 2L)
  expect_true(all(networks$node_count == 5L))
  expect_true(all(networks$edge_count == 10L))
  expect_true(file.exists(attr(networks, "manifest_file")))

  divisions <- module_division(
    output_base_path = division_dir,
    network_method = "String",
    module_method = "Louvain",
    network_manifest_file = attr(networks, "manifest_file"),
    subtypes = c("C1", "C2")
  )
  expect_equal(nrow(divisions), 2L)
  expect_true(all(divisions$module_count == 1L))
  expect_true(file.exists(attr(divisions, "manifest_file")))
  expect_true(all(file.exists(file.path(division_dir, divisions$node_module_file))))
})

test_that("the configured runner chains steps 1 through 4", {
  root <- tempfile("runner-upstream-")
  dir.create(root)
  files <- write_upstream_fixture(root)
  config_file <- file.path(root, "config.yml")
  config <- yaml::read_yaml(system.file("config", "default.yml", package = "MLSnpDR"))
  config$project$target_subtypes <- c("C1", "C2")
  config$pipeline$network_methods <- "String"
  config$pipeline$module_methods <- "Louvain"
  config$paths$expression_file <- files$expression_file
  config$paths$subtype_file <- files$phenotype_file
  config$paths$ppi_index_file <- files$index_file
  config$paths$differential_expression_dir <- file.path(root, "01_dep")
  config$paths$network_construction_dir <- file.path(root, "02_network")
  config$paths$module_division_dir <- file.path(root, "03_division")
  config$paths$module_selection_dir <- file.path(root, "04_selection")
  config$differential_expression$p_adjust_threshold <- 1
  config$differential_expression$write_legacy_excel <- TRUE
  config$module_prefilter$min_size_exclusive <- 2
  yaml::write_yaml(config, config_file)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  run <- run_ML_SnpDR(
    config_file,
    from = "deps",
    to = "module_selection",
    dry_run = FALSE
  )
  expect_identical(run$status, "completed")
  expect_identical(names(run$outputs), c("deps", "network_construction", "module_division", "module_selection"))
  expect_equal(nrow(run$outputs$module_selection), 2L)
  expect_true(file.exists(attr(run$outputs$module_selection, "manifest_file")))
})

test_that("step 2 applies score filters and explicit ID mappings", {
  root <- tempfile("ppi-mapping-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  diff <- data.frame(
    subtype = rep("C1", 3),
    gene = c("G1", "G2", "G3"),
    label = rep("up", 3)
  )
  raw_edges <- data.frame(
    protein1 = c("P1", "P1"),
    protein2 = c("P2", "P3"),
    combined_score = c(900, 200)
  )
  mapping <- data.frame(protein_id = c("P1", "P2", "P3"), symbol = c("G1", "G2", "G3"))
  utils::write.table(diff, file.path(root, "diff.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(raw_edges, file.path(root, "raw.txt"), sep = " ", quote = FALSE, row.names = FALSE)
  utils::write.table(mapping, file.path(root, "map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  index <- data.frame(
    network = "String",
    edge_file = "raw.txt",
    node1_column = "protein1",
    node2_column = "protein2",
    score_column = "combined_score",
    score_min = 400,
    delimiter = "whitespace",
    mapping_file = "map.tsv",
    mapping_id_column = "protein_id",
    mapping_symbol_column = "symbol",
    mapping_delimiter = "tab"
  )
  utils::write.table(index, file.path(root, "index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  result <- run_network_construction(
    diff_file = file.path(root, "diff.tsv"),
    ppi_index_file = file.path(root, "index.tsv"),
    output_dir = file.path(root, "networks"),
    ppi_method = "String",
    subtypes = "C1"
  )
  edges <- utils::read.delim(file.path(attr(result, "output_dir"), result$ppi_file))
  expect_equal(nrow(edges), 1L)
  expect_identical(edges$node1, "G1")
  expect_identical(edges$node2, "G2")
})

test_that("WF produces deterministic disjoint module assignments", {
  root <- tempfile("wf-module-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  clique <- function(nodes) {
    pairs <- utils::combn(nodes, 2)
    data.frame(node1 = pairs[1, ], node2 = pairs[2, ], stringsAsFactors = FALSE)
  }
  ppi <- rbind(
    clique(paste0("A", 1:6)),
    clique(paste0("B", 1:6)),
    data.frame(node1 = "A1", node2 = "B1")
  )
  ppi_file <- file.path(root, "ppi.tsv")
  output <- file.path(root, "out")
  utils::write.table(ppi, ppi_file, sep = "\t", quote = FALSE, row.names = FALSE)
  result <- subtype_module(ppi_file, output, "String", "WF", seed = 123)
  nodes <- utils::read.delim(result$node_module_file_abs, stringsAsFactors = FALSE)
  expect_identical(anyDuplicated(nodes$node), 0L)
  expect_gte(result$module_count, 1L)
})
