write_selected_handoff_fixture <- function(root) {
  uid <- "string__louvain__C1__M1"
  relative_dir <- file.path("selected", "C1", uid)
  module_dir <- file.path(root, relative_dir)
  dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(node = c("N01", "N02"), module = "M1"),
    file.path(module_dir, "nodes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node1 = "N01", node2 = "N02", module = "M1"),
    file.path(module_dir, "edges.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      protein = c("N01", "N02"),
      drug = c("DrugA", "DrugB"),
      p_value = c(0.01, 0.02),
      p_adjust = c(0.01, 0.02),
      significant = TRUE
    ),
    file.path(module_dir, "drn.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node = c("N01", "N02", "DrugA", "DrugB"), type = c("protein", "protein", "drug", "drug")),
    file.path(module_dir, "drn_info.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  selected <- data.frame(
    module_uid = uid,
    network = "String",
    method = "Louvain",
    subtype = "C1",
    module = "M1",
    module_size = 2L,
    primary_drug_panel = "PRISM",
    node_file = file.path(relative_dir, "nodes.tsv"),
    edge_file = file.path(relative_dir, "edges.tsv"),
    drn_file = file.path(relative_dir, "drn.tsv"),
    drn_info_file = file.path(relative_dir, "drn_info.tsv"),
    selection_rank = 1L,
    selection_reason = "test",
    stringsAsFactors = FALSE
  )
  selected_file <- file.path(root, "selected_modules.tsv")
  utils::write.table(selected, selected_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(uid = uid, file = selected_file)
}

test_that("steps 7-9 chain selected modules without directory rediscovery", {
  selected_root <- tempfile("selected-handoff-")
  seq_output <- tempfile("seq-output-")
  binding_output <- tempfile("binding-output-")
  perturbation_output <- tempfile("perturbation-output-")
  protein_file <- tempfile(fileext = ".tsv")
  smiles_file <- tempfile(fileext = ".tsv")
  dir.create(selected_root)
  on.exit(
    unlink(
      c(selected_root, seq_output, binding_output, perturbation_output, protein_file, smiles_file),
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )
  selected <- write_selected_handoff_fixture(selected_root)
  utils::write.table(
    data.frame(node = c("N01", "N02"), sequence = c("MAAAA", "MBBBB")),
    protein_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node = c("DrugA", "DrugB"), SMILES = c("CCO", "CCN")),
    smiles_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  seq_manifest <- run_SEQCre(
    input_base = selected$file,
    output_base = seq_output,
    protein_sequence_file = protein_file,
    drug_smiles_file = smiles_file
  )
  expect_equal(nrow(seq_manifest), 1L)
  expect_identical(seq_manifest$module_uid, selected$uid)
  expect_identical(seq_manifest$matched_protein_number, 2L)
  expect_identical(seq_manifest$matched_drug_number, 2L)

  binding <- predict_BA(
    selected_modules_file = selected$file,
    seq_smiles_manifest_file = attr(seq_manifest, "manifest_file"),
    output_base = binding_output,
    predictor = function(dpi) seq_len(nrow(dpi)) / 10
  )
  expect_equal(nrow(binding), 2L)
  expect_true(all(binding$module_uid == selected$uid))
  expect_identical(binding$binding_rank, 1:2)

  perturbation <- process_prs_dti(
    selected_modules_file = selected$file,
    binding_scores_file = attr(binding, "scores_file"),
    output_base = perturbation_output,
    sensitivity_function = function(edges, module) {
      data.frame(
        target_name = sort(unique(c(edges$node1, edges$node2))),
        sensitivity = c(2, 3),
        stringsAsFactors = FALSE
      )
    },
    top_n = 2L
  )
  expect_equal(nrow(perturbation), 2L)
  expect_true(all(perturbation$module_uid == selected$uid))
  expect_equal(sort(perturbation$perturbation_score), c(0.2, 0.6))
  final <- utils::read.delim(attr(perturbation, "final_file"), check.names = FALSE)
  expect_equal(nrow(final), 2L)
  expect_identical(final$final_rank_in_subtype, 1:2)
})

test_that("step 7 strict mode rejects incomplete lookup coverage", {
  selected_root <- tempfile("selected-handoff-")
  output <- tempfile("seq-output-")
  protein_file <- tempfile(fileext = ".tsv")
  smiles_file <- tempfile(fileext = ".tsv")
  dir.create(selected_root)
  on.exit(unlink(c(selected_root, output, protein_file, smiles_file), recursive = TRUE, force = TRUE), add = TRUE)
  selected <- write_selected_handoff_fixture(selected_root)
  utils::write.table(
    data.frame(node = "N01", sequence = "MAAAA"),
    protein_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node = c("DrugA", "DrugB"), SMILES = c("CCO", "CCN")),
    smiles_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  expect_error(
    run_SEQCre(
      input_base = selected$file,
      output_base = output,
      protein_sequence_file = protein_file,
      drug_smiles_file = smiles_file,
      strict = TRUE
    ),
    "missing proteins=1"
  )
  expect_false(dir.exists(output))
})
