write_triage_division_fixture <- function(root) {
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
    data.frame(node1 = sprintf("N%02d", 1:9), node2 = sprintf("N%02d", 2:10), module = 1L),
    file.path(leaf, "edges_String_Louvain.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_triage_inputs <- function(root, module_uid, survival_direction = "High_score_worse") {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  ml_top_file <- file.path(root, "ml_top1.tsv")
  survival_file <- file.path(root, "survival.tsv")
  drug_dir <- file.path(root, "drug")
  drn_dir <- file.path(drug_dir, "modules", module_uid, "PRISM")
  dir.create(drn_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(
      module_uid = module_uid,
      target_subtype = "C1",
      target_subtype_probability = 0.9,
      probability_margin = 0.8,
      rank_in_subtype = 1L,
      ml_gate = TRUE,
      ml_gate_reason = "passed_configured_ml_gates"
    ),
    ml_top_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      module_id = "String|Louvain|C1|M1",
      n_samples = 100L,
      n_events = 30L,
      matched_gene_fraction = 1,
      cox_HR_per_1SD = 1.5,
      cox_P = 0.01,
      cox_FDR = 0.04,
      optimal_cutpoint = 0.2,
      optimal_logrank_P = 0.01,
      optimal_logrank_FDR = 0.04,
      optimal_HR_high_vs_low = 2,
      optimal_direction = survival_direction,
      stringsAsFactors = FALSE
    ),
    survival_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(protein = "N01", drug = "DrugA", p_value = 0.01, p_adjust = 0.01, significant = TRUE),
    file.path(drn_dir, "drn.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(node = c("N01", "DrugA"), type = c("protein", "drug")),
    file.path(drn_dir, "drn_info.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  drug_summary_file <- file.path(drug_dir, "drug_response_summary.tsv")
  utils::write.table(
    data.frame(
      module_uid = module_uid,
      drug_panel = "PRISM",
      drug_number = 5L,
      tested_drug_number = 10L,
      drug_response_density = 0.5,
      drn_edge_number = 1L,
      drn_file = file.path("modules", module_uid, "PRISM", "drn.tsv"),
      drn_info_file = file.path("modules", module_uid, "PRISM", "drn_info.tsv"),
      analysis_status = "completed",
      status_reason = "",
      stringsAsFactors = FALSE
    ),
    drug_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(ml_top = ml_top_file, survival = survival_file, drug = drug_summary_file)
}

test_that("step 6C writes one self-contained selected module per subtype", {
  division <- tempfile("triage-division-")
  selection <- tempfile("triage-selection-")
  inputs <- tempfile("triage-inputs-")
  output <- tempfile("triage-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, inputs, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_triage_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  files <- write_triage_inputs(inputs, manifest$module_uid)

  selected <- triage_modules(
    module_manifest_file = attr(manifest, "manifest_file"),
    ml_top_file = files$ml_top,
    survival_file = files$survival,
    drug_response_summary_file = files$drug,
    output_dir = output,
    min_module_size = 10L
  )

  expect_equal(nrow(selected), 1L)
  expect_identical(selected$module_uid, manifest$module_uid)
  expect_identical(selected$selection_rank, 1L)
  expect_identical(selected$primary_drug_panel, "PRISM")
  expect_true(file.exists(attr(selected, "selected_file")))
  expect_true(file.exists(attr(selected, "evidence_file")))
  expect_true(file.exists(attr(selected, "stepwise_file")))
  expect_true(all(file.exists(file.path(
    output,
    c(selected$node_file, selected$edge_file, selected$drn_file, selected$drn_info_file)
  ))))
  stepwise <- utils::read.delim(attr(selected, "stepwise_file"), check.names = FALSE)
  expect_identical(stepwise$gate_name, c("ml_gate", "prognosis_gate", "module_size_gate", "drug_response_gate"))
  expect_true(all(stepwise$cumulative_pass))
})

test_that("step 6C strict mode reports a subtype with no passing candidate", {
  division <- tempfile("triage-division-")
  selection <- tempfile("triage-selection-")
  inputs <- tempfile("triage-inputs-")
  output <- tempfile("triage-output-")
  dir.create(division)
  on.exit(unlink(c(division, selection, inputs, output), recursive = TRUE, force = TRUE), add = TRUE)
  write_triage_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  files <- write_triage_inputs(inputs, manifest$module_uid, survival_direction = "High_score_better")

  expect_error(
    triage_modules(
      module_manifest_file = attr(manifest, "manifest_file"),
      ml_top_file = files$ml_top,
      survival_file = files$survival,
      drug_response_summary_file = files$drug,
      output_dir = output,
      min_module_size = 10L
    ),
    "No module passed all triage gates for subtype C1"
  )
  expect_false(dir.exists(output))
})
