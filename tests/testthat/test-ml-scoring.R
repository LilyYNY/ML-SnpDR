write_ml_division_fixture <- function(root) {
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

ml_feature_fixture <- function() {
  result <- data.frame(
    module_id = "String|Louvain|C1|M1",
    module_label = "M1_C1",
    Subtype = "C1",
    Network = "String",
    Method = "Louvain",
    module_size = 10,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (index in 2:34) result[[sprintf("feature_%02d", index)]] <- index
  result
}

make_ml_feature_output <- function(division, selection, feature_source, feature_output) {
  write_ml_division_fixture(division)
  manifest <- module_selection(
    base_input_path = division,
    base_output_path = selection,
    network_method = "String",
    module_method = "Louvain",
    numberCutoff = 9,
    subtypes = "C1"
  )
  utils::write.table(
    ml_feature_fixture(),
    feature_source,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  features <- prepare_module_features(
    module_manifest_file = attr(manifest, "manifest_file"),
    feature_file = feature_source,
    output_dir = feature_output
  )
  list(manifest = manifest, features = features)
}

test_that("step 6B validates all-module probabilities and writes Top-K", {
  division <- tempfile("ml-division-")
  selection <- tempfile("ml-selection-")
  feature_source <- tempfile(fileext = ".tsv")
  feature_output <- tempfile("ml-features-")
  score_file <- tempfile(fileext = ".tsv")
  output <- tempfile("ml-output-")
  dir.create(division)
  on.exit(
    unlink(c(division, selection, feature_source, feature_output, score_file, output), recursive = TRUE, force = TRUE),
    add = TRUE
  )
  built <- make_ml_feature_output(division, selection, feature_source, feature_output)
  scores <- data.frame(
    module_id = "String|Louvain|C1|M1",
    score_type = "nested_5x5fold_oof_probability_after_inner_hyperopt",
    selected_model = "GradientBoosting",
    prob_C1 = 0.8,
    prob_C2 = 0.1,
    prob_C3 = 0.05,
    prob_C4 = 0.05,
    stringsAsFactors = FALSE
  )
  utils::write.table(scores, score_file, sep = "\t", quote = FALSE, row.names = FALSE)

  result <- prepare_ml_scores(
    module_features_file = attr(built$features, "feature_file"),
    score_file = score_file,
    output_dir = output,
    top_k = 1L
  )

  expect_equal(nrow(result), 1L)
  expect_identical(result$module_uid, built$manifest$module_uid)
  expect_identical(result$predicted_subtype, "C1")
  expect_equal(result$target_subtype_probability, 0.8)
  expect_equal(result$probability_margin, 0.7)
  expect_identical(result$rank_in_subtype, 1L)
  expect_identical(result$score_type, "nested_oof_probability")
  top <- utils::read.delim(attr(result, "top_file"), check.names = FALSE)
  expect_equal(nrow(top), 1L)
  expect_identical(top$target_subtype, "C1")
  expect_true(top$ml_gate)
})

test_that("step 6B strict mode requires normalized probabilities", {
  division <- tempfile("ml-division-")
  selection <- tempfile("ml-selection-")
  feature_source <- tempfile(fileext = ".tsv")
  feature_output <- tempfile("ml-features-")
  score_file <- tempfile(fileext = ".tsv")
  output <- tempfile("ml-output-")
  dir.create(division)
  on.exit(
    unlink(c(division, selection, feature_source, feature_output, score_file, output), recursive = TRUE, force = TRUE),
    add = TRUE
  )
  built <- make_ml_feature_output(division, selection, feature_source, feature_output)
  scores <- data.frame(
    module_uid = built$manifest$module_uid,
    prob_C1 = 0.8,
    prob_C2 = 0.2,
    prob_C3 = 0.2,
    prob_C4 = 0.1,
    stringsAsFactors = FALSE
  )
  utils::write.table(scores, score_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_ml_scores(
      module_features_file = attr(built$features, "feature_file"),
      score_file = score_file,
      output_dir = output,
      top_k = 1L
    ),
    "sum to 1"
  )
})

test_that("step 6B reports insufficient eligible modules", {
  division <- tempfile("ml-division-")
  selection <- tempfile("ml-selection-")
  feature_source <- tempfile(fileext = ".tsv")
  feature_output <- tempfile("ml-features-")
  score_file <- tempfile(fileext = ".tsv")
  output <- tempfile("ml-output-")
  dir.create(division)
  on.exit(
    unlink(c(division, selection, feature_source, feature_output, score_file, output), recursive = TRUE, force = TRUE),
    add = TRUE
  )
  built <- make_ml_feature_output(division, selection, feature_source, feature_output)
  scores <- data.frame(
    module_uid = built$manifest$module_uid,
    prob_C1 = 0.1,
    prob_C2 = 0.7,
    prob_C3 = 0.1,
    prob_C4 = 0.1,
    stringsAsFactors = FALSE
  )
  utils::write.table(scores, score_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_ml_scores(
      module_features_file = attr(built$features, "feature_file"),
      score_file = score_file,
      output_dir = output,
      top_k = 1L,
      require_predicted_subtype_match = TRUE
    ),
    "only 0 eligible"
  )
})
