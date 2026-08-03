test_that("paper configuration covers all subtypes and protects selected-only downstream analysis", {
  path <- system.file("config", "paper_luad.yml", package = "MLSnpDR")
  expect_true(nzchar(path))
  config <- read_mlsnpdr_config(path)
  expect_identical(config$project$target_subtypes, c("C1", "C2", "C3", "C4"))
  expect_identical(config$ml$top_k_per_subtype, 10L)
  expect_identical(config$ml$expected_feature_count, 34L)
  expect_identical(
    config$pipeline$network_methods,
    c("String", "physicalPPIN", "chengF")
  )
  expect_identical(config$drug_response$analysis_scope, "all_manifest_modules")
  expect_identical(config$candidate_selection$select_n_per_subtype, 1L)
  expect_identical(config$annotation$databases, c("GO_BP", "KEGG", "HALLMARK"))
  expect_identical(config$annotation$top_n, 15L)
  expect_true(config$post_selection$selected_only)
})

test_that("selected-only cannot be disabled", {
  path <- system.file("config", "paper_luad.yml", package = "MLSnpDR")
  expect_true(nzchar(path))
  config <- yaml::read_yaml(path)
  config$post_selection$selected_only <- FALSE
  expect_error(validate_mlsnpdr_config(config), "selected_only must be true")
})

test_that("stage registry preserves subnetDR numbering with ML inserted after step 6", {
  stages <- mlsnpdr_stage_registry()
  expect_identical(
    stages$stage,
    c("01", "02", "03", "04", "05", "06", "06A", "06B", "06C", "07", "08", "09")
  )
  expect_identical(stages$name[stages$stage == "06A"], "module_features")
  expect_identical(stages$name[stages$stage == "06B"], "ml_scoring")
  expect_identical(stages$name[stages$stage == "06C"], "module_triage")
  expect_true(all(stages$scope[stages$stage %in% c("07", "08", "09")] ==
    "one_best_module_per_subtype"))
  expect_true(all(stages$implemented[stages$stage %in% c("04", "05")]))
  expect_true(stages$implemented[stages$stage == "06"])
  expect_true(stages$implemented[stages$stage == "06A"])
  expect_true(stages$implemented[stages$stage == "06B"])
  expect_true(stages$implemented[stages$stage == "06C"])
  expect_true(all(stages$implemented[stages$stage %in% c("07", "08", "09")]))
})

test_that("complete runner selects the step-4 through step-9 chain", {
  path <- system.file("config", "default.yml", package = "MLSnpDR")
  run <- run_ML_SnpDR(path, dry_run = TRUE)
  expect_identical(
    run$plan$name,
    c(
      "module_selection", "module_annotation", "drug_response",
      "module_features", "ml_scoring", "module_triage",
      "sequence_smiles", "binding_score", "perturbation_score"
    )
  )
  expect_identical(run$status, "planned")
})
