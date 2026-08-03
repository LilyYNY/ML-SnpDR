test_that("paper configuration protects selected-only downstream analysis", {
  path <- system.file("config", "paper_luad.yml", package = "MLSnpDR")
  expect_true(nzchar(path))
  config <- read_mlsnpdr_config(path)
  expect_identical(config$project$target_subtypes, "C3")
  expect_identical(config$ml$top_k_per_subtype, 10L)
  expect_true(config$post_selection$selected_only)
})

test_that("selected-only cannot be disabled", {
  path <- system.file("config", "paper_luad.yml", package = "MLSnpDR")
  expect_true(nzchar(path))
  config <- yaml::read_yaml(path)
  config$post_selection$selected_only <- FALSE
  expect_error(validate_mlsnpdr_config(config), "selected_only must be true")
})
