test_that("module UID normalizes legacy spelling", {
  expect_identical(
    make_module_uid("PhysicalPPIN", "louvain", "c3", "m10"),
    "physicalppin__louvain__C3__M10"
  )
  expect_identical(
    make_module_uid("String", "WF", "C1", "M2"),
    "string__wf__C1__M2"
  )
})

test_that("module UID round-trips", {
  parsed <- parse_module_uid("physicalppin__louvain__C3__M10")
  expect_identical(parsed$legacy_module_id, "physicalPPIN|Louvain|C3|M10")
  expect_identical(parsed$module_uid, "physicalppin__louvain__C3__M10")
})

test_that("unknown identifiers fail loudly", {
  expect_error(make_module_uid("unknown", "WF", "C1", "M1"), "Unknown network")
  expect_error(make_module_uid("String", "other", "C1", "M1"), "Unknown module method")
})

