test_that("Codex review regressions remain fixed", {
  direct <- names(formals(recover_dropout_expression))
  expect_equal(direct[5:10], c("neighbor_k", "neighbor_sigma", "min_positive_neighbors",
                               "neighbor_positive_only", "cap_quantile", "return_details"))
  workflow <- names(formals(dropout_killer))
  i <- match("min_negative", workflow)
  expect_equal(workflow[(i + 1L):(i + 7L)], c("neighbor_k", "neighbor_sigma",
                                              "min_positive_neighbors", "neighbor_positive_only",
                                              "cap_quantile", "seed", "return_score"))
  expect_gt(match("recovery_method", workflow), match("return_score", workflow))
})
