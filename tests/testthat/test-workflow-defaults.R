test_that("high-level detector defaults advertise scGACL first", {
  f <- formals(dropout_killer)
  expect_equal(
    eval(f$detection_method),
    c("scgacl_gamma_normal", "alra_global", "alra_global_by_group",
      "eb_zero_null", "alra_quantile")
  )
  expect_equal(eval(f$scgacl_dropout_threshold), 0.5)
})

test_that("legacy public positional slots before detection controls are unchanged", {
  f <- names(formals(dropout_killer))
  i <- match("min_negative", f)
  expect_equal(
    f[(i + 1L):(i + 7L)],
    c("neighbor_k", "neighbor_sigma", "min_positive_neighbors",
      "neighbor_positive_only", "cap_quantile", "seed", "return_score")
  )
  expect_gt(match("detection_method", f), match("min_target_observed", f))
  expect_gt(match("scgacl_dropout_threshold", f), match("detection_method", f))
})
