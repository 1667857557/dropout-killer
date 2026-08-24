test_that("default recovery uses masked factor and preserves observed non-zero values", {
  x <- matrix(c(0, 0, 4, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:4)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  d <- recover_dropout_expression(x, mask, rep(1, 4), factor_rank = 1,
                                  factor_features = 2, min_feature_observed = 3,
                                  min_target_observed = 10, return_details = TRUE)
  expect_equal(d$expression[x != 0], x[x != 0])
  expect_gt(d$expression[1, 1], 0)
  expect_equal(d$events$recovery_method[1], "membership_mean")
  expect_true(is.finite(d$events$prediction_sd[1]))
  expect_true(d$uncertainty_available)
})

test_that("legacy direct-recovery positional arguments retain their slots", {
  f <- names(formals(recover_dropout_expression))
  expect_equal(f[5:10], c("neighbor_k", "neighbor_sigma", "min_positive_neighbors",
                          "neighbor_positive_only", "cap_quantile", "return_details"))
  x <- matrix(c(0, 8, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- weighted_neighbor_prediction(x, z, c(1, 1, 1), mask, k = 2)
  d <- recover_dropout_expression(x, mask, c(1, 1, 1), z,
                                  2, NULL, 1L, TRUE, NULL, TRUE,
                                  recovery_method = "neighbor")
  expect_equal(d$expression[1, 1], as.numeric(p[1, 1]))
  expect_equal(d$events$recovery_method[1], "neighbor")
  expect_false(d$uncertainty_available)
})

test_that("main workflow preserves the pre-0.4 positional recovery layout", {
  f <- names(formals(dropout_killer))
  i <- match("min_negative", f)
  expect_equal(f[(i + 1L):(i + 7L)], c("neighbor_k", "neighbor_sigma",
                                       "min_positive_neighbors", "neighbor_positive_only",
                                       "cap_quantile", "seed", "return_score"))
  expect_gt(match("recovery_method", f), match("return_score", f))
})

test_that("recovery refuses a mask over non-zero data", {
  x <- matrix(1, 2, 2)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  expect_error(recover_dropout_expression(x, mask, c(1, 1)), "non-zero")
})

test_that("sparse input remains sparse and observed nonzeros remain exact", {
  x <- Matrix::Matrix(matrix(c(0, 0, 4, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE), sparse = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:4)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  y <- recover_dropout_expression(x, mask, rep(1, 4), factor_rank = 1,
                                  factor_features = 2, min_feature_observed = 3,
                                  min_target_observed = 10)
  expect_true(inherits(y, "sparseMatrix"))
  expect_equal(as.numeric(y[x != 0]), as.numeric(x[x != 0]))
  expect_gt(as.numeric(y[1, 1]), 0)
})
