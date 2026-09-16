test_that("default direct recovery uses P1 stabilized state and preserves observations", {
  x <- matrix(c(0, 0, 4, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:4)
  z <- matrix(c(0, 0.2, 0.8, 1), ncol = 1,
              dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  d <- recover_dropout_expression(
    x, mask, rep(1, 4), z,
    factor_rank = 1, factor_features = 2,
    min_feature_observed = 2, min_target_observed = 10,
    return_details = TRUE
  )
  expect_equal(d$expression[x != 0], x[x != 0])
  expect_gt(d$expression[1, 1], 0)
  expect_equal(d$events$recovery_method[1], "positive_membership_mean")
  expect_equal(d$events$target_mode[1], "positive")
  expect_true(is.finite(d$events$prediction_sd[1]))
  expect_true(d$uncertainty_available)
})

test_that("neighbor comparison engine remains available through expert recovery API", {
  x <- matrix(c(0, 8, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  d <- recover_dropout_expression(
    x, mask, c(1, 1, 1), z,
    neighbor_k = 2, return_details = TRUE,
    recovery_method = "neighbor"
  )
  expect_gt(d$expression[1, 1], 2)
  expect_lt(d$expression[1, 1], 8)
  expect_equal(d$events$recovery_method[1], "neighbor")
  expect_false(d$uncertainty_available)
})

test_that("production recovery defaults match the current P1 contract", {
  f <- formals(recover_dropout_expression)
  expect_equal(eval(f$recovery_method)[1], "p1_stabilized_state")
  expect_equal(eval(f$factor_ridge), 2)
  expect_equal(eval(f$min_target_observed), 8L)
  expect_equal(eval(f$factor_crossfit_folds), 5L)
  expect_true(eval(f$support_adaptive_rank))
  expect_equal(eval(f$bias_kappa), 10)
})

test_that("recovery refuses a mask over non-zero data", {
  x <- matrix(1, 2, 2)
  colnames(x) <- c("c1", "c2")
  z <- matrix(c(0, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  expect_error(recover_dropout_expression(x, mask, c(1, 1), z), "non-zero")
})

test_that("sparse input remains sparse and observed nonzeros remain exact", {
  x <- Matrix::Matrix(
    matrix(c(0, 0, 4, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE),
    sparse = TRUE
  )
  rownames(x) <- c("g1", "g2"); colnames(x) <- paste0("c", 1:4)
  z <- matrix(c(0, 0.2, 0.8, 1), ncol = 1,
              dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  y <- recover_dropout_expression(
    x, mask, rep(1, 4), z,
    factor_rank = 1, factor_features = 2,
    min_feature_observed = 2, min_target_observed = 10
  )
  expect_true(inherits(y, "sparseMatrix"))
  expect_equal(as.numeric(y[x != 0]), as.numeric(x[x != 0]))
  expect_gt(as.numeric(y[1, 1]), 0)
})
