test_that("selective recovery equals weighted-neighbor prediction and preserves observed non-zero values", {
  x <- matrix(c(0, 8, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- weighted_neighbor_prediction(x, z, c(1, 1, 1), mask, k = 2)
  d <- recover_dropout_expression(x, mask, c(1, 1, 1), z, neighbor_k = 2, return_details = TRUE)
  expect_equal(d$expression[x != 0], x[x != 0])
  expect_equal(d$expression[1, 1], as.numeric(p[1, 1]))
  expect_gt(d$expression[1, 1], 0)
  expect_false(any(c("prior_prediction", "prior_available", "prior_support") %in% names(d$events)))
})

test_that("recovery refuses a mask over non-zero data", {
  x <- matrix(1, 2, 2); z <- matrix(1:4, 2, 2)
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  expect_error(recover_dropout_expression(x, mask, c(1, 1), z), "non-zero")
})

test_that("sparse input remains sparse and nonzero observations remain exact", {
  x <- Matrix::Matrix(matrix(c(0, 8, 2, 5, 5, 5), nrow = 2, byrow = TRUE), sparse = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  y <- recover_dropout_expression(x, mask, c(1, 1, 1), z, neighbor_k = 2)
  expect_true(inherits(y, "sparseMatrix"))
  expect_equal(as.numeric(y[x != 0]), as.numeric(x[x != 0]))
  expect_gt(as.numeric(y[1, 1]), 0)
})
