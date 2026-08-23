test_that("selective recovery preserves every observed non-zero", {
  x <- matrix(c(0, 8, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.2, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  y <- recover_dropout_expression(x, mask, c(1, 1, 1), z, alpha = 1, neighbor_k = 2)
  expect_equal(y[x != 0], x[x != 0])
  expect_gt(y[1, 1], 0)
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
  y <- recover_dropout_expression(x, mask, c(1, 1, 1), z, alpha = 1, neighbor_k = 2)
  expect_true(inherits(y, "sparseMatrix"))
  expect_equal(as.numeric(y[x != 0]), as.numeric(x[x != 0]))
  expect_gt(as.numeric(y[1, 1]), 0)
})
