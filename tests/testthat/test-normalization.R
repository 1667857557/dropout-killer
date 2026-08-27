test_that("ALRA library+log normalization matches the reference formula", {
  x <- matrix(c(1, 3, 0, 2, 4, 0, 5, 5, 10, 0, 0, 10), nrow = 3)
  rownames(x) <- paste0("g", 1:3); colnames(x) <- paste0("c", 1:4)
  lib <- colSums(x)
  expected <- log1p(sweep(x, 2L, lib, FUN = "/") * 1e4)
  got <- DropoutKiller:::.dk_alra_library_log(x)
  expect_equal(dim(got), dim(expected))
  expect_equal(dimnames(got), dimnames(expected))
  expect_equal(as.numeric(got), as.numeric(expected), tolerance = 1e-12)
  expect_equal(got[x == 0], rep(0, sum(x == 0)))
  expect_equal(attr(got, "dropoutkiller_normalization")$scale_factor, 1e4)
  expect_equal(attr(got, "dropoutkiller_normalization")$method,
               "ALRA_library_size_log1p")
})

test_that("sparse ALRA normalization matches dense normalization", {
  x <- matrix(c(1, 3, 0, 2, 4, 0, 5, 5, 10, 0, 0, 10), nrow = 3)
  rownames(x) <- paste0("g", 1:3); colnames(x) <- paste0("c", 1:4)
  dense <- DropoutKiller:::.dk_alra_library_log(x)
  sparse <- DropoutKiller:::.dk_alra_library_log(Matrix::Matrix(x, sparse = TRUE))
  expect_true(inherits(sparse, "sparseMatrix"))
  expect_equal(dim(sparse), dim(dense))
  expect_equal(dimnames(sparse), dimnames(dense))
  expect_equal(as.numeric(sparse), as.numeric(dense), tolerance = 1e-12)
  expect_equal(attr(sparse, "dropoutkiller_normalization")$scale_factor, 1e4)
})

test_that("zero-library cells are rejected instead of silently breaking alignment", {
  x <- matrix(c(1, 2, 0, 0), nrow = 2)
  colnames(x) <- c("c1", "c2")
  expect_error(DropoutKiller:::.dk_alra_library_log(x), "positive library size")
})

test_that("high-level workflow normalizes once and validates against raw input", {
  x <- matrix(c(
    1, 0, 3, 2,
    0, 4, 1, 2,
    5, 2, 0, 1,
    2, 1, 3, 0,
    4, 2, 1, 3
  ), nrow = 5, byrow = TRUE)
  rownames(x) <- paste0("g", 1:5); colnames(x) <- paste0("c", 1:4)
  z <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  expected <- DropoutKiller:::.dk_alra_library_log(x)
  fit <- dropout_killer(
    x, z, membership = rep(1, 4), min_cells = 10,
    normalize = TRUE, normalization_scale_factor = 1e4
  )
  expect_equal(as.matrix(fit$expression), as.matrix(expected), tolerance = 1e-12)
  expect_true(fit$settings$normalize)
  expect_equal(fit$settings$normalization, "ALRA_library_size_log1p")
  expect_equal(fit$settings$normalization_scale_factor, 1e4)
  expect_true(validate_dropout_result(fit, x)$valid)
})

test_that("normalize FALSE preserves the supplied working scale", {
  x <- matrix(c(1, 0, 2, 3, 4, 1, 0, 2, 1, 3, 2, 4), nrow = 3)
  rownames(x) <- paste0("g", 1:3); colnames(x) <- paste0("c", 1:4)
  z <- matrix(seq_len(8), nrow = 4,
              dimnames = list(colnames(x), c("PC1", "PC2")))
  fit <- dropout_killer(x, z, membership = rep(1, 4), min_cells = 10, normalize = FALSE)
  expect_equal(as.matrix(fit$expression), x)
  expect_equal(fit$settings$normalization, "none")
})
