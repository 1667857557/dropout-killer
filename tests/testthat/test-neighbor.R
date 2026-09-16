test_that("neighbor comparison uses local positive donors", {
  x <- matrix(c(0, 10, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.1, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- recover_dropout_expression(
    x, mask, c(1, 1, 1), z,
    neighbor_k = 2, neighbor_sigma = 1,
    recovery_method = "neighbor"
  )
  expect_true(as.numeric(p[1, 1]) > 2)
  expect_true(as.numeric(p[1, 1]) < 10)
})

test_that("neighbor comparison can include local zeros when requested", {
  x <- matrix(c(0, 10, 0), nrow = 1,
              dimnames = list("g1", c("c1", "c2", "c3")))
  z <- matrix(c(0, 0.1, 0.2), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p_pos <- recover_dropout_expression(
    x, mask, rep(1, 3), z, neighbor_k = 2, neighbor_sigma = 1,
    neighbor_positive_only = TRUE, recovery_method = "neighbor"
  )
  p_all <- recover_dropout_expression(
    x, mask, rep(1, 3), z, neighbor_k = 2, neighbor_sigma = 1,
    neighbor_positive_only = FALSE, recovery_method = "neighbor"
  )
  expect_gt(as.numeric(p_pos[1, 1]), as.numeric(p_all[1, 1]))
})
