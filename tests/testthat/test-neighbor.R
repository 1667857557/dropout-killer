test_that("weighted borrowing uses local positive donors", {
  x <- matrix(c(0, 10, 2, 5, 5, 5), nrow = 2, byrow = TRUE)
  rownames(x) <- c("g1", "g2"); colnames(x) <- c("c1", "c2", "c3")
  z <- matrix(c(0, 0.1, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- weighted_neighbor_prediction(x, z, c(1, 1, 1), mask, k = 2, sigma = 1)
  expect_true(as.numeric(p[1, 1]) > 2)
  expect_true(as.numeric(p[1, 1]) < 10)
})

test_that("all-neighbor borrowing can include local zeros", {
  x <- matrix(c(0, 10, 0), nrow = 1, dimnames = list("g1", c("c1", "c2", "c3")))
  z <- matrix(c(0, 0.1, 0.2), ncol = 1, dimnames = list(colnames(x), "PC1"))
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p_pos <- weighted_neighbor_prediction(x, z, rep(1, 3), mask, k = 2, sigma = 1, positive_only = TRUE)
  p_all <- weighted_neighbor_prediction(x, z, rep(1, 3), mask, k = 2, sigma = 1, positive_only = FALSE)
  expect_gt(as.numeric(p_pos[1, 1]), as.numeric(p_all[1, 1]))
})
