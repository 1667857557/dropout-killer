test_that("gene prior aligns and predicts masked targets", {
  x <- matrix(c(0, 1, 2, 2, 3, 4, 4, 3, 2), nrow = 3, byrow = TRUE)
  rownames(x) <- paste0("g", 1:3); colnames(x) <- paste0("c", 1:3)
  A <- matrix(0, 3, 3, dimnames = list(rownames(x), rownames(x)))
  A[1, 2:3] <- 1
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- gene_prior_prediction(x, A, rep(1, 3), mask)
  expect_true(is.finite(as.numeric(p[1, 1])))
  expect_gte(as.numeric(p[1, 1]), 0)
})

test_that("combined networks retain expression gene dimensions", {
  genes <- paste0("g", 1:4)
  a <- diag(4); dimnames(a) <- list(genes, genes)
  b <- a[c(1, 3), c(2, 4), drop = FALSE]
  rownames(b) <- genes[c(1, 3)]; colnames(b) <- genes[c(2, 4)]
  W <- combine_gene_prior(list(a, b), genes)
  expect_equal(dim(W), c(4, 4))
  expect_true(all(diag(as.matrix(W)) == 0))
})

test_that("edge-list network uses target-by-source orientation", {
  e <- data.frame(tf = c("g1", "g2"), target = c("g3", "g3"), w = c(2, 1))
  A <- gene_network_from_edges(e, genes = paste0("g", 1:3), source = "tf", target = "target", weight = "w")
  expect_equal(as.numeric(A[3, 1]), 2)
  expect_equal(as.numeric(A[3, 2]), 1)
  expect_equal(as.numeric(A[1, 3]), 0)
})
