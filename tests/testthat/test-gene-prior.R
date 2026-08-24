test_that("PPI/pathway prior aligns and predicts masked targets", {
  x <- matrix(c(0, 1, 2, 2, 3, 4, 4, 3, 2), nrow = 3, byrow = TRUE)
  rownames(x) <- paste0("g", 1:3); colnames(x) <- paste0("c", 1:3)
  A <- matrix(0, 3, 3, dimnames = list(rownames(x), rownames(x)))
  A[1, 2:3] <- 1
  mask <- Matrix::sparseMatrix(i = 1, j = 1, x = TRUE, dims = dim(x))
  p <- gene_prior_prediction(x, A, rep(1, 3), mask)
  expect_true(is.finite(as.numeric(p[1, 1])))
  expect_gte(as.numeric(p[1, 1]), 0)
})

test_that("PPI and pathway priors retain expression gene dimensions", {
  genes <- paste0("g", 1:4)
  ppi <- matrix(0, 4, 4, dimnames = list(genes, genes))
  ppi[1, 2] <- ppi[2, 1] <- 1
  pathway <- matrix(0, 4, 4, dimnames = list(genes, genes))
  pathway[3, 4] <- 2
  W <- combine_gene_prior(genes = genes, ppi = ppi, pathway = pathway,
                          weights = c(ppi = 0.6, pathway = 0.4))
  expect_equal(dim(W), c(4, 4))
  expect_true(all(diag(as.matrix(W)) == 0))
  expect_gt(as.numeric(W[1, 2]), 0)
  expect_gt(as.numeric(W[3, 4]), 0)
})

test_that("PPI edge-list prior is symmetric by default", {
  e <- data.frame(a = c("g1", "g2"), b = c("g2", "g3"), w = c(2, 1))
  A <- gene_prior_from_edges(e, genes = paste0("g", 1:3), source = "a", target = "b",
                             weight = "w", prior_type = "ppi")
  expect_equal(as.numeric(A[2, 1]), 2)
  expect_equal(as.numeric(A[1, 2]), 2)
  expect_equal(as.numeric(A[3, 2]), 1)
  expect_equal(as.numeric(A[2, 3]), 1)
})

test_that("pathway edge-list prior can preserve direction", {
  e <- data.frame(source = c("g1", "g2"), target = c("g3", "g3"), w = c(2, 1))
  A <- gene_prior_from_edges(e, genes = paste0("g", 1:3), weight = "w",
                             prior_type = "pathway")
  expect_equal(as.numeric(A[3, 1]), 2)
  expect_equal(as.numeric(A[3, 2]), 1)
  expect_equal(as.numeric(A[1, 3]), 0)
})

test_that("GRN prior type is not accepted", {
  e <- data.frame(source = "g1", target = "g2")
  expect_error(gene_prior_from_edges(e, genes = c("g1", "g2"), prior_type = "grn"))
})
