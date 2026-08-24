test_that("Codex review regressions remain fixed", {
  direct <- names(formals(recover_dropout_expression))
  expect_equal(direct[5:10], c("neighbor_k", "neighbor_sigma", "min_positive_neighbors",
                               "neighbor_positive_only", "cap_quantile", "return_details"))
  workflow <- names(formals(dropout_killer))
  i <- match("min_negative", workflow)
  expect_equal(workflow[(i + 1L):(i + 7L)], c("neighbor_k", "neighbor_sigma",
                                              "min_positive_neighbors", "neighbor_positive_only",
                                              "cap_quantile", "seed", "return_score"))
  expect_gt(match("recovery_method", workflow), match("return_score", workflow))
})

test_that("sparse Matrix summary dispatch does not require Matrix attachment", {
  was_attached <- "package:Matrix" %in% search()
  if (was_attached) detach("package:Matrix", character.only = TRUE)
  on.exit(if (was_attached) suppressPackageStartupMessages(library(Matrix)), add = TRUE)

  x <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(2L, 1L), dims = c(2L, 2L))
  dimnames(x) <- list(c("g1", "g2"), c("c1", "c2"))
  empty_mask <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = dim(x)) > 0
  expect_silent(recover_dropout_expression(x, empty_mask, c(1, 1)))

  xn <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(2L, 1L), x = c(2, 3), dims = c(2L, 2L),
                             dimnames = dimnames(x))
  mask <- Matrix::sparseMatrix(i = 1L, j = 1L, x = TRUE, dims = dim(xn))
  d <- recover_dropout_expression(xn, mask, c(1, 1), min_target_observed = 10,
                                  return_details = TRUE)
  expect_equal(nrow(d$events), 1L)
  expect_equal(d$events$i, 1L)
  expect_equal(d$events$j, 1L)

  score <- Matrix::sparseMatrix(i = 1L, j = 2L, x = 0.99, dims = dim(xn), dimnames = dimnames(xn))
  selected <- select_dropout_mask(score, threshold = 0.95)
  expect_true(as.logical(selected[1, 2]))

  res <- list(expression = x, mask = empty_mask)
  class(res) <- "DropoutKillerResult"
  chk <- validate_dropout_result(res, x)
  expect_true(chk$valid)
})
