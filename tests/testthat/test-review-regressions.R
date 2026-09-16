test_that("deprecated detector and compatibility APIs are not exported", {
  exports <- getNamespaceExports("DropoutKiller")
  deprecated <- c(
    "DropoutKiller", "run_dropout_killer", "dropout_membership",
    "local_alra_detect", "local_alra_score", "select_dropout_mask",
    "masked_factor_prediction", "weighted_neighbor_prediction",
    "recovery_architecture_prediction"
  )
  expect_false(any(deprecated %in% exports))
  expect_equal(
    eval(formals(dropout_killer)$detection_method),
    c("Supercell_hierarchy_Lean_membership",
      "Supercell_hierarchy_Lean_membership_WNN")
  )
  old_args <- c("membership", "split_by", "max_rank", "min_negative",
                "variance_prior_df", "alra_K", "alra_noise_start",
                "alra_choose_q", "alra_svd_q")
  expect_false(any(old_args %in% names(formals(dropout_killer))))
})

test_that("sparse Matrix paths do not require Matrix attachment", {
  was_attached <- "package:Matrix" %in% search()
  if (was_attached) detach("package:Matrix", character.only = TRUE)
  on.exit(if (was_attached) suppressPackageStartupMessages(library(Matrix)), add = TRUE)

  x <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(2L, 1L), dims = c(2L, 2L))
  dimnames(x) <- list(c("g1", "g2"), c("c1", "c2"))
  empty_mask <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = dim(x)) > 0
  z <- matrix(c(0, 1), ncol = 1, dimnames = list(colnames(x), "PC1"))
  expect_silent(recover_dropout_expression(x, empty_mask, c(1, 1), z))

  xn <- Matrix::sparseMatrix(
    i = c(1L, 2L), j = c(2L, 1L), x = c(2, 3), dims = c(2L, 2L),
    dimnames = dimnames(x)
  )
  mask <- Matrix::sparseMatrix(i = 1L, j = 1L, x = TRUE, dims = dim(xn))
  d <- recover_dropout_expression(
    xn, mask, c(1, 1), recovery_method = "masked_factor",
    min_target_observed = 10, return_details = TRUE
  )
  expect_equal(nrow(d$events), 1L)
  expect_equal(d$events$i, 1L)
  expect_equal(d$events$j, 1L)

  res <- list(
    expression = x, mask = empty_mask,
    settings = list(normalize = FALSE)
  )
  class(res) <- "DropoutKillerResult"
  chk <- validate_dropout_result(res, x)
  expect_true(chk$valid)
})
