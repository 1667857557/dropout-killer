test_that("end-to-end result satisfies selective invariants", {
  set.seed(11)
  x <- matrix(rexp(240), 20, 12); x[x < 0.7] <- 0
  rownames(x) <- paste0("g", 1:20); colnames(x) <- paste0("c", 1:12)
  z <- matrix(rnorm(36), 12, 3, dimnames = list(colnames(x), paste0("PC", 1:3)))
  fit <- dropout_killer(x, z, membership = rep(1, 12), rank = 2, max_rank = 4,
                        quantile_prob = 0.1, threshold = 0.5, min_cells = 4,
                        neighbor_k = 5)
  chk <- validate_dropout_result(fit, x)
  expect_true(chk$valid)
  expect_equal(dim(fit$expression), dim(x))
  expect_false(any(c("alpha", "ppi_prior", "pathway_prior", "prior_weights", "prior_standardize") %in% names(fit$settings)))
  expect_false(any(c("prior_prediction", "prior_available", "prior_support") %in% names(fit$events)))
})
