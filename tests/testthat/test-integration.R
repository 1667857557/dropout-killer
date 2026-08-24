test_that("end-to-end result satisfies selective and uncertainty invariants", {
  set.seed(11)
  x <- matrix(rexp(240), 20, 12); x[x < 0.7] <- 0
  rownames(x) <- paste0("g", 1:20); colnames(x) <- paste0("c", 1:12)
  z <- matrix(rnorm(36), 12, 3, dimnames = list(colnames(x), paste0("PC", 1:3)))
  fit <- dropout_killer(x, z, membership = rep(1, 12), rank = 2, max_rank = 4,
                        quantile_prob = 0.1, threshold = 0.5, min_cells = 4,
                        factor_rank = 2, factor_features = 10,
                        min_feature_observed = 4, min_target_observed = 4)
  chk <- validate_dropout_result(fit, x)
  expect_true(chk$valid)
  expect_equal(dim(fit$expression), dim(x))
  expect_equal(fit$settings$recovery_method, "masked_factor")
  expect_true(inherits(fit$predictive_variance, "sparseMatrix"))
  expect_true(all(c("prediction_sd", "predictability", "recovery_method") %in% names(fit$events)))
  expect_false(any(c("alpha", "ppi_prior", "pathway_prior", "prior_weights", "prior_standardize") %in% names(fit$settings)))
})
