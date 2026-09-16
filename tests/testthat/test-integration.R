integration_lean_model <- function(threshold = 1e10) {
  set.seed(1001)
  f <- data.frame(z = rnorm(600), local_hierarchy = rnorm(600), H = rnorm(600))
  truth <- runif(600) < plogis(f$z + 0.5 * f$local_hierarchy + 0.25 * f$H)
  m <- fit_lean_detector(f, truth)
  m$threshold <- threshold
  m
}

test_that("P1 stabilized-state default satisfies selective invariants", {
  set.seed(11)
  x <- matrix(rpois(20 * 12, 2), 20, 12)
  x[, colSums(x) == 0] <- 1
  rownames(x) <- paste0("g", 1:20); colnames(x) <- paste0("c", 1:12)
  z <- matrix(rnorm(36), 12, 3, dimnames = list(colnames(x), paste0("PC", 1:3)))
  group <- rep(c("A", "B"), each = 6)
  fit <- dropout_killer(
    x, z, group = group, lean_model = integration_lean_model(),
    gamma = 3, k_knn = 3, rank = 2,
    factor_rank = 2, factor_features = 10,
    min_feature_observed = 2, min_target_observed = 2
  )
  chk <- validate_dropout_result(fit, x)
  expect_true(chk$valid)
  expect_equal(dim(fit$expression), dim(x))
  expect_equal(fit$settings$detection_method, "Supercell_hierarchy_Lean_membership")
  expect_equal(fit$settings$recovery_method, "p1_stabilized_state")
  expect_equal(fit$settings$factor_crossfit_folds, 5L)
  expect_true(fit$settings$support_adaptive_rank)
  expect_equal(fit$settings$bias_kappa, 10)
  expect_equal(fit$settings$predictor_smoothing, 0.25)
  expect_true(fit$uncertainty_available)
  expect_true(inherits(fit$predictive_variance, "sparseMatrix"))
})

test_that("neighbor comparison results do not encode unknown variance as zero", {
  set.seed(19)
  x <- matrix(rpois(20 * 12, 2), 20, 12)
  x[, colSums(x) == 0] <- 1
  rownames(x) <- paste0("g", 1:20); colnames(x) <- paste0("c", 1:12)
  z <- matrix(rnorm(36), 12, 3, dimnames = list(colnames(x), paste0("PC", 1:3)))
  group <- rep(c("A", "B"), each = 6)
  fit <- dropout_killer(
    x, z, group = group, lean_model = integration_lean_model(),
    gamma = 3, k_knn = 3, rank = 2,
    neighbor_k = 5, recovery_method = "neighbor"
  )
  expect_false(fit$uncertainty_available)
  expect_null(fit$predictive_variance)
  expect_error(sample_dropout_expression(fit), "unavailable")
})
